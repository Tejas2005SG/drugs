import ProteinStructure from '../models/protienstructure.model.js'; // Fixed typo in import
import { User } from '../models/auth.model.js'; // Added User model import
import axios from 'axios';
import dotenv from 'dotenv';

dotenv.config();

const MOLMIM_API_KEY = process.env.MOLMIM_API_KEY;
const MOLMIM_API_URL = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate";
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

export const getProteinStructure = async (req, res) => {
  try {
    const { id } = req.params; // Changed from req.query to req.params
    
    if (!id) {
      return res.status(400).json({ message: 'User ID is required' });
    }

    // Find user and populate their protein structures
    const user = await User.findById(id).populate('protienStructures');
    if (!user) {
      return res.status(404).json({ message: 'User not found' });
    }

    const proteinStructures = user.protienStructures
      .sort((a, b) => b.created - a.created) // Sort by creation date descending
      .slice(0, 20); // Limit to 20

    return res.status(200).json(proteinStructures);
  } catch (error) {
    console.error('Error fetching protein structures:', error);
    return res.status(500).json({ message: 'Server error', error: error.message });
  }
};

export const postProteinStructure = async (req, res) => {
  try {
    const { id } = req.params; // Changed from req.body to req.params
    const { 
      name, 
      smiles, 
      algorithm = "CMA-ES", 
      numMolecules = 30, 
      propertyName = "QED", 
      minimize = false, 
      minSimilarity = 0.3, 
      particles = 30, 
      iterations = 10 
    } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }

    if (!id) {
      return res.status(400).json({ message: 'User ID is required' });
    }

    // Validate algorithm
    const validAlgorithms = ["CMA-ES", "SSD"];
    if (!validAlgorithms.includes(algorithm)) {
      return res.status(400).json({ message: `Invalid algorithm: ${algorithm}. Supported algorithms are ${validAlgorithms.join(', ')}.` });
    }

    let molecules;
    if (algorithm === "SSD") {
      console.log('Using Sampling Standard Deviation (mock implementation)');
      molecules = Array.from({ length: numMolecules }, (_, i) => ({
        sample: smiles,
        score: Math.random(),
        similarity: Math.random() * (1 - minSimilarity) + minSimilarity,
      }));
    } else {
      const requestPayload = {
        algorithm,
        num_molecules: numMolecules,
        property_name: propertyName,
        minimize,
        min_similarity: minSimilarity,
        particles,
        iterations,
        smi: smiles
      };
      console.log('Sending request to MolMIM API:', JSON.stringify(requestPayload, null, 2));

      const response = await axios.post(
        MOLMIM_API_URL,
        requestPayload,
        {
          headers: {
            'Authorization': `Bearer ${MOLMIM_API_KEY}`,
            'Accept': 'application/json',
            'Content-Type': 'application/json',
          }
        }
      );

      molecules = typeof response.data.molecules === 'string' 
        ? JSON.parse(response.data.molecules) 
        : response.data.molecules;

      if (!Array.isArray(molecules)) {
        throw new Error('Parsed "molecules" is not an array');
      }
    }

    let moleculeInfo = '';
    try {
      if (!GEMINI_API_KEY) {
        throw new Error('Gemini API key is not set');
      }

      const prompt = `Provide detailed information about the molecule with the SMILES string "${smiles}". Include its chemical properties, potential uses, and any relevant biological or medicinal information.`;
      
      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        {
          contents: [{ parts: [{ text: prompt }] }]
        },
        {
          headers: { 'Content-Type': 'application/json' }
        }
      );

      moleculeInfo = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text || 
        'No content returned from Gemini API';
    } catch (geminiError) {
      console.error('Error calling Gemini API:', geminiError.message);
      moleculeInfo = 'Failed to retrieve detailed information about the molecule.';
    }

    // Create new protein structure
    const newProteinStructure = new ProteinStructure({
      name: name || 'Untitled Structure',
      smiles,
      properties: { 
        algorithm,
        propertyName,
        minimize,
        minSimilarity 
      },
      generatedStructures: molecules.map(mol => ({
        smiles: mol.sample,
        properties: {
          qed: mol.score,
          logp: mol.logp || null
        },
        similarity: mol.similarity || 0
      })),
      information: moleculeInfo,
      userId: id // Link to user
    });
    
    await newProteinStructure.save();

    // Update user's proteinStructures array
    await User.findByIdAndUpdate(
      id,
      { $push: { protienStructures: newProteinStructure._id } },
      { new: true }
    );

    return res.status(201).json({
      message: 'Protein structure generated successfully',
      structure: newProteinStructure
    });
  } catch (error) {
    console.error('Error generating protein structure:', error);
    
    if (error.response) {
      return res.status(error.response.status).json({ 
        message: 'MolMIM API error',
        error: error.message,
        status: error.response.status
      });
    }
    
    return res.status(500).json({ 
      message: 'Server error', 
      error: error.message 
    });
  }
};

// getVariantInfo remains unchanged as it doesn't interact with user data
export const getVariantInfo = async (req, res) => {
  try {
    const { smiles } = req.query;
    
    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }

    if (!GEMINI_API_KEY) {
      return res.status(500).json({ message: 'Gemini API key is not configured on the server' });
    }

    const prompt = `Provide detailed information about the molecule with the SMILES string "${smiles}". Include its chemical properties, potential uses, and any relevant biological or medicinal information.`;
    
    const geminiResponse = await axios.post(
      `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
      {
        contents: [{ parts: [{ text: prompt }] }]
      },
      {
        headers: { 'Content-Type': 'application/json' }
      }
    );

    if (geminiResponse.data.candidates?.length > 0) {
      const info = geminiResponse.data.candidates[0].content.parts[0].text;
      return res.status(200).json({ information: info });
    }
    
    return res.status(404).json({ message: 'No content returned from Gemini API' });
  } catch (error) {
    console.error('Error fetching variant information:', error);
    let errorMessage = 'Server error';
    let errorDetails = error.message;
    
    if (error.response) {
      errorMessage = `API Error (${error.response.status})`;
      errorDetails = JSON.stringify(error.response.data);
    }
    
    return res.status(500).json({ 
      message: errorMessage, 
      error: errorDetails 
    });
  }
};