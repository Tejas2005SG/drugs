import ProteinStructure from '../models/protienstructure.model.js';
import axios from 'axios';
import dotenv from 'dotenv';

dotenv.config();

const MOLMIM_API_KEY = process.env.MOLMIM_API_KEY;
const MOLMIM_API_URL = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate";
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

export const getProteinStructure = async (req, res) => {
  try {
    const { id } = req.query;
    
    if (id) {
      const proteinStructure = await ProteinStructure.findById(id);
      if (!proteinStructure) {
        return res.status(404).json({ message: 'Protein structure not found' });
      }
      return res.status(200).json(proteinStructure);
    } else {
      const proteinStructures = await ProteinStructure.find()
        .sort({ created: -1 })
        .limit(20);
      return res.status(200).json(proteinStructures);
    }
  } catch (error) {
    console.error('Error fetching protein structures:', error);
    return res.status(500).json({ message: 'Server error', error: error.message });
  }
};

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
        contents: [
          {
            parts: [
              {
                text: prompt
              }
            ]
          }
        ]
      },
      {
        headers: {
          'Content-Type': 'application/json',
        }
      }
    );

    if (geminiResponse.data.candidates && geminiResponse.data.candidates.length > 0) {
      const info = geminiResponse.data.candidates[0].content.parts[0].text;
      return res.status(200).json({ information: info });
    } else {
      return res.status(404).json({ message: 'No content returned from Gemini API' });
    }
  } catch (error) {
    console.error('Error fetching variant information:', error);
    
    // Provide more detailed error information for debugging
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

export const postProteinStructure = async (req, res) => {
  try {
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

    // Validate algorithm
    const validAlgorithms = ["CMA-ES", "SSD"];
    if (!validAlgorithms.includes(algorithm)) {
      return res.status(400).json({ message: `Invalid algorithm: ${algorithm}. Supported algorithms are ${validAlgorithms.join(', ')}.` });
    }

    let molecules;
    if (algorithm === "SSD") {
      // Mock SSD implementation (simplified example)
      console.log('Using Sampling Standard Deviation (mock implementation)');
      molecules = Array.from({ length: numMolecules }, (_, i) => ({
        sample: smiles, // Use the input SMILES as a placeholder
        score: Math.random(), // Mock QED score
        similarity: Math.random() * (1 - minSimilarity) + minSimilarity, // Mock similarity
      }));
    } else {
      // Log the request payload for MolMIM API (without showing full API key)
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
      console.log('Using MolMIM API Key:', MOLMIM_API_KEY ? 'Key is configured' : 'Not set');

      // Make API call to NVIDIA MolMIM
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

      console.log('MolMIM API Response Status:', response.status);

      // Parse molecules
      if (typeof response.data.molecules === 'string') {
        molecules = JSON.parse(response.data.molecules);
      } else if (Array.isArray(response.data.molecules)) {
        molecules = response.data.molecules;
      } else {
        throw new Error('Invalid API response: "molecules" is neither a string nor an array');
      }

      if (!Array.isArray(molecules)) {
        throw new Error('Parsed "molecules" is not an array');
      }
    }

    // Call Gemini API to get detailed information about the molecule
    let moleculeInfo = '';
    try {
      if (!GEMINI_API_KEY) {
        throw new Error('Gemini API key is not set in environment variables');
      }

      const prompt = `Provide detailed information about the molecule with the SMILES string "${smiles}". Include its chemical properties, potential uses, and any relevant biological or medicinal information.`;
      
      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        {
          contents: [
            {
              parts: [
                {
                  text: prompt
                }
              ]
            }
          ]
        },
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );

      console.log('Gemini API Response Status:', geminiResponse.status);

      // Extract the generated text from the Gemini API response
      if (geminiResponse.data.candidates && geminiResponse.data.candidates.length > 0) {
        moleculeInfo = geminiResponse.data.candidates[0].content.parts[0].text;
      } else {
        throw new Error('No content returned from Gemini API');
      }
    } catch (geminiError) {
      console.error('Error calling Gemini API:', geminiError.message);
      moleculeInfo = 'Failed to retrieve detailed information about the molecule.';
    }

    // Create a new protein structure in the database
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
      information: moleculeInfo // Save the Gemini API response
    });
    
    await newProteinStructure.save();
    
    return res.status(201).json({
      message: 'Protein structure generated successfully',
      structure: newProteinStructure
    });
  } catch (error) {
    console.error('Error generating protein structure:', error);
    
    if (error.response) {
      console.error('API Error Status:', error.response.status);
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