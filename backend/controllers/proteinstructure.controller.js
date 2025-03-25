import ProteinStructure from '../models/protienstructure.model.js'; // Fixed typo in import
import { User } from '../models/auth.model.js'; // Added User model import
import GeneratenewMolecule from '../models/generatenew.model.js'; // Added GeneratenewMolecule model import
import axios from 'axios';
import dotenv from 'dotenv';
import OpenAI from 'openai';

dotenv.config();

const QWQ_API_KEY = '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';
const MOLMIM_API_KEY = process.env.MOLMIM_API_KEY;
const MOLMIM_API_URL = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate";
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

const openai = new OpenAI({
  apiKey: QWQ_API_KEY,
  baseURL: 'https://integrate.api.nvidia.com/v1',
});
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

      const prompt = `
     Analyze the drug molecule represented by the SMILES string "${smiles}" and provide a detailed report that includes the following sections:

### 1. **Structural Analysis**
- **Core Structure Identification**: Describe the main structural features of the molecule, including any fused ring systems, heterocycles, and functional groups.
- **Stereochemistry**: Identify and explain the stereochemical configurations at all chiral centers (e.g., [C@@], [C@H]) and their relevance to biological activity.
- **Substituents**: List and characterize significant substituents attached to the core structure, explaining their potential impact on the molecule's properties.

### 2. **Chemical Properties**
- **Molecular Weight**: Calculate and provide the exact molecular weight of the compound.
- **Physicochemical Properties**: Include values for:
  - LogP (partition coefficient)
  - Polar Surface Area (PSA)
  - Hydrogen Bond Donors/Acceptors
  - Rotatable Bonds
  - pKa values for ionizable groups
- **Solubility Profile**: Predict solubility in various solvents (aqueous, organic) based on functional groups.
- **Melting Point and Boiling Point**: Provide estimated values based on structural analogs.

### 3. **ADMET Profile**
- **Absorption**: Discuss bioavailability and factors affecting absorption.
- **Distribution**: Predict volume of distribution and blood-brain barrier permeability.
- **Metabolism**: Identify metabolic pathways and potential metabolites, including enzyme interactions (e.g., CYP450).
- **Excretion**: Discuss elimination routes and half-life.
- **Toxicity Predictions**: Highlight any known or predicted toxicological concerns, including hERG inhibition or cytotoxicity.

### 4. **Biological Activity**
- **Target Interaction**: Identify potential biological targets (receptors, enzymes) and predict binding affinities. Include:
  - Mechanism of action for known targets.
  - Any relevant SAR (Structure-Activity Relationship) data.
- **Therapeutic Applications**: Discuss potential uses in therapy based on structural similarity to known drugs or biological activity.

### 5. **Synthesis Pathways**
- Outline possible synthetic routes to obtain the compound, referencing established methods in literature.
- Discuss any challenges or considerations in synthesis.

### 6. **Clinical Context**
- Provide a summary of any clinical data available regarding this compound, including:
  - Approved indications (if applicable).
  - Clinical trial results or ongoing studies.
  - Comparison with similar drugs in terms of efficacy and safety profiles.

### Output Requirements
- Format the response using clear subheadings for each section.
- Include tables or bullet points where appropriate for clarity.
- Cite relevant databases (PubChem, ChEMBL, DrugBank) for additional context or data sources.

### Conclusion
Summarize the overall potential of this molecule as a therapeutic agent, highlighting any critical research questions that remain unanswered or areas for further investigation.



      `;
      
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

export const generatenewmolecule = async (req, res) => {
  try {
    const { id } = req.params;
    const { smilesoffirst, smilesofsecond, newmoleculetitle } = req.body;

    console.log(`Processing request for user ID: ${id}`);
    console.log("Request body:", { smilesoffirst, smilesofsecond, newmoleculetitle });

    // Validation
    if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
      return res.status(400).json({ message: "SMILES strings and molecule title are required" });
    }
    if (!id) {
      return res.status(400).json({ message: "User ID is required" });
    }

    const prompt = `
You are an expert cheminformatics AI designed to combine two SMILES strings into a new molecule. Follow these steps strictly and ensure your response ends with a valid JSON object. The fields in your response will be extracted and saved to a MongoDB database, so they must be clearly labeled and structured.

1. Analyze the Input SMILES Strings:
   - Parse SMILES 1: "${smilesoffirst}" and describe its chemical structure (e.g., functional groups, atoms, bonds).
   - Parse SMILES 2: "${smilesofsecond}" and describe its chemical structure similarly.
   - Identify key reactive sites or functional groups in each molecule that could be used for bonding.
2. Combine into a New Molecule:
   - Propose a chemically reasonable bond to link the two molecules (e.g., an ether linkage -O-, an ester -COO-, a carbon-carbon bond).
   - Explain the rationale for choosing this bond type (e.g., based on functional group compatibility, stability).
   - Detail the step-by-step process of how the atoms are connected, including which atoms are involved and any adjustments needed.
3. Generate the New Molecule:
   - Provide the SMILES string for the resulting molecule.
   - Determine and provide the IUPAC name for the new molecule.
4. Therapeutic Potential:
   - Analyze the new molecule’s structure to hypothesize potential biological activity (e.g., similarity to known drug classes).
   - Suggest one or more diseases this molecule might target, with a brief explanation linking chemical features to mechanisms.
   - Note this is speculative, not experimental.
5. Return the Result:
   - At the end of your response, provide a valid JSON object with the following structure. This JSON object will be parsed and its fields (newSmiles, newIupacName, conversionDetails, potentialDiseases) will be saved to a MongoDB database.
   - The JSON object must be the last part of your response, with no additional text, comments, or whitespace after it:
     {
       "newSmiles": "[NEW_SMILES]",
       "newIupacName": "[NEW_IUPAC_NAME]",
       "conversionDetails": "[DETAILED_EXPLANATION]",
       "potentialDiseases": "[LIST_OF_DISEASES_WITH_EXPLANATIONS]"
     }
   - Ensure the JSON object is properly formatted, contains all required fields, and has valid values (e.g., no empty strings or null values).
   - If you cannot generate a valid SMILES string or IUPAC name, provide a placeholder value (e.g., "Unable to generate SMILES" or "Unknown IUPAC Name") and explain why in the conversionDetails.
   - If you cannot determine potential diseases, provide a placeholder value (e.g., "Unknown potential diseases") and explain why in the potentialDiseases field.
`;

    console.log("Sending streaming request to NVIDIA NIM API...");
    res.setHeader('Content-Type', 'text/event-stream');
    res.setHeader('Cache-Control', 'no-cache');
    res.setHeader('Connection', 'keep-alive');

    const completion = await openai.chat.completions.create({
      model: "qwen/qwq-32b", // Consider switching to "gpt-3.5-turbo" if the issue persists
      messages: [{ role: "user", content: prompt }],
      temperature: 0.6,
      top_p: 0.7,
      max_tokens: 4096,
      stream: true,
    });

    let fullResponse = '';
    for await (const chunk of completion) {
      const content = chunk.choices[0]?.delta?.content || '';
      fullResponse += content;
      res.write(`data: ${content}\n\n`);
    }
    res.write('data: [DONE]\n\n');
    res.end();

    // Log the full response for debugging
    console.log("Full QWQ response:", fullResponse);

    // Extract JSON-like object from the streaming response
    let newSmiles = "";
    let newIupacName = "";
    let conversionDetails = "";
    let potentialDiseases = "";

    // Try to extract JSON object
    const jsonMatch = fullResponse.match(/{[\s\S]*?}$/); // Match the last JSON object in the response
    if (jsonMatch) {
      try {
        const result = JSON.parse(jsonMatch[0]);
        newSmiles = result.newSmiles || "";
        newIupacName = result.newIupacName || "";
        conversionDetails = result.conversionDetails || "";
        potentialDiseases = result.potentialDiseases || "";
        console.log("Extracted fields from JSON:", { newSmiles, newIupacName, conversionDetails, potentialDiseases });
      } catch (parseError) {
        console.error("Failed to parse JSON from response:", parseError);
      }
    }

    // Fallback: Parse the fields from the text if JSON is not found or parsing fails
    if (!jsonMatch || !newSmiles) {
      console.log("Falling back to text parsing...");
      const lines = fullResponse.split('\n');

      for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        console.log(`Parsing line ${i}:`, line);

        // Try different formats for field extraction
        if (line.toLowerCase().includes("newsmiles")) {
          const match = line.match(/(?:newSmiles|newsmiles)\s*[:=]\s*(.+)/i);
          if (match) {
            newSmiles = match[1].trim();
            console.log("Found newSmiles:", newSmiles);
          }
        }
        if (line.toLowerCase().includes("newiupacname")) {
          const match = line.match(/(?:newIupacName|newiupacname)\s*[:=]\s*(.+)/i);
          if (match) {
            newIupacName = match[1].trim();
            console.log("Found newIupacName:", newIupacName);
          }
        }
        if (line.toLowerCase().includes("conversiondetails")) {
          const match = line.match(/(?:conversionDetails|conversiondetails)\s*[:=]\s*(.+)/i);
          if (match) {
            conversionDetails = match[1].trim();
            console.log("Found conversionDetails:", conversionDetails);
          }
        }
        if (line.toLowerCase().includes("potentialdiseases")) {
          const match = line.match(/(?:potentialDiseases|potentialdiseases)\s*[:=]\s*(.+)/i);
          if (match) {
            potentialDiseases = match[1].trim();
            console.log("Found potentialDiseases:", potentialDiseases);
          }
        }
      }
    }

    // Validate and set fallback values if fields are empty
    newSmiles = newSmiles || `${smilesoffirst}.${smilesofsecond}`;
    newIupacName = newIupacName || "Unknown IUPAC Name";
    conversionDetails = conversionDetails || "Unable to parse conversion details from the response.";
    potentialDiseases = potentialDiseases || "Unknown potential diseases due to lack of information in the response.";

    console.log("Final extracted fields:", { newSmiles, newIupacName, conversionDetails, potentialDiseases });

    // Save the new molecule to the database
    console.log("Saving new molecule to database with raw streaming response...");
    const newMolecule = new GeneratenewMolecule({
      smilesoffirst,
      smilesofsecond,
      newmoleculetitle,
      newSmiles,
      newIupacName,
      conversionDetails,
      potentialDiseases,
      information: fullResponse, // Save the raw streaming response
      userId: id,
    });
    await newMolecule.save();

    console.log("Updating user with new molecule ID:", newMolecule._id);
    await User.findByIdAndUpdate(id, { $push: { proteinStructures: newMolecule._id } }, { new: true });
  } catch (error) {
    console.error("Error generating new molecule:", error.message, error.stack);
    res.write(`data: {"error": "${error.message}"}\n\n`);
    res.end();
  }
};

export const getgeneratednewmolecule = async (req, res) => {
  try {
    const userId = req.user.id; // From protectRoute middleware

    console.log(`Request received for /api/protein/generatednewmolecule for user: ${userId}`);
    const molecules = await GeneratenewMolecule.find({ userId }).sort({ created: -1 });

    if (!molecules || molecules.length === 0) {
      return res.status(404).json({ message: "No generated molecules found for this user" });
    }

    const moleculeData = molecules.map((molecule) => ({
      id: molecule._id,
      newSmiles: molecule.newSmiles,
      newIupacName: molecule.newIupacName,
      newmoleculetitle: molecule.newmoleculetitle,
      conversionDetails: molecule.conversionDetails,
      potentialDiseases: molecule.potentialDiseases,
      information: molecule.information, // Raw streaming response as text
      created: molecule.created,
    }));

    res.status(200).json({
      molecules: moleculeData,
      total: moleculeData.length,
    });
  } catch (error) {
    console.error("Error retrieving generated molecules:", error);
    res.status(500).json({ message: "Failed to retrieve molecules", error: error.message });
  }
};;


// export const generatenewmolecule = async (req, res) => {   
//   try {
//     const { id } = req.params;
//     const { smilesoffirst, smilesofsecond, newmoleculetitle } = req.body;

//     console.log(`Processing request for user ID: ${id}`);
//     console.log("Request body:", { smilesoffirst, smilesofsecond, newmoleculetitle });

//     // Validation
//     if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
//       return res.status(400).json({ message: "SMILES strings and molecule title are required" });
//     }

//     if (!id) {
//       return res.status(400).json({ message: "User ID is required" });
//     }

//     if (!GEMINI_API_KEY) {
//       console.error("GEMINI_API_KEY is not set in environment variables");
//       return res.status(500).json({ message: "Gemini API key is not configured" });
//     }

//     const prompt = `
//       You are an expert cheminformatics AI designed to combine two SMILES strings into a new molecule. Follow these steps:
// 1. Analyze the two input SMILES strings to understand their chemical structures.
// 2. Combine them into a new molecule by forming a chemically reasonable bond (e.g., an ether linkage -O-, an ester, or a carbon-carbon bond) between appropriate atoms.
// 3. Generate the SMILES string for the new molecule and its corresponding IUPAC name.
// 4. Return the result in this JSON format:
// {
//   "newSmiles": "[NEW_SMILES]",
//   "newIupacName": "[NEW_IUPAC_NAME]"
// }

// Example:
// - Input SMILES 1: "CCO" (ethanol)
// - Input SMILES 2: "c1ccccc1" (benzene)
// - Combined: Link the oxygen of ethanol to a carbon of benzene to form an ether.
// - Output: {
//   "newSmiles": "CCOc1ccccc1",
//   "newIupacName": "Ethoxybenzene"
// }

// Input:
// - SMILES 1: "${smilesoffirst}"
// - SMILES 2: "${smilesofsecond}"
//     `;

//     console.log("Sending request to Gemini API...");
//     const geminiResponse = await axios.post(
//       `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
//       { contents: [{ parts: [{ text: prompt }] }] },
//       { headers: { "Content-Type": "application/json" } }
//     );

//     const geminiContent = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
//     console.log("Gemini API response:", geminiContent);
//     if (!geminiContent) {
//       throw new Error("No content returned from Gemini API");
//     }

//     let result;
//     try {
//       result = JSON.parse(geminiContent);
//     } catch (parseError) {
//       console.error("Failed to parse Gemini response as JSON:", geminiContent);
//       const smilesMatch = geminiContent.match(/"newSmiles":\s*"([^"]+)"/);
//       const iupacMatch = geminiContent.match(/"newIupacName":\s*"([^"]+)"/);
//       if (smilesMatch && iupacMatch) {
//         result = { newSmiles: smilesMatch[1], newIupacName: iupacMatch[1] };
//       } else {
//         result = {
//           newSmiles: `${smilesoffirst}${smilesofsecond}`,
//           newIupacName: "Generated Molecule (Fallback Name)",
//         };
//         console.warn("Using fallback molecule due to invalid Gemini response");
//       }
//     }

//     const { newSmiles, newIupacName } = result;
//     if (!newSmiles || !newIupacName) {
//       throw new Error("Gemini response missing required fields");
//     }

//     const information = JSON.stringify({
//       newSmiles,
//       newIupacName,
//       description: "Generated by combining two SMILES strings",
//     });

//     console.log("Saving new molecule to database...");
//     const newMolecule = new GeneratenewMolecule({
//       smilesoffirst,
//       smilesofsecond,
//       newmoleculetitle, // Add the new field
//       information,
//       userId: id,
//     });
//     await newMolecule.save();

//     console.log("Updating user with new molecule ID:", newMolecule._id);
//     const userUpdate = await User.findByIdAndUpdate(
//       id,
//       { $push: { proteinStructures: newMolecule._id } }, // Corrected typo: "protienStructures" -> "proteinStructures"
//       { new: true }
//     );
//     if (!userUpdate) {
//       throw new Error(`User with ID ${id} not found`);
//     }

//     res.status(201).json({
//       newSmiles,
//       newIupacName,
//       newmoleculetitle, // Include in response
//       id: newMolecule._id,
//     });
//   } catch (error) {
//     console.error("Error generating new molecule:", error.message, error.stack);
//     res.status(500).json({ message: "Failed to generate new molecule", error: error.message });
//   }
// };


// export const generatenewmolecule = async (req, res) => {
//   try {
//     const { id } = req.params;
//     const { smilesoffirst, smilesofsecond, newmoleculetitle } = req.body;

//     console.log(`Processing request for user ID: ${id}`);
//     console.log("Request body:", { smilesoffirst, smilesofsecond, newmoleculetitle });

//     // Validation
//     if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
//       return res.status(400).json({ message: "SMILES strings and molecule title are required" });
//     }

//     if (!id) {
//       return res.status(400).json({ message: "User ID is required" });
//     }

//     const prompt = `
//      You are an expert cheminformatics AI designed to combine two SMILES strings into a new molecule. Follow these steps:

// 1. **Analyze the Input SMILES Strings**:
//    - Parse SMILES 1: "${smilesoffirst}" and describe its chemical structure (e.g., functional groups, atoms, bonds).
//    - Parse SMILES 2: "${smilesofsecond}" and describe its chemical structure similarly.
//    - Identify key reactive sites or functional groups in each molecule that could be used for bonding.

// 2. **Combine into a New Molecule**:
//    - Propose a chemically reasonable bond to link the two molecules (e.g., an ether linkage -O-, an ester -COO-, a carbon-carbon bond, or another suitable connection).
//    - Explain the rationale for choosing this bond type (e.g., based on functional group compatibility, stability, or common chemical synthesis patterns).
//    - Detail the step-by-step process of how the atoms are connected, including which atoms from each SMILES string are involved and any adjustments (e.g., hydrogen removal) needed for a valid structure.

// 3. **Generate the New Molecule**:
//    - Provide the SMILES string for the resulting molecule.
//    - Determine and provide the IUPAC name for the new molecule.

// 4. **Therapeutic Potential**:
//    - Analyze the new molecule’s structure and functional groups to hypothesize potential biological activity (e.g., similarity to known drug classes like analgesics, antibiotics, or anticancer agents).
//    - Suggest one or more diseases this molecule might target based on its structure (e.g., inflammation, bacterial infections, cancer), with a brief explanation linking the chemical features to potential mechanisms of action.
//    - Note that this is a speculative analysis based on cheminformatics reasoning, not experimental data.

// 5. **Return the Result**:
//    - Format the output as a JSON object with the following structure:
//      {
//        "newSmiles": "[NEW_SMILES]",
//        "newIupacName": "[NEW_IUPAC_NAME]",
//        "conversionDetails": "[DETAILED_STEP_BY_STEP_EXPLANATION]",
//        "potentialDiseases": "[LIST_OF_DISEASES_WITH_BRIEF_EXPLANATIONS]"
//      }

// Example:
// - Input SMILES 1: "CCO" (ethanol)
// - Input SMILES 2: "c1ccccc1" (benzene)
// - Process:
//   - SMILES 1: "CCO" represents ethanol (CH3-CH2-OH), with a hydroxyl group (-OH) as a reactive site.
//   - SMILES 2: "c1ccccc1" represents benzene (C6H6), an aromatic ring with carbons available for substitution.
//   - Bond choice: Form an ether linkage by connecting the oxygen of ethanol to a carbon on the benzene ring, removing one hydrogen from benzene and the hydrogen from the -OH group.
//   - Resulting structure: CH3-CH2-O-C6H5.
// - Output:
//   {
//     "newSmiles": "CCOc1ccccc1",
//     "newIupacName": "Ethoxybenzene",
//     "conversionDetails": "Ethanol (CCO) has a hydroxyl group (-OH) on the second carbon. Benzene (c1ccccc1) is an aromatic ring. The oxygen from ethanol's -OH replaces a hydrogen on one of benzene's carbons, forming an ether linkage (-O-). The hydrogen from the -OH is removed to balance the valency, yielding CCOc1ccccc1.",
//     "potentialDiseases": "Anxiety (Ether compounds like ethoxybenzene may have sedative properties similar to some anxiolytics); Inflammation (Aromatic ethers can sometimes act as mild anti-inflammatory agents by interacting with specific enzymes)."
//   }

// Input:
// - SMILES 1: "${smilesoffirst}"
// - SMILES 2: "${smilesofsecond}"
//     `;

//     console.log("Sending request to NVIDIA NIM API...");
//     const completion = await openai.chat.completions.create({
//       model: "qwen/qwq-32b",
//       messages: [{ role: "user", content: prompt }],
//       temperature: 0.6,
//       top_p: 0.7,
//       max_tokens: 4096,
//     });

//     const qwqResponse = completion.choices[0]?.message?.content;
//     console.log("NVIDIA NIM API response:", qwqResponse);
//     if (!qwqResponse) {
//       throw new Error("No content returned from NVIDIA NIM API");
//     }

//     let result;
//     try {
//       result = JSON.parse(qwqResponse);
//     } catch (parseError) {
//       console.error("Failed to parse QWQ response as JSON:", qwqResponse);
//       const smilesMatch = qwqResponse.match(/"newSmiles":\s*"([^"]+)"/);
//       const iupacMatch = qwqResponse.match(/"newIupacName":\s*"([^"]+)"/);
//       if (smilesMatch && iupacMatch) {
//         result = { newSmiles: smilesMatch[1], newIupacName: iupacMatch[1] };
//       } else {
//         result = {
//           newSmiles: `${smilesoffirst}${smilesofsecond}`,
//           newIupacName: "Generated Molecule (Fallback Name)",
//         };
//         console.warn("Using fallback molecule due to invalid QWQ response");
//       }
//     }

//     const { newSmiles, newIupacName } = result;
//     if (!newSmiles || !newIupacName) {
//       throw new Error("QWQ response missing required fields");
//     }

//     const information = JSON.stringify({
//       newSmiles,
//       newIupacName,
//       description: "Generated by combining two SMILES strings",
//     });

//     console.log("Saving new molecule to database...");
//     const newMolecule = new GeneratenewMolecule({
//       smilesoffirst,
//       smilesofsecond,
//       newmoleculetitle,
//       information,
//       userId: id,
//     });
//     await newMolecule.save();

//     console.log("Updating user with new molecule ID:", newMolecule._id);
//     const userUpdate = await User.findByIdAndUpdate(
//       id,
//       { $push: { proteinStructures: newMolecule._id } },
//       { new: true }
//     );
//     if (!userUpdate) {
//       throw new Error(`User with ID ${id} not found`);
//     }

//     res.status(201).json({
//       newSmiles,
//       newIupacName,
//       newmoleculetitle,
//       id: newMolecule._id,
//     });
//   } catch (error) {
//     console.error("Error generating new molecule:", error.message, error.stack);
//     res.status(500).json({ message: "Failed to generate new molecule", error: error.message });
//   }
// };


// export const getgeneratednewmolecule = async (req, res) => {
//   try {
//     const userId = req.user.id; // From protectRoute middleware

//     console.log(`Request received for /api/protein/generatednewmolecule for user: ${userId}`);

//     // Fetch all molecules for the user
//     const molecules = await GeneratenewMolecule.find({ userId }).sort({ created: -1 }); // Sort by creation date, descending

//     if (!molecules || molecules.length === 0) {
//       return res.status(404).json({ message: "No generated molecules found for this user" });
//     }

//     // Parse the information field for each molecule and include newmoleculetitle
//     const moleculeData = molecules.map((molecule) => {
//       const info = JSON.parse(molecule.information);
//       return {
//         id: molecule._id,
//         newSmiles: info.newSmiles,
//         newIupacName: info.newIupacName,
//         newmoleculetitle: molecule.newmoleculetitle, // Add the newmoleculetitle field
//         created: molecule.created,
//       };
//     });

//     res.status(200).json({
//       molecules: moleculeData,
//       total: moleculeData.length,
//     });
//   } catch (error) {
//     console.error("Error retrieving generated molecules:", error);
//     res.status(500).json({ message: "Failed to retrieve molecules", error: error.message });
//   }
// };





// New endpoint to fetch the latest molecule for a user
// export const getLatestGeneratedMolecule = async (req, res) => {
//   try {
//     const { userId } = req.params; // User ID from route
//     const authUserId = req.user.id; // From protectRoute middleware

//     console.log(`Request received for /api/protein/generatednewmolecule/latest/${userId}`);

//     if (userId !== authUserId) {
//       return res.status(403).json({ message: "Unauthorized to access this user’s molecules" });
//     }

//     const latestMolecule = await GeneratenewMolecule
//       .findOne({ userId })
//       .sort({ createdAt: -1 }); // Sort by creation date, descending

//     if (!latestMolecule) {
//       return res.status(404).json({ message: "No generated molecules found for this user" });
//     }

//     const info = JSON.parse(latestMolecule.information);

//     res.status(200).json({
//       newSmiles: info.newSmiles,
//       newIupacName: info.newIupacName,
//       id: latestMolecule._id,
//     });
//   } catch (error) {
//     console.error("Error retrieving latest generated molecule:", error);
//     res.status(500).json({ message: "Failed to retrieve latest molecule", error: error.message });
//   }
// };