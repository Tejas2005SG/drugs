import Toxicity from '../models/toxicity.model.js';

// Mock ProTox-II Prediction (replace with real API if available)
const predictToxicity = async (smiles) => {
  // Simulated response based on ProTox-II output structure
  return {
    smiles,
    acuteToxicity: {
      LD50: `${Math.floor(Math.random() * 5000 + 100)} mg/kg`, // Mock LD50
      toxicityClass: `Class ${Math.floor(Math.random() * 5) + 1}`, // Mock class (1-5)
    },
    endpoints: {
      hepatotoxicity: Math.random() > 0.5 ? 'Active' : 'Inactive',
      carcinogenicity: Math.random() > 0.5 ? 'Active' : 'Inactive',
    },
  };
};

// Predict Toxicity
export const predictToxicityController = async (req, res) => {
  try {
    const { smiles } = req.body;
    const userId = req.user._id; // Assuming auth middleware sets req.user

    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }

    // Get toxicity prediction
    const toxicityResult = await predictToxicity(smiles);
    
    // Save prediction to database
    const toxicityEntry = new Toxicity({ 
      smiles, 
      toxicityResult, 
      userId 
    });
    
    await toxicityEntry.save();

    // Return successful response
    res.status(200).json({
      message: 'Toxicity predicted successfully',
      result: toxicityResult,
    });
  } catch (error) {
    console.error('Error predicting toxicity:', error);
    res.status(500).json({ 
      message: 'Failed to predict toxicity', 
      error: error.message 
    });
  }
};

// Get Toxicity History
export const getToxicityHistory = async (req, res) => {
  try {
    const userId = req.user._id; // Assuming auth middleware sets req.user
    
    // Find all toxicity entries for the current user, sorted by creation date
    const history = await Toxicity.find({ userId }).sort({ created: -1 });
    
    res.status(200).json({ history });
  } catch (error) {
    console.error('Error fetching toxicity history:', error);
    res.status(500).json({ 
      message: 'Failed to fetch history', 
      error: error.message 
    });
  }
};