// controllers/toxicity.controllers.js
import Toxicity from '../models/toxicity.model.js';

const predictToxicity = async (smiles) => {
  return {
    smiles,
    acuteToxicity: {
      LD50: `${Math.floor(Math.random() * 5000 + 100)} mg/kg`,
      toxicityClass: `Class ${Math.floor(Math.random() * 5) + 1}`,
    },
    endpoints: {
      hepatotoxicity: Math.random() > 0.5 ? 'Active' : 'Inactive',
      carcinogenicity: Math.random() > 0.5 ? 'Active' : 'Inactive',
    },
  };
};

export const predictToxicityController = async (req, res) => {
  try {
    const { smiles } = req.body;
    const userId = req.user?._id;

    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }
    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const toxicityResult = await predictToxicity(smiles);

    const toxicityEntry = new Toxicity({
      smiles,
      toxicityResult,
      userId,
    });

    await toxicityEntry.save();
    console.log('Saved toxicity entry:', toxicityEntry);

    res.status(200).json({
      message: 'Toxicity predicted successfully',
      result: toxicityResult,
    });
  } catch (error) {
    console.error('Error predicting toxicity:', {
      message: error.message,
      stack: error.stack,
      body: req.body,
    });
    res.status(500).json({
      message: 'Failed to predict toxicity',
      error: error.message,
    });
  }
};

export const getToxicityHistory = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const history = await Toxicity.find({ userId }).sort({ created: -1 });
    console.log('Fetched history for user:', userId, 'Count:', history.length);

    res.status(200).json({ history });
  } catch (error) {
    console.error('Error fetching toxicity history:', {
      message: error.message,
      stack: error.stack,
    });
    res.status(500).json({
      message: 'Failed to fetch history',
      error: error.message,
    });
  }
};

export const saveGeminiAnalysis = async (req, res) => {
  try {
    const { smiles, geminiAnalysis } = req.body;
    const userId = req.user?._id;

    if (!smiles || !geminiAnalysis) {
      return res.status(400).json({ message: 'SMILES and Gemini analysis are required' });
    }

    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const toxicityEntry = await Toxicity.findOne({ smiles, userId }).sort({ created: -1 });

    if (!toxicityEntry) {
      return res.status(404).json({ message: 'No matching toxicity prediction found' });
    }

    toxicityEntry.geminiAnalysis = geminiAnalysis;
    await toxicityEntry.save();
    console.log('Updated toxicity entry with Gemini analysis:', toxicityEntry);

    res.status(200).json({ message: 'Gemini analysis saved successfully' });
  } catch (error) {
    console.error('Error saving Gemini analysis:', {
      message: error.message,
      stack: error.stack,
      body: req.body,
    });
    res.status(500).json({ message: 'Failed to save Gemini analysis', error: error.message });
  }
};