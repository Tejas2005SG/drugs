import { v4 as uuidv4 } from 'uuid';
import GeneratenewMolecule from '../models/generatenew.model.js';
import CostEstimation from '../models/costestimination.model.js';
import DrugName from '../models/drugName.js';
import ResearchPaper from '../models/researchPapers.model.js';
import { SavedSearch } from '../models/savedSearch.model.js';
import Toxicity from '../models/toxicity.model.js';
import GeneratedResearchPaper from '../models/GeneratedResearchPaper.js';
import SummarySavedItem from '../models/summarySavedItems.js';

// Valid types for validation
const VALID_TYPES = [
  'newDrug',
  'costEstimation',
  'drugName',
  'researchPaper',
  'generatedResearchPaper',
  'targetPrediction',
  'toxicityResult',
];

// Get summary data (no auto-saving)
export const getSummary = async (req, res) => {
  try {
    const userId = req.user?._id;
    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    // Fetch all relevant data for the user
    const [
      newDrugs,
      costEstimations,
      drugNames,
      researchPapers,
      targetPredictions,
      toxicityResults,
      GeneratedResearchPaperDetails,
    ] = await Promise.all([
      GeneratenewMolecule.find({ userId }).sort({ createdAt: -1 }).lean(),
      CostEstimation.find({ userId }).sort({ createdAt: -1 }).lean(),
      DrugName.find({ userId }).sort({ createdAt: -1 }).lean(),
      ResearchPaper.find({ userId }).sort({ createdAt: -1 }).lean(),
      SavedSearch.find({ userId }).sort({ createdAt: -1 }).lean(),
      Toxicity.find({ userId }).sort({ createdAt: -1 }).lean(),
      GeneratedResearchPaper.find({ userId }).sort({ createdAt: -1 }).lean(),
    ]);

    res.status(200).json({
      success: true,
      data: {
        newDrugs,
        costEstimations,
        drugNames,
        researchPapers,
        targetPredictions,
        toxicityResults,
        GeneratedResearchPaperDetails,
      },
    });
  } catch (error) {
    console.error('Error fetching summary:', error);
    res.status(500).json({
      success: false,
      message: 'Failed to fetch summary data',
      error: error.message,
    });
  }
};

// Save a single summary item
export const saveItem = async (req, res) => {
  try {
    const { item, type, moleculeId } = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    if (!item || !type) {
      return res.status(400).json({ success: false, message: 'Item and type are required' });
    }

    if (!VALID_TYPES.includes(type)) {
      return res.status(400).json({ success: false, message: `Invalid type: ${type}` });
    }

    const savedItem = await SummarySavedItem.create({
      userId,
      type,
      moleculeId: moleculeId || uuidv4(), // Generate a new moleculeId if not provided
      data: item,
    });

    res.status(201).json({ success: true, data: savedItem });
  } catch (error) {
    console.error('Error saving item:', error);
    res.status(500).json({ success: false, message: 'Failed to save item', error: error.message });
  }
};

// Get saved items (grouped by moleculeId)
export const getSavedItems = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    const savedItems = await SummarySavedItem.find({ userId })
      .sort({ createdAt: -1 })
      .lean();

    // Group items by moleculeId
    const groupedItems = savedItems.reduce((acc, item) => {
      const key = item.moleculeId || 'ungrouped';
      if (!acc[key]) {
        acc[key] = [];
      }
      acc[key].push(item);
      return acc;
    }, {});

    res.status(200).json({ success: true, data: groupedItems });
  } catch (error) {
    console.error('Error fetching saved items:', error);
    res.status(500).json({ success: false, message: 'Failed to fetch saved items', error: error.message });
  }
};

// Bulk save summary items
export const bulkSaveItems = async (req, res) => {
  try {
    const items = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    if (!Array.isArray(items) || items.length === 0) {
      return res.status(400).json({ success: false, message: 'Items array is required and cannot be empty' });
    }

    // Validate all items
    for (const { item, type, moleculeId } of items) {
      if (!item || !type) {
        return res.status(400).json({ success: false, message: 'Each item must have item and type' });
      }
      if (!VALID_TYPES.includes(type)) {
        return res.status(400).json({ success: false, message: `Invalid type: ${type}` });
      }
    }

    // Prepare items for bulk insert
    const itemsToSave = items.map(({ item, type, moleculeId }) => ({
      userId,
      type,
      moleculeId: moleculeId || uuidv4(),
      data: item,
    }));

    const savedItems = await SummarySavedItem.insertMany(itemsToSave);

    res.status(201).json({ success: true, data: savedItems });
  } catch (error) {
    console.error('Error bulk saving items:', error);
    res.status(500).json({ success: false, message: 'Failed to bulk save items', error: error.message });
  }
};

// Get new molecule generation progress
export const getNewMoleculeProgress = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    const latestMolecule = await GeneratenewMolecule.findOne({ userId })
      .sort({ createdAt: -1 })
      .lean();

    const progress = {
      status: latestMolecule ? 'Completed' : 'Idle',
      percentage: latestMolecule ? 100 : 0,
    };

    res.status(200).json({
      success: true,
      data: {
        progress,
        data: latestMolecule || null,
      },
    });
  } catch (error) {
    console.error('Error fetching new molecule progress:', error);
    res.status(500).json({
      success: false,
      message: 'Failed to fetch new molecule progress',
      error: error.message,
    });
  }
};