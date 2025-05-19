import GeneratenewMolecule from '../models/generatenew.model.js';
import CostEstimation from '../models/costestimination.model.js';
import DrugName from '../models/drugName.js';
import ResearchPaper from '../models/researchPapers.model.js';
import SavedSearch from '../models/savedSearch.model.js';
import Toxicity from '../models/toxicity.model.js';

export const getSummary = async (req, res) => {
  try {
    const userId = req.user?._id;
    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    // Fetch all relevant data for the user
    const [
      newDrugs,
      costEstimations,
      drugNames,
      researchPapers,
      targetPredictions,
      toxicityResults
    ] = await Promise.all([
      GeneratenewMolecule.find({ userId }).sort({ created: -1 }),
      CostEstimation.find({ userId }).sort({ created: -1 }),
      DrugName.find({ userId }).sort({ createdAt: -1 }),
      ResearchPaper.find({ userId }).sort({ createdAt: -1 }),
      SavedSearch.find({ userId }).sort({ createdAt: -1 }),
      Toxicity.find({ userId }).sort({ created: -1 })
    ]);

    res.status(200).json({
      success: true,
      data: {
        newDrugs,
        costEstimations,
        drugNames,
        researchPapers,
        targetPredictions,
        toxicityResults
      }
    });
  } catch (error) {
    console.error('Error fetching summary:', error);
    res.status(500).json({
      success: false,
      message: 'Failed to fetch summary data',
      error: error.message
    });
  }
};