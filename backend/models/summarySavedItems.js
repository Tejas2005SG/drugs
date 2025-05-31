import mongoose from 'mongoose';

const summarySavedItemSchema = new mongoose.Schema({
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    required: true,
    ref: 'User',
  },
  type: {
    type: String,
    required: true,
    enum: [
      'newDrug',
      'costEstimation',
      'drugName',
      'researchPaper',
      'generatedResearchPaper',
      'targetPrediction',
      'toxicityResult',
    ],
  },
  moleculeId: {
    type: String, // Store UUID for grouping items under the same molecule
    required: true,
  },
  data: {
    type: mongoose.Schema.Types.Mixed,
    required: true,
  },
  createdAt: {
    type: Date,
    default: Date.now,
  },
});

export default mongoose.model('SummarySavedItem', summarySavedItemSchema);