import mongoose from 'mongoose';

const toxicitySchema = new mongoose.Schema({
  smiles: { type: String, required: true },
  symptoms: { type: String, required: true },
  toxicityResult: {
    smiles: String,
    symptoms: String,
    acuteToxicity: {
      LD50: String,
      toxicityClass: String,
    },
    endpoints: {
      hepatotoxicity: String,
      carcinogenicity: String,
      cardiotoxicity: String,
      neurotoxicity: String,
      nephrotoxicity: String,
      skin_sensitization: String,
    },
    sideEffects: [
      {
        name: String,
        description: String,
        severity: String,
      },
    ],
  },
  geminiAnalysis: {
    isInputValid: { type: Boolean, default: false },
    inputError: { type: String, default: null },
    inferredCompound: { type: String, default: null },
    inferredSmiles: { type: String, default: null },
    errorDetails: {
      reason: { type: String, default: null },
      suggestions: [{ type: String }],
    },
    mechanisms: [
      {
        feature: String,
        pathway: String,
      },
    ],
    structureConcerns: [
      {
        substructure: String,
        similarCompound: String,
        concern: String,
      },
    ],
    safetyRecommendations: [
      {
        test: String,
        modification: String,
      },
    ],
  },
  userId: { type: mongoose.Schema.Types.ObjectId, ref: 'User', required: true },
  created: { type: Date, default: Date.now },
});


export default mongoose.model('Toxicity', toxicitySchema);
