import mongoose from 'mongoose';

const toxicitySchema = new mongoose.Schema({
  smiles: { type: String, required: true },
  toxicityResult: {
    acuteToxicity: {
      LD50: String,
      toxicityClass: String,
    },
    endpoints: {
      hepatotoxicity: String,
      carcinogenicity: String,
    },
  },
  userId: { type: mongoose.Schema.Types.ObjectId, ref: 'User', required: true },
  created: { type: Date, default: Date.now },
});

export default mongoose.model('Toxicity', toxicitySchema);