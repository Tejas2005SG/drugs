import mongoose from "mongoose";

const savedSearchSchema = new mongoose.Schema({
  userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
  smiles: { type: String, required: true },
  targets: [{
    protein: { type: String, required: true },
    confidence: { type: Number, required: true },
    moa: { type: String, required: true },
    pathways: [{ type: String }],
    diseases: [{ type: String }],
    knownInteractions: { type: Object, default: null },
  }],
  research: [{
    title: { type: String, required: true },
    authors: { type: String, required: true },
    journal: { type: String },
    year: { type: String },
    abstract: { type: String },
    doi: { type: String },
    url: { type: String },
  }],
  docking: {
    bindingEnergy: { type: String },
    pose: { type: String },
    details: { type: String },
  },
  createdAt: { type: Date, default: Date.now },
});

export const SavedSearch = mongoose.model("SavedSearch", savedSearchSchema);