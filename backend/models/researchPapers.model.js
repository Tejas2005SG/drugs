// researchPapers.model.js
import mongoose from "mongoose";

const researchPaperSchema = new mongoose.Schema({
  userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
  molecule: {
    title: { type: String, required: true },
    smiles: { type: String, required: true },
  },
  papers: [
    {
      title: String,
      authors: String,
      year: String,
      abstract: String,
      doi: String,
      url: String,
    },
  ],
  createdAt: { type: Date, default: Date.now },
});

const ResearchPaper = mongoose.model("ResearchPaper", researchPaperSchema);

export default ResearchPaper;