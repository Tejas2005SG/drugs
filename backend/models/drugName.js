// src/models/DrugName.js
import mongoose from "mongoose";

const DrugNameSchema = new mongoose.Schema({
  moleculeTitle: {
    type: String,
    required: true,
  },
  smiles: {
    type: String,
    required: true,
  },
  suggestedName: {
    type: String,
    required: true,
  },
  namingDetails: {
    type: String, // Details about how the name was derived
    default: "",
  },
  createdAt: {
    type: Date,
    default: Date.now,
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
    required: true,
  },
});

export default mongoose.model("DrugName", DrugNameSchema);