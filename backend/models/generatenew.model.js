// GeneratenewMoleculeSchema.js
import mongoose from "mongoose";

const GeneratenewMoleculeSchema = new mongoose.Schema({
  smilesoffirst: {
    type: String,
    required: true,
  },
  smilesofsecond: {
    type: String,
    required: true,
  },
  newmoleculetitle: {
    type: String,
    required: true,
  },
  newSmiles: {
    type: String,
    default: "",
  },
  newIupacName: {
    type: String,
    default: "",
  },
  conversionDetails: {
    type: String,
    default: "",
  },
  potentialDiseases: {
    type: String,
    default: "",
  },
  information: {
    type: String, // Stores the raw streaming response as text
    default: "",
  },
  created: {
    type: Date,
    default: Date.now,
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
  },
});

export default mongoose.model("GeneratenewMolecule", GeneratenewMoleculeSchema);