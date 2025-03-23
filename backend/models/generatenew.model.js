import mongoose from "mongoose";

const  GeneratenewMoleculeSchema= new mongoose.Schema({
  smilesoffirst: {
    type: String,
    required: true,
  },
  smilesofsecond: {
    type: String,
    required: true,
  },

  information: {
    type: String, // Store the detailed information as a string
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
