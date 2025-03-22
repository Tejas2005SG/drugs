import mongoose from "mongoose";

const userSchema = new mongoose.Schema(
  {
    firstName: {
      type: String,
      required: [true, "Name is required"],
    },
    lastName: {
      type: String,
      required: [true, "Name is required"],
    },
    email: {
      type: String,
      required: true,
      unique: true,
    },
    password: {
      type: String,
      required: true,
    },
    confirmPassword: {
      type: String,
      required: true,
    },
    
    lastLogin: {
      type: Date,
      default: Date.now,
    },
    protienStructures: [
      {
        type: mongoose.Schema.Types.ObjectId,
        ref: "ProteinStructure",
      },
    ],
   
   
  },
  { timestamps: true }
);

export const User = mongoose.model("User", userSchema);
