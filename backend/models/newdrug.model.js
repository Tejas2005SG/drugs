import mongoose from "mongoose";

const predictDiseaseSchema = new mongoose.Schema({
  symptoms: {
    type: [String],
    required: true,
  },
  predictedDiseases: [
    {
      diseaseName: {
        type: String,
        required: true,
      },
      DiseaseMatchness: {
        type: String,
        required: true,
      },
      diseaseCautions: {
        type: [String],
        required : true,
      },
    },
  ],
  DiseaseAnalysis: {
    likelyDiseasesOrDisorders: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    associatedBiologicalPathways: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    affectedBiologicalProcesses: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    involvedOrgansAndTissues: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    relevantCellTypes: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    molecularCellularMechanismDisruptions: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    associatedGenes: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    knownOrPotentialTargetProteins: [
      {
        protein: String,
        source: String,
        summary: String,
        url: String,
      },
    ],
    biomarkers: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    relevantTherapeuticClassesOrDrugExamples: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
  },
    userId: {
      type: mongoose.Schema.Types.ObjectId,
      ref: 'User',
    },
  createdAt: { type: Date, default: Date.now },
});

export const PredictDisease = mongoose.model(
  "PredictDisease",
  predictDiseaseSchema
);

// Updated TargetProteinSchema
// Updated TargetProteinSchema
const TargetProteinSchema = new mongoose.Schema({
  TargetProteins: [
    {
      proteinName: {
        type: String,
        required: true,
      },
      proteinFunction: {
        type: String,
        required: true,
      },
      associatedDiseases: {
        type: [String],
        required: true,
      },
      ProtienDiscription: {
        type: String,
        required: true,
      },
      proteinDetailedDiscription: {
        type: String,
        required: true,
      },
      ProteinSmile: {
        type: String,
        required: true,
      },
      proteinUniport: {
        type: String,
        required: true,
      },
    },
  ],
  TargetLigands: [
    {
      ligandName: {
        type: String,
        required: true,
      },
      ligandFunction: {
        type: String,
        required: true,
      },
      associatedDiseases: {
        type: [String],
        required: true,
      },
      LigandDiscription: {
        type: String,
        required: true,
      },
      LigandDetailedDiscription: {
        type: String,
        required: true,
      },
      LigandSmile: {
        type: String,
        required: true,
      },
      ligandDrugBankID: { // Renamed from ligandUniport to reflect the correct identifier
        type: String,
        required: true,
      },
    },
  ],
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: 'User',
  },
  createdAt: { type: Date, default: Date.now },
});

export const TargetProtein = mongoose.model("TargetProtein", TargetProteinSchema);
