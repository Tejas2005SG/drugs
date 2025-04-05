import express from "express";
import {
  getProteinStructure,
  postProteinStructure,
  generatenewmolecule,
  getgeneratednewmolecule,
  proxyGeminiRequest,
  saveResearchPapers,
  getSavedResearchPapers,
  checkSavedPapers,
  saveGeneratedResearchPaper,
  getSavedGeneratedResearchPapers,
  checkSavedGeneratedPapers,
  convertFileToSmiles, // New controller for converting MOL/SDF to SMILES
  getFingerprints, // New controller for fingerprint extraction (mocked)
  performDocking, // New controller for molecular docking (mocked)
  saveSearch, // New controller for saving searches
  getSavedSearches, // New controller for retrieving saved searches
  checkSavedSearches,
  generateDrugName,
  getSavedDrugNames,
  checkSavedDrugName,
  acceptDrugName,
  savePendingDrugName, // New controller for checking if a search exists
} from "../controllers/proteinstructure.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";
import drugName from "../models/drugName.js";

const router = express.Router();

// Existing Routes
router.get("/getproteinstructure/:id", protectRoute, getProteinStructure);
router.post("/postproteinstructure/:id", protectRoute, postProteinStructure);
router.post("/generatenewmolecule/:id", protectRoute, generatenewmolecule);
router.get("/generatednewmolecule", protectRoute, getgeneratednewmolecule);

// research paper
router.post("/proxy/gemini", proxyGeminiRequest);
router.post("/save-research-papers", protectRoute, saveResearchPapers);
router.get("/saved-research-papers", protectRoute, getSavedResearchPapers);
router.get("/check-saved-papers", protectRoute, checkSavedPapers);
router.post("/save-generated-research-paper", protectRoute, saveGeneratedResearchPaper);
router.get("/saved-generated-research-papers", protectRoute, getSavedGeneratedResearchPapers);
router.get("/check-saved-generated-papers", protectRoute, checkSavedGeneratedPapers);

// New Routes for AI-Driven Target Prediction
router.post("/convert-file-to-smiles", protectRoute, convertFileToSmiles); // Convert MOL/SDF to SMILES
router.post("/rdkit-fingerprints", protectRoute, getFingerprints); // Extract fingerprints (mocked)
router.post("/docking", protectRoute, performDocking); // Perform molecular docking (mocked)
router.post("/save-search", protectRoute, saveSearch); // Save search results
router.get("/saved-searches", protectRoute, getSavedSearches); // Retrieve saved searches
router.get("/check-saved-searches", protectRoute, checkSavedSearches); // Check if a search exists



// ai-naming
router.post("/generate-drug-name/:id", protectRoute, generateDrugName);
router.post("/accept-drug-name/:id", protectRoute, acceptDrugName);
router.post("/save-pending-drug-name/:id", protectRoute, savePendingDrugName);
router.get("/saved-drug-names", protectRoute, getSavedDrugNames);
router.get("/check-saved-drug-name", protectRoute, checkSavedDrugName);
router.delete("/delete-drug-name/:id", protectRoute, async (req, res) => {
  try {
    const drugNameId = req.params.id;
    await drugName.findByIdAndDelete(drugNameId);
    res.status(200).json({ message: "Drug name deleted successfully" });
  } catch (error) {
    res.status(500).json({ message: "Error deleting drug name", error: error.message });
  }
});
export default router;