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

} from "../controllers/proteinstructure.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";


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



export default router;