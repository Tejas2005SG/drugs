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
    
  } from "../controllers/proteinstructure.controller.js";
  import { protectRoute } from "../middleware/auth.middleware.js";

  const router = express.Router();

  router.get("/getproteinstructure/:id", protectRoute,getProteinStructure);
  router.post("/postproteinstructure/:id",protectRoute, postProteinStructure);
  router.post("/generatenewmolecule/:id",protectRoute, generatenewmolecule);
  router.get("/generatednewmolecule", protectRoute, getgeneratednewmolecule);
  router.get("/generatednewmolecule",protectRoute, getgeneratednewmolecule);
  // router.get('/generatednewmolecule/latest/:userId', protectRoute, getLatestGeneratedMolecule);

router.post("/proxy/gemini", proxyGeminiRequest);
router.post("/save-research-papers", protectRoute, saveResearchPapers);
router.get("/saved-research-papers", protectRoute, getSavedResearchPapers);

router.get("/check-saved-papers", protectRoute, checkSavedPapers);

router.post("/save-generated-research-paper", protectRoute, saveGeneratedResearchPaper);
router.get("/saved-generated-research-papers", protectRoute, getSavedGeneratedResearchPapers);
router.get("/check-saved-generated-papers", protectRoute, checkSavedGeneratedPapers);
  export default router;
