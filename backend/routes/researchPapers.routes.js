import express from "express";

const router = express.Router();
import {
    proxyGeminiRequest,
    saveResearchPapers,
    getSavedResearchPapers,
    checkSavedPapers,
    saveGeneratedResearchPaper,
    getSavedGeneratedResearchPapers,
    checkSavedGeneratedPapers,
} from "../controllers/researchPapers.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";



// research paper
router.post("/proxy/gemini", proxyGeminiRequest);
router.post("/save-research-papers", saveResearchPapers);
router.get("/saved-research-papers", getSavedResearchPapers);
router.get("/check-saved-papers", checkSavedPapers);
router.post("/save-generated-research-paper", saveGeneratedResearchPaper);
router.get("/saved-generated-research-papers", getSavedGeneratedResearchPapers);
router.get("/check-saved-generated-papers",  checkSavedGeneratedPapers);

export default router;