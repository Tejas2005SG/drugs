import express from 'express';
import { protectRoute } from '../middleware/auth.middleware.js';
import {
  getToxicityHistory,
  predictToxicityController,
  saveGeminiAnalysis,
  getGeminiAnalysis,
} from '../controllers/toxicity.controllers.js';

const router = express.Router();

router.post('/predict',  predictToxicityController);
router.get('/history',  getToxicityHistory);
router.post('/save-analysis',  saveGeminiAnalysis);
router.post('/gemini-analysis',  getGeminiAnalysis);

export default router;