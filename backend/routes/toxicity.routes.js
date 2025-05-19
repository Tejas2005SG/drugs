import express from 'express';

import { protectRoute } from '../middleware/auth.middleware.js'; // Assuming you have this from your auth setup
import { getToxicityHistory, predictToxicityController, saveGeminiAnalysis } from '../controllers/toxicity.controllers.js';

const router = express.Router();

router.post('/predict', protectRoute, predictToxicityController);
router.get('/history', protectRoute, getToxicityHistory);
router.post('/save-analysis', protectRoute, saveGeminiAnalysis);
export default router;