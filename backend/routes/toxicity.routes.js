import express from 'express';

import { protectRoute } from '../middleware/auth.middleware.js'; // Assuming you have this from your auth setup
import { getToxicityHistory, predictToxicityController } from '../controllers/toxicity.controllers.js';

const router = express.Router();

router.post('/predict', protectRoute, predictToxicityController);
router.get('/history', protectRoute, getToxicityHistory);

export default router;