import express from 'express';
import { getSummary } from '../controllers/summary.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';

const router = express.Router();

router.get('/summary', protectRoute, getSummary);

export default router;