// routes/summary.routes.js
import express from 'express';
import { bulkSaveItems, getMoleculeProgress, getSavedItems, getSummary, saveItem } from '../controllers/summary.controllers.js';
import { protectRoute } from '../middleware/auth.middleware.js';

const router = express.Router();

router.get('/summary', getSummary);
router.post('/save', saveItem);
router.get('/saved',  getSavedItems);
router.post('/bulk-save', bulkSaveItems);
router.get('/molecule-progress',getMoleculeProgress);

export default router;