import express from 'express';
import { bulkSaveItems, getNewMoleculeProgress, getSavedItems, getSummary, saveItem } from '../controllers/summary.controllers.js';
import { protectRoute } from '../middleware/auth.middleware.js';

const router = express.Router();

router.get('/summary', protectRoute, getSummary);
router.post('/save', protectRoute, saveItem);

// Route to get all saved items
router.get('/saved', protectRoute, getSavedItems);
router.post('/bulk-save', protectRoute, bulkSaveItems);
router.get('/new-molecule-progress', protectRoute, getNewMoleculeProgress);
export default router;