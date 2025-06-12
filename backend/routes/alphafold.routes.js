import express from 'express';
const router = express.Router();
import {
  submitPrediction,
  getUniprotSummary,
  getUniprotAnnotations,
  getJobStatus,
  getPreviousJobs,
  // getPDBFile
} from '../controllers/alphafold.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';

router.post('/predict', submitPrediction);
router.get('/status/:jobId',  getJobStatus);
router.get('/uniprot/summary/:uniprotId',getUniprotSummary);
router.get('/uniprot/annotations/:uniprotId',  getUniprotAnnotations);
router.get('/previous-jobs',  getPreviousJobs);
// router.get('/pdb/:jobId',  getPDBFile);

router.get('/health', (req, res) => {
  res.status(200).json({ status: 'healthy' });
});

export default router;