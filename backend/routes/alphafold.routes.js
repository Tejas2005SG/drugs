import express from 'express';
const router = express.Router();
import {
  submitPrediction,
  getUniprotSummary,
  getUniprotAnnotations,
  getJobStatus,
} from '../controllers/alphafold.controller.js';

router.post('/predict', submitPrediction);
router.get('/status/:jobId', getJobStatus);
router.get('/uniprot/summary/:uniprotId', getUniprotSummary);
router.get('/uniprot/annotations/:uniprotId', getUniprotAnnotations);
router.get('/health', (req, res) => {
  res.status(200).json({ status: 'healthy' });
});

export default router;