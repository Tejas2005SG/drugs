// src/routes/react.routes.ts
import express from 'express';
import { handleReaction } from '../controllers/rdkit.controller.js';

const router = express.Router();

router.post('/react', handleReaction);

export default router;
