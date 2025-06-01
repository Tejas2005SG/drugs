import express from 'express';
import { predictDisease,predictTargetProtein } from '../controllers/newdrug.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();

router.post("/predictDisease",predictDisease);
router.post("/predictTargetProtein",predictTargetProtein);

export default router;