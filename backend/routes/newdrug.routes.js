import express from 'express';
import { predictDisease,predictTargetProtein } from '../controllers/newdrug.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();

router.post("/predictDisease",protectRoute,predictDisease);
router.post("/predictTargetProtein",protectRoute,predictTargetProtein);

export default router;