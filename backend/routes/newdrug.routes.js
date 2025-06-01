import express from 'express';
import { predictDisease,predictTargetProtein } from '../controllers/newdrug.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();

<<<<<<< HEAD
router.post("/predictDisease",protectRoute,predictDisease);
router.post("/predictTargetProtein",protectRoute,predictTargetProtein);
=======
router.post("/predictDisease",predictDisease);
router.post("/predictTargetProtein",predictTargetProtein);
>>>>>>> 7a57edbf7e80b95e0950d9d102d6807d155c80e2

export default router;