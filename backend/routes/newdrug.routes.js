import express from 'express';
import { predictDisease,predictTargetProtein, getnewdrug,getSymptoms } from '../controllers/newdrug.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();

router.post("/predictDisease/:id",protectRoute,predictDisease);
router.post("/predictTargetProtein",protectRoute,predictTargetProtein);
router.get('/getnewdrug/:id',protectRoute, getnewdrug);
router.get('/symptoms/:id',protectRoute, getSymptoms);


export default router;