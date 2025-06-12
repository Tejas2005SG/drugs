import  express from  'express';
import {getTopHeadlines, getEverything, getSources, getSavedNews} from '../controllers/news.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();


router.get('/top-headlines', getTopHeadlines);
router.get('/everything', getEverything);
router.get('/sources',getSources);
router.get('/saved', getSavedNews);


export default router;