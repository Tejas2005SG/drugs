import express from 'express';
import { login, logout, signup, getProfile } from '../controllers/auth.controller.js'; // Add getProfile
import { protectRoute } from '../middleware/auth.middleware.js'; // Import protectRoute

const router = express.Router();

router.post('/signup', signup);
router.post('/login', login);
router.post('/logout', logout);
router.get('/profile', protectRoute, getProfile); // Add this line

export default router;