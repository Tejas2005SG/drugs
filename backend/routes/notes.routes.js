import express from 'express';
import { getNotes, createNotes, updateNote, deleteNote } from '../controllers/notes.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();

router.post("/createnotes", protectRoute, createNotes);
router.get("/getnotes", protectRoute, getNotes);
router.put("/updatenotes/:id", protectRoute, updateNote);
router.delete("/deletenotes/:id", protectRoute, deleteNote);

export default router;