import express from "express";
import { sendMessage, getMessages, getUsers } from "../controllers/message.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";

const router = express.Router();

router.post("/send", protectRoute, sendMessage);
router.get("/:userId", protectRoute, getMessages);
router.get("/users/list", protectRoute, getUsers);

export default router;