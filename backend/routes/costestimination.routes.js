// routes/costEstimationRoutes.js
import express from "express";
import { protectRoute } from "../middleware/auth.middleware.js";
import { postCostEstimation, getCostEstimation } from "../controllers/costestimination.controller.js";

const router = express.Router();

router.post("/cost-estimation", protectRoute, postCostEstimation);
router.get("/getcostestimation/:userId", protectRoute, getCostEstimation);

export default router;