// routes/costEstimationRoutes.js
import express from "express";
import { protectRoute } from "../middleware/auth.middleware.js";
import {
  postCostEstimation,
  getCostEstimation,
  // getSymptomsAndAllProducts
} from "../controllers/costestimination.controller.js";

const router = express.Router();

router.post("/cost-estimation",  postCostEstimation);
router.get("/getcostestimation/:userId",  getCostEstimation);
// router.get("/getsymptoms-product/:id", protectRoute, getSymptomsAndAllProducts);


export default router;