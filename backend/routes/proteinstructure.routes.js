  import express from "express";
  import {
    getProteinStructure,
    postProteinStructure,
  } from "../controllers/proteinstructure.controller.js";
  import { protectRoute } from "../middleware/auth.middleware.js";

  const router = express.Router();

  router.get("/getproteinstructure/:id", protectRoute,getProteinStructure);
  router.post("/postproteinstructure/:id",protectRoute, postProteinStructure);

  export default router;
