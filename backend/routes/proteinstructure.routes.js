  import express from "express";
  import {
    getProteinStructure,
    postProteinStructure,
  } from "../controllers/proteinstructure.controller.js";

  const router = express.Router();

  router.get("/getproteinstructure", getProteinStructure);
  router.post("/postproteinstructure", postProteinStructure);

  export default router;
