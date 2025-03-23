  import express from "express";
  import {
    getProteinStructure,
    postProteinStructure,
    generatenewmolecule,
    getgeneratednewmolecule,
    getLatestGeneratedMolecule,
  } from "../controllers/proteinstructure.controller.js";
  import { protectRoute } from "../middleware/auth.middleware.js";

  const router = express.Router();

  router.get("/getproteinstructure/:id", protectRoute,getProteinStructure);
  router.post("/postproteinstructure/:id",protectRoute, postProteinStructure);
  router.post("/generatenewmolecule/:id",protectRoute, generatenewmolecule);
  router.get("/generatednewmolecule/:id",protectRoute, getgeneratednewmolecule);
  router.get('/generatednewmolecule/latest/:userId', protectRoute, getLatestGeneratedMolecule);

  export default router;
