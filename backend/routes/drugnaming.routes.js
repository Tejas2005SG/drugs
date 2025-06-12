import express from "express";
import {
  generateDrugName,
  getSavedDrugNames,
  checkSavedDrugName,
  acceptDrugName,
  savePendingDrugName, // New controller for checking if a search exists
  deleteDrugName,
} from "../controllers/drugnaming.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";


const router = express.Router();



// ai-naming
router.post("/generate-drug-name/:id", generateDrugName);
router.post("/accept-drug-name/:id",  acceptDrugName);
router.post("/save-pending-drug-name/:id", savePendingDrugName);
router.get("/saved-drug-names",getSavedDrugNames);
router.get("/check-saved-drug-name", checkSavedDrugName);
router.delete("/delete-drug-name/:id",  deleteDrugName);
export default router;