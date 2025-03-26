import express from "express";
import dotenv from "dotenv";
import cors from "cors";
import cookieParser from "cookie-parser";
import path from "path";
import { fileURLToPath } from "url";
import { connectionDb } from "./lib/db.js";
import multer from "multer"; // Added multer for file uploads

import proteinRoutes from "./routes/proteinstructure.routes.js";
import authRoutes from "./routes/auth.routes.js";
import costestiminationRoutes from "./routes/costestimination.routes.js";
import newsRoutes from "./routes/news.routes.js";

// Configure environment variables
dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Multer setup for file uploads (store in memory)
const upload = multer({ storage: multer.memoryStorage() }); // Added multer configuration

// Middleware
const corsOptions = {
  origin: "http://localhost:5173",
  credentials: true,
};
app.use(cors(corsOptions));
app.use(express.json());
app.use(express.urlencoded({ extended: true }));
app.use(cookieParser());

// API Routes
app.use("/api/protein", upload.single("file"), proteinRoutes); // Added upload middleware for protein routes
app.use("/api/auth", authRoutes);
app.use("/api/costestimation", costestiminationRoutes);
app.use('/api/news', newsRoutes);

// Serve static files in production
if (process.env.NODE_ENV === "production") {
  app.use(express.static(path.join(__dirname, "../client/build")));
  app.get("*", (req, res) => {
    res.sendFile(path.resolve(__dirname, "../client/build", "index.html"));
  });
}

// Database connection

// Start server
app.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
});