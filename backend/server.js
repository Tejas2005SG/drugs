import express from "express";
import dotenv from "dotenv";
import cors from "cors";
import cookieParser from "cookie-parser";
import path from "path";
import { connectionDb } from "./lib/db.js";
import multer from "multer";

import proteinRoutes from "./routes/proteinstructure.routes.js";
import authRoutes from "./routes/auth.routes.js";
import costestiminationRoutes from "./routes/costestimination.routes.js";
import newsRoutes from "./routes/news.routes.js";
import alphafoldRoutes from "./routes/alphafold.routes.js";
import toxicityRoutes from "./routes/toxicity.routes.js"
import summaryRoutes from "./routes/summary.routes.js";
dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __dirname = path.resolve();

const upload = multer({ storage: multer.memoryStorage() });

// Dynamic CORS configuration for Render
const allowedOrigins = [
  process.env.CLIENT_URL,
  process.env.RENDER_EXTERNAL_URL, // Render provides this automatically
  'http://localhost:5173' // For local development
].filter(Boolean); // Remove any undefined values

const corsOptions = {
  origin: function (origin, callback) {
    // Allow requests with no origin (like mobile apps or curl requests)
    if (!origin) return callback(null, true);
    
    if (allowedOrigins.some(allowedOrigin => 
      origin === allowedOrigin || 
      origin.startsWith(allowedOrigin) ||
      origin.includes('render.com')
    ) ){
      callback(null, true);
    } else {
      callback(new Error('Not allowed by CORS'));
    }
  },
  credentials: true,
  methods: ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
  allowedHeaders: ['Content-Type', 'Authorization']
};

// Apply CORS middleware
app.use(cors(corsOptions));

// Handle preflight requests for all routes
app.options('*', cors(corsOptions));

app.use(express.json());
app.use(express.urlencoded({ extended: true }));
app.use(cookieParser());

// Serve static files from the 'jobs' directory under '/results'
app.use('/results', express.static(path.join(__dirname, 'jobs')));

// API Routes
app.use("/api/protein", upload.single("file"), proteinRoutes);
app.use("/api/auth", authRoutes);
app.use("/api/costestimation", costestiminationRoutes);
app.use("/api/news", newsRoutes);
app.use("/api/alphafold", alphafoldRoutes);
app.use("/api/toxicity", toxicityRoutes)
app.use("/api/summary", summaryRoutes);

if (process.env.NODE_ENV === "production") {
  app.use(express.static(path.join(__dirname, "/frontend/dist")));
  app.get("*", (req, res) => {
    res.sendFile(path.resolve(__dirname, "frontend", "dist", "index.html"));
  });
}

app.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
  console.log(`Allowed origins: ${allowedOrigins.join(', ')}`);
});