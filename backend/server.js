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
import noteRoutes from "./routes/note.routes.js";
import newdrugRoutes from "./routes/newdrug.routes.js";
import getsymptomproductRoutes from "./routes/getsymptomproduct.routes.js";
import drugNameRoutes from "./routes/drugnaming.routes.js";
import researchPaperRoutes from "./routes/researchPapers.routes.js";

dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __dirname = path.resolve();

const upload = multer({ storage: multer.memoryStorage() });

// FIXED: Correct CORS configuration
// These should be your FRONTEND URLs (where your React app is hosted)
const allowedOrigins = [
  process.env.CLIENT_URL, // Add this to your environment variables
  'https://drugs-1-8gub.onrender.com',  // Your FRONTEND URL (React app)
  'http://localhost:5173',  // Local development
  'http://localhost:3000',  // Alternative local development
].filter(Boolean);

const corsOptions = {
  origin: function (origin, callback) {
    console.log('Origin:', origin); // For debugging
    
    // Allow requests with no origin (like mobile apps, curl, or same-origin requests)
    if (!origin) return callback(null, true);
    
    // Check if origin is in allowed origins
    if (allowedOrigins.includes(origin)) {
      callback(null, true);
    } else {
      console.log('CORS blocked origin:', origin);
      callback(new Error('Not allowed by CORS'));
    }
  },
  credentials: true,
  methods: ["GET", "POST", "PUT", "DELETE", "OPTIONS", "PATCH"],
  allowedHeaders: [
    'Content-Type', 
    'Authorization', 
    'X-Requested-With',
    'Accept',
    'Origin'
  ],
  optionsSuccessStatus: 200,
  preflightContinue: false
};

// Apply CORS middleware
app.use(cors(corsOptions));

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
app.use("/api/notes", noteRoutes);
app.use("/api/newdrug", newdrugRoutes);
app.use("/api/getdata", getsymptomproductRoutes);
app.use("/api/drugname", drugNameRoutes)
app.use("/api/researchPaper", researchPaperRoutes);

app.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
  console.log(`Allowed origins: ${allowedOrigins.join(', ')}`);
});