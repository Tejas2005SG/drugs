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
// import jarvisRoutes from "./routes/jarvis.routes.js"
import researchPaperRoutes from "./routes/researchPapers.routes.js";
dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __dirname = path.resolve();

const upload = multer({ storage: multer.memoryStorage() });

// Dynamic CORS configuration for Render
const allowedOrigins = [
  
   // React default port
  'https://drugs-1-8gub.onrender.com',  // Vite default port
  'http://localhost:5173',  // Vite alternate port'
  // Remove the backend URL - it shouldn't be in allowed origins
].filter(Boolean); // Remove any undefined values

const corsOptions = {
  origin: function (origin, callback) {
    console.log('Origin:', origin); // For debugging
    
    // Allow requests with no origin (like mobile apps, curl, or same-origin requests)
    if (!origin) return callback(null, true);
    
    // Check if origin is in allowed origins
    if (allowedOrigins.includes(origin)) {
      callback(null, true);
    } else {
      // For additional flexibility, you can also check for specific patterns
      const isAllowed = allowedOrigins.some(allowedOrigin => {
        // Check for exact match or subdomain match
        return origin === allowedOrigin || 
               (allowedOrigin && origin.endsWith(allowedOrigin.replace(/^https?:\/\//, '')));
      });
      
      if (isAllowed) {
        callback(null, true);
      } else {
        console.log('CORS blocked origin:', origin);
        callback(new Error('Not allowed by CORS'));
      }
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
  // Add these for better compatibility
  optionsSuccessStatus: 200, // For legacy browser support
  preflightContinue: false
};

// Apply CORS middleware
app.use(cors(corsOptions));

// Handle preflight requests is already handled by the cors middleware above
// You don't need the additional app.options('*', cors(corsOptions));

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
app.use("/api/toxicity", toxicityRoutes);
app.use("/api/notes", noteRoutes);
app.use("/api/newdrug", newdrugRoutes);
app.use("/api/getdata", getsymptomproductRoutes);
app.use("/api/drugname", drugNameRoutes)
app.use("/api/researchPaper", researchPaperRoutes);

// // app.use("/api/jarvis",jarvisRoutes);

// if (process.env.NODE_ENV === "production") {
//   app.use(express.static(path.join(__dirname, "/frontend/dist")));
//   app.get("*", (req, res) => {
//     res.sendFile(path.resolve(__dirname, "frontend", "dist", "index.html"));
//   });


app.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
  console.log(`Allowed origins: ${allowedOrigins.join(', ')}`);
});