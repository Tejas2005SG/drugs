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

// Setup allowed origins for CORS
// const allowedOrigins = [
//   process.env.CLIENT_URL,                          // From environment variable (e.g. Vercel frontend)
//                           // Local dev (Vite default)
//                             // Alternate local dev (e.g. CRA)
// ].filter(Boolean); // Removes any undefined/null

// CORS options
// const corsOptions = {
//   origin: (origin, callback) => {
//     console.log('Incoming Origin:', origin);

//     // Allow server-to-server, Postman, curl, SSR
//     if (!origin) return callback(null, true);

//     // Allow only if origin is in the allowed list
//     if (allowedOrigins.includes(origin)) {
//       callback(null, true);
//     } else {
//       console.warn('🚫 Blocked by CORS:', origin);
//       callback(new Error('Not allowed by CORS'));
//     }
//   },
//   credentials: true, // ✅ Required for cookies
//   methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS', 'PATCH'],
//   allowedHeaders: [
//     'Content-Type',
//     'Authorization',
//     'X-Requested-With',
//     'Accept',
//     'Origin'
//   ],
//   optionsSuccessStatus: 200
// };

// Apply to all incoming requests


  app.use(cors({
  origin: ['https://drugs-1-8gub.onrender.com',             // Render frontend fallback
  'http://localhost:5173'], // or use env variable
  credentials: true,
}));

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
  
});