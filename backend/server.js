import express from "express";
import dotenv from "dotenv";
import cors from "cors";
import cookieParser from "cookie-parser";
import path from "path";
import { connectionDb } from "./lib/db.js";
import multer from "multer";
import { Server } from "socket.io";
import http from "http";
import { setupSocket } from "./controllers/message.controller.js";

import proteinRoutes from "./routes/proteinstructure.routes.js";
import authRoutes from "./routes/auth.routes.js";
import costestiminationRoutes from "./routes/costestimination.routes.js";
import newsRoutes from "./routes/news.routes.js";
import messageRoutes from "./routes/message.routes.js";
import alphafoldRoutes from "./routes/alphafold.routes.js";

dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __dirname = path.resolve();

const server = http.createServer(app);

// Updated Socket.io configuration to handle both development and production
export const io = new Server(server, {
  cors: {
    origin: process.env.NODE_ENV === "production"
      ? ["https://drugs-assistant.onrender.com", process.env.CLIENT_URL]
      : [process.env.CLIENT_URL, "http://localhost:5173"],
    methods: ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    credentials: true,
  },
});

io.on("connection", (socket) => {
  console.log("A user connected:", socket.id);
  setupSocket(socket);
  socket.on("disconnect", () => {
    console.log("User disconnected:", socket.id);
  });
});

const upload = multer({ storage: multer.memoryStorage() });

// Updated CORS options to handle both development and production
const corsOptions = {
  origin: process.env.NODE_ENV === "production"
    ? ["https://drugs-assistant.onrender.com", process.env.CLIENT_URL]
    : [process.env.CLIENT_URL, "http://localhost:5173"],
  credentials: true,
  methods: ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
};

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
app.use("/api/message", messageRoutes);
app.use("/api/alphafold", alphafoldRoutes);

if (process.env.NODE_ENV === "production") {
  app.use(express.static(path.join(__dirname, "/frontend/dist")));
  app.get("*", (req, res) => {
    res.sendFile(path.resolve(__dirname, "frontend", "dist", "index.html"));
  });
}

server.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
});