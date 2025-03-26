import express from "express";
import dotenv from "dotenv";
import cors from "cors";
import cookieParser from "cookie-parser";
import path from "path";
import { fileURLToPath } from "url";
import { connectionDb } from "./lib/db.js";
import { Server } from "socket.io";
import http from "http";
import { setupSocket } from "./controllers/message.controller.js";

import proteinRoutes from "./routes/proteinstructure.routes.js";
import authRoutes from "./routes/auth.routes.js";
import costestiminationRoutes from "./routes/costestimination.routes.js";
import newsRoutes from "./routes/news.routes.js";
import messageRoutes from "./routes/message.routes.js";

dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Create HTTP server
const server = http.createServer(app);
// Initialize Socket.IO
export const io = new Server(server, {
  cors: {
    origin: process.env.CLIENT_URL||"http://localhost:5173", // Match your frontend URL (Vite default port)
    methods: ["GET", "POST"],
    credentials: true,
  },
});

// Socket.IO connection
io.on("connection", (socket) => {
  console.log("A user connected:", socket.id);
  setupSocket(socket); // Set up socket events from controller

  socket.on("disconnect", () => {
    console.log("User disconnected:", socket.id);
  });
});

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
app.use("/api/protein", proteinRoutes);
app.use("/api/auth", authRoutes);
app.use("/api/costestimation", costestiminationRoutes);
app.use("/api/news", newsRoutes);
app.use("/api/message", messageRoutes);

// Production static file serving
if (process.env.NODE_ENV === "production") {
  app.use(express.static(path.join(__dirname, "../client/build")));
  app.get("*", (req, res) => {
    res.sendFile(path.resolve(__dirname, "../client/build", "index.html"));
  });
}

// Start server
server.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
});