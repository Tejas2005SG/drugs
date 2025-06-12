import jwt from "jsonwebtoken";
import { User } from "../models/auth.model.js";

export const protectRoute = async (req, res, next) => {
  let token;

  token = req.cookies.token;
  // console.log("Token from cookies:", token);

  if (!token) {
    return res.status(401).json({ message: "Not authorized, no token" });
  }

  try {
    const decoded = jwt.verify(token, "SECRET_KEY_PLEXADUBAI");
    console.log("Decoded Token:", decoded);

    req.userId = decoded.userId;
    req.user = await User.findById(decoded.userId);
    console.log("User Found:", req.user);

    if (!req.user) {
      return res.status(401).json({ message: "User not found" });
    }

    next();
  } catch (error) {
    console.error("Token Verification Error:", error);
    res.status(401).json({ message: "Token is not valid" });
  }
};