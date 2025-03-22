import bcrypt from 'bcrypt';
import jwt from 'jsonwebtoken';
import { User } from '../models/auth.model.js';
import { generateTokenAndCookie } from '../utils/generateTokenAndCookie.js';

export const signup = async (req, res) => {
  const { firstName, lastName, email, password, confirmPassword } = req.body;

  try {
    if (!firstName || !lastName || !email || !password || !confirmPassword) {
      throw new Error('All fields are required');
    }

    const userAlreadyExists = await User.findOne({ email });
    if (userAlreadyExists) {
      return res.status(400).json({ success: false, message: 'User already exists' });
    }

    if (password !== confirmPassword) {
      throw new Error('Passwords do not match.');
    }

    const hashedPassword = await bcrypt.hash(password, 10);
    const user = new User({ firstName, lastName, email, password: hashedPassword, confirmPassword: hashedPassword });

    await user.save();
    generateTokenAndCookie(res, user._id);

    res.status(201).json({
      success: true,
      message: 'User created successfully',
      user: { ...user._doc, password: undefined, confirmPassword: undefined },
    });
  } catch (error) {
    res.status(400).json({ success: false, message: error.message });
  }
};

export const login = async (req, res) => {
  const { email, password } = req.body;

  try {
    const user = await User.findOne({ email });
    if (!user) return res.status(400).json({ success: false, message: 'Invalid Credentials' });

    const isPasswordValid = await bcrypt.compare(password, user.password);
    if (!isPasswordValid) return res.status(400).json({ success: false, message: 'Invalid Credentials' });

    generateTokenAndCookie(res, user._id);
    user.lastLogin = new Date();
    await user.save();

    res.status(200).json({
      success: true,
      message: 'Logged in Successfully',
      user: { ...user._doc, password: undefined, confirmPassword: undefined },
    });
  } catch (error) {
    res.status(400).json({ success: false, message: `Login Failed: ${error.message}` });
  }
};

export const logout = async (req, res) => {
  res.clearCookie('token');
  res.status(200).json({ success: true, message: 'Logged out successfully' });
};

export const getProfile = async (req, res) => {
  try {
    // req.user is set by protectRoute middleware
    if (!req.user) {
      return res.status(401).json({ success: false, message: 'Not authenticated' });
    }
    res.status(200).json({
      success: true,
      user: { ...req.user._doc, password: undefined, confirmPassword: undefined },
    });
  } catch (error) {
    console.error('Get profile error:', error);
    res.status(500).json({ success: false, message: 'Server error', error: error.message });
  }
};