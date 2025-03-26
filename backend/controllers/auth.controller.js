import bcrypt from 'bcrypt';
import jwt from 'jsonwebtoken';
import { User } from '../models/auth.model.js'; // Update the import to use user.model.js
import { generateTokenAndCookie } from '../utils/generateTokenAndCookie.js';

export const signup = async (req, res) => {
  const { firstName, lastName, username, email, password, confirmPassword } = req.body;

  try {
    if (!firstName || !lastName || !username || !email || !password || !confirmPassword) {
      throw new Error('All fields are required');
    }

    const userAlreadyExists = await User.findOne({ $or: [{ email }, { username }] });
    if (userAlreadyExists) {
      return res.status(400).json({
        success: false,
        message: userAlreadyExists.email === email ? 'Email already exists' : 'Username already exists',
      });
    }

    if (password !== confirmPassword) {
      throw new Error('Passwords do not match.');
    }

    const hashedPassword = await bcrypt.hash(password, 10);
    const user = new User({
      firstName,
      lastName,
      username,
      email,
      password: hashedPassword,
    });

    await user.save();
    generateTokenAndCookie(res, user._id);

    res.status(201).json({
      success: true,
      message: 'User created successfully',
      user: { ...user._doc, password: undefined },
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
      user: { ...user._doc, password: undefined },
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
    if (!req.user) {
      return res.status(401).json({ success: false, message: 'Not authenticated' });
    }
    res.status(200).json({
      success: true,
      user: { ...req.user._doc, password: undefined },
    });
  } catch (error) {
    console.error('Get profile error:', error);
    res.status(500).json({ success: false, message: 'Server error', error: error.message });
  }
};