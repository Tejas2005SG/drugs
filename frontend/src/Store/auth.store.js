// auth.store.js
import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import axios from 'axios';
import { toast } from 'react-hot-toast';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode === "development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

export const useAuthStore = create(
  persist(
    (set, get) => ({
      user: null,
      loading: false,
      checkingAuth: true,
      phoneNumber: null, // Added to store phoneNumber during signup

      signup: async ({ firstName, lastName, username, email, phoneNumber, password, confirmPassword }) => {
        set({ loading: true });
        if (password !== confirmPassword) {
          set({ loading: false });
          return toast.error('Passwords do not match');
        }
        try {
          const res = await axiosInstance.post(`${API_BASE_URL}/auth/signup`, { 
            firstName, 
            lastName, 
            username, 
            email, 
            phoneNumber, 
            password, 
            confirmPassword 
          });
          set({ 
            phoneNumber: res.data.phoneNumber, // Store phoneNumber from response
            loading: false 
          });
          console.log('Signup - Updated State:', get());
          toast.success('OTP sent to your phone number');
          return res.data; // Return response for navigation in component
        } catch (error) {
          set({ loading: false });
          toast.error(error.response?.data?.message || 'An error occurred during signup');
          throw error; // Throw error to handle navigation in component
        }
      },

      verifyPhone: async ({ phoneNumber, otp }) => {
        set({ loading: true });
        try {
          const res = await axiosInstance.post(`${API_BASE_URL}/auth/verify-phone`, { 
            phoneNumber, 
            otp 
          });
          set({ 
            user: res.data.user, // Set user after successful verification
            phoneNumber: null, // Clear phoneNumber after verification
            loading: false 
          });
          console.log('VerifyPhone - Updated State:', get());
          toast.success('Phone number verified successfully');
          return res.data; // Return response for navigation in component
        } catch (error) {
          set({ loading: false });
          toast.error(error.response?.data?.message || 'OTP verification failed');
          throw error; // Throw error to handle navigation in component
        }
      },

      login: async ({ email, password }) => {
        set({ loading: true });
        try {
          const res = await axiosInstance.post(`${API_BASE_URL}/auth/login`, { email, password });
          console.log('Login Response:', res.data);
          set({ user: res.data.user, loading: false });
          console.log('Login - Updated State:', get().user);
          toast.success('Logged in successfully');
        } catch (error) {
          console.error('Login Error:', error.response?.data);
          set({ loading: false });
          toast.error(error.response?.data?.message || 'No User found');
        }
      },

      logout: async () => {
        try {
          await axiosInstance.post(`${API_BASE_URL}/auth/logout`);
          set({ user: null, phoneNumber: null }); // Clear phoneNumber on logout
          console.log('Logout - Updated State:', get().user);
          toast.success('Logged out successfully');
        } catch (error) {
          toast.error(error.response?.data?.message || 'An error occurred during logout');
        }
      },

      checkAuth: async () => {
        set({ checkingAuth: true });
        try {
          const response = await axiosInstance.get(`${API_BASE_URL}/auth/profile`);
          console.log('checkAuth Response:', response.data);
          set({ user: response.data.user, checkingAuth: false });
          console.log('checkAuth - Updated State:', get().user);
          return response.data.user;
        } catch (error) {
          console.error('checkAuth Error:', error.response?.data || error.message);
          set({ checkingAuth: false, user: null });
          console.log('checkAuth - State after error:', get().user);
          return null;
        }
      },
    }),
    {
      name: 'auth-storage', // Name for localStorage key
      partialize: (state) => ({ user: state.user }), // Only persist user, not phoneNumber
    }
  )
);