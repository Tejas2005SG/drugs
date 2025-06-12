import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import axios from 'axios';
import { toast } from 'react-hot-toast';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'https://js-backend-6k7s.onrender.com';

// IMPORTANT: Configure axios to send cookies with requests
axios.defaults.withCredentials = true;

console.log('API_BASE_URL:', API_BASE_URL); // Debug log

export const useAuthStore = create(
  persist(
    (set, get) => ({
      user: null,
      loading: false,
      checkingAuth: true,
      phoneNumber: null,

      setCheckingAuth: (status) => set({ checkingAuth: status }),

      signup: async ({ firstName, lastName, username, email, phoneNumber, password, confirmPassword }) => {
        set({ loading: true });
        if (password !== confirmPassword) {
          set({ loading: false });
          return toast.error('Passwords do not match');
        }
        try {
          console.log('Making signup request to:', `${API_BASE_URL}/api/auth/signup`);
          const res = await axios.post(`${API_BASE_URL}/api/auth/signup`, { 
            firstName, 
            lastName, 
            username, 
            email, 
            phoneNumber, 
            password, 
            confirmPassword 
          });
          set({ 
            phoneNumber: res.data.phoneNumber,
            loading: false 
          });
          console.log('Signup - Updated State:', get());
          toast.success('OTP sent to your phone number');
          return res.data;
        } catch (error) {
          console.error('Signup error:', error.response?.data || error.message);
          set({ loading: false });
          toast.error(error.response?.data?.message || 'An error occurred during signup');
          throw error;
        }
      },

      verifyPhone: async ({ phoneNumber, otp }) => {
        set({ loading: true });
        try {
          console.log('Making verify phone request to:', `${API_BASE_URL}/api/auth/verify-phone`);
          const res = await axios.post(`${API_BASE_URL}/api/auth/verify-phone`, { 
            phoneNumber, 
            otp 
          });
          set({ 
            user: res.data.user,
            phoneNumber: null,
            loading: false 
          });
          console.log('VerifyPhone - Updated State:', get());
          toast.success('Phone number verified successfully');
          return res.data;
        } catch (error) {
          console.error('Verify phone error:', error.response?.data || error.message);
          set({ loading: false });
          toast.error(error.response?.data?.message || 'OTP verification failed');
          throw error;
        }
      },

      login: async ({ email, password }) => {
        set({ loading: true });
        try {
          console.log('Making login request to:', `${API_BASE_URL}/api/auth/login`);
          const res = await axios.post(`${API_BASE_URL}/api/auth/login`, { email, password });
          console.log('Login Response:', res.data);
          set({ user: res.data.user, loading: false });
          console.log('Login - Updated State:', get().user);
          toast.success('Logged in successfully');
        } catch (error) {
          console.error('Login Error:', error.response?.data || error.message);
          set({ loading: false });
          toast.error(error.response?.data?.message || 'No User found');
        }
      },

      logout: async () => {
        try {
          console.log('Making logout request to:', `${API_BASE_URL}/api/auth/logout`);
          await axios.post(`${API_BASE_URL}/api/auth/logout`);
          set({ user: null, phoneNumber: null });
          console.log('Logout - Updated State:', get().user);
          toast.success('Logged out successfully');
        } catch (error) {
          console.error('Logout error:', error.response?.data || error.message);
          toast.error(error.response?.data?.message || 'An error occurred during logout');
        }
      },

      checkAuth: async () => {
        set({ checkingAuth: true });
        try {
          console.log('Making checkAuth request to:', `${API_BASE_URL}/api/auth/profile`);
          const response = await axios.get(`${API_BASE_URL}/api/auth/profile`);
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
      name: 'auth-storage',
      partialize: (state) => ({ user: state.user }),
    }
  )
);