import { create } from "zustand";
import axios from "axios";
import { toast } from "react-hot-toast";

const API_URL = "http://localhost:5000/api/auth";
axios.defaults.withCredentials = true;

export const useAuthStore = create((set, get) => ({
  user: null,
  token: null,
  isAuthenticated: false,
  error: null,
  isLoading: false,
  isCheckingAuth: true,

  // Persist user in local storage
  persistUser: (userData, token) => {
    if (!userData) return;
    try {
      const authData = { user: userData, token: token || userData.token };
      localStorage.setItem("authData", JSON.stringify(authData));

      if (authData.token) {
        axios.defaults.headers.common["Authorization"] = `Bearer ${authData.token}`;
      }

      set({ user: userData, token: authData.token, isAuthenticated: true });
    } catch (error) {
      console.error("Failed to persist auth data:", error);
    }
  },

  // Restore auth state from localStorage
  restoreAuth: () => {
    const storedAuth = localStorage.getItem("authData");
    if (storedAuth) {
      try {
        const { user, token } = JSON.parse(storedAuth);
        set({ user, token, isAuthenticated: true });
        axios.defaults.headers.common["Authorization"] = `Bearer ${token}`;
      } catch (error) {
        console.error("Error restoring auth:", error);
        localStorage.removeItem("authData");
      }
    }
    set({ isCheckingAuth: false });
  },

  // Signup function
  signup: async (firstName, lastName, email, password, confirmPassword) => {
    set({ isLoading: true, error: null });

    if (password !== confirmPassword) {
      set({ isLoading: false });
      toast.error("Passwords do not match");
      return;
    }

    try {
      const response = await axios.post(`${API_URL}/signup`, {
        firstName, lastName, email, password, confirmPassword,
      });

      get().persistUser(response.data.user, response.data.token);
      toast.success("Signup successful!");
    } catch (error) {
      set({ error: error.response?.data?.message || "Signup failed", isLoading: false });
      toast.error(error.response?.data?.message || "Signup failed");
    } finally {
      set({ isLoading: false });
    }
  },

  // Login function
  login: async (email, password) => {
    set({ isLoading: true, error: null });

    try {
      const response = await axios.post(`${API_URL}/login`, { email, password });

      get().persistUser(response.data.user, response.data.token);
      toast.success("Logged in successfully!");
    } catch (error) {
      set({ error: error.response?.data?.message || "Login failed", isLoading: false });
      toast.error(error.response?.data?.message || "Login failed");
    } finally {
      set({ isLoading: false });
    }
  },

  // Logout function
  logout: async () => {
    try {
      await axios.post(`${API_URL}/logout`);
      localStorage.removeItem("authData");
      set({ user: null, token: null, isAuthenticated: false });
      delete axios.defaults.headers.common["Authorization"];
      toast.success("Logged out successfully");
    } catch (error) {
      toast.error("Logout failed",error);
    }
  },
}));

// Auto-restore authentication on store initialization
useAuthStore.getState().restoreAuth();

export default useAuthStore;
