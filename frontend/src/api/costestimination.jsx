// src/api/costEstimation.js
import axios from 'axios';

const API_URL = 'http://localhost:5000/api/costestimation'; // Matches server.js PORT and /api prefix

const axiosInstance = axios.create({
  baseURL: API_URL,
  withCredentials: true, // Enable cookies/credentials for CORS
});

export const postCostEstimation = async (smiles) => {
  const response = await axiosInstance.post('/cost-estimation', { smiles });
  return response.data;
};

export const getCostEstimations = async (userId) => {
  const response = await axiosInstance.get(`/getcostestimation/${userId}`);
  return response.data;
};