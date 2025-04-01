import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Spinner } from '../../components/Spinner.jsx';
import StructureForm from '../../components/StructureForm.jsx';
import StructureList from '../../components/StructureList.jsx';
import StructureDetails from '../../components/StructureDetails.jsx';
import { toast } from 'react-hot-toast';
import { useAuthStore } from '../../Store/auth.store.js';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode==="development" ? API_BASE_URL : '/api',
  withCredentials: true,
});
// baseURL: import.meta.mode === "development" ? "http://localhost:5000/api" : "/api",
const ProteinStructureApp = () => {
  const [structures, setStructures] = useState([]);
  const [selectedStructure, setSelectedStructure] = useState(null);
  const [loading, setLoading] = useState(false);
  const [rdkitLoaded, setRdkitLoaded] = useState(false);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const loadRDKit = async () => {
      try {
        if (window.RDKit) {
          setRdkitLoaded(true);
          return;
        }
        await new Promise((resolve, reject) => {
          const script = document.createElement('script');
          script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
          script.async = true;
          script.onload = () => resolve();
          script.onerror = () => reject(new Error('Failed to load RDKit script'));
          document.head.appendChild(script);
        });
        const rdkitModule = await window.initRDKitModule();
        window.RDKit = rdkitModule;
        setRdkitLoaded(true);
      } catch (err) {
        setError('Failed to load RDKit: ' + err.message);
      }
    };

    const fetchStructures = async () => {
      try {
        setLoading(true);
        const response = await axiosInstance.get(`/protein/getproteinstructure/${user._id}`);
        setStructures(response.data);
      } catch (err) {
        setError(err.response?.data?.message || 'Failed to fetch structures');
      } finally {
        setLoading(false);
      }
    };

    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError('Authentication failed. Please log in.');
        return;
      }
      await loadRDKit();
      await fetchStructures();
    };

    initializeApp();
  }, [checkAuth]);

  const handleSubmit = async (formData) => {
    if (!user?._id) {
      setError('Please log in to generate structures');
      return;
    }
    try {
      setLoading(true);
      setError(null);
      const response = await axiosInstance.post(`/protein/postproteinstructure/${user._id}`, formData);
      setStructures((prev) => [response.data.structure, ...prev]);
      setSelectedStructure(response.data.structure);
      toast.success('Structure generated successfully!');
    } catch (err) {
      setError(err.response?.data?.message || 'Failed to generate structure');
    } finally {
      setLoading(false);
    }
  };

  const selectStructure = (structure) => {
    setSelectedStructure(structure);
  };

  if (checkingAuth || (!rdkitLoaded && !error)) {
    return (
      <div className="flex items-center justify-center h-screen">
        <Spinner />
        <p className="mt-4 text-gray-600">
          {checkingAuth ? 'Verifying authentication...' : 'Initializing molecular viewer...'}
        </p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="text-center">
          <p className="text-gray-600">Please log in to access the Protein Structure Generator</p>
          <button
            className="mt-4 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            onClick={() => window.location.href = '/login'}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="container mx-auto px-4 py-8">
      <h1 className="text-3xl text-center font-bold text-gray-800 mb-8">Protein Structure Generator

      <p className="text-xs  p-1 text-blue-700 font-semibold">(Powered by Gemini and MolMIM Nvidia Model)</p>

      </h1>
      {error && (
        <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
          <p>{error}</p>
          <button className="text-red-700 underline ml-2" onClick={() => setError(null)}>
            Dismiss
          </button>
        </div>
      )}
      <div className="grid grid-cols-1 gap-8">
        {/* Section 1: Generate New Structure */}
        <div>
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h2 className="text-xl font-semibold mb-4">Generate New Structure</h2>
            <StructureForm onSubmit={handleSubmit} loading={loading} />
          </div>
        </div>

        {/* Section 2: Select a structure to view details */}
        <div>
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h2 className="text-xl font-semibold mb-4">Structure Details</h2>
            {selectedStructure ? (
              <StructureDetails structure={selectedStructure} rdkitLoaded={rdkitLoaded} />
            ) : (
              <div className="flex flex-col items-center justify-center h-96">
                <p className="text-gray-500 mb-4">Select a structure to view details</p>
                <p className="text-gray-400 text-sm">or generate a new one</p>
              </div>
            )}
          </div>
        </div>

        {/* Section 3: Saved Structures */}
        <div>
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h2 className="text-xl font-semibold mb-4">Saved Structures</h2>
            <StructureList
              structures={structures}
              onSelect={selectStructure}
              selected={selectedStructure}
              loading={loading}
            />
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureApp;