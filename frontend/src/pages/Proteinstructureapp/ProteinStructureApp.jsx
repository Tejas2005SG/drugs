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
  baseURL: import.meta.mode === "development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

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
      toast.success('Structure generated successfully!', {
        style: {
          background: 'var(--color-success)',
          color: 'var(--color-primary)',
        },
      });
    } catch (err) {
      setError(err.response?.data?.message || 'Failed to generate structure');
      toast.error(err.response?.data?.message || 'Failed to generate structure', {
        style: {
          background: 'var(--color-error)',
          color: 'white',
        },
      });
    } finally {
      setLoading(false);
    }
  };

  const selectStructure = (structure) => {
    setSelectedStructure(structure);
  };

  if (checkingAuth || (!rdkitLoaded && !error)) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-primary">
        <Spinner className="text-accent" />
        <p className="mt-4 text-text-primary font-body">
          {checkingAuth ? 'Verifying authentication...' : 'Initializing molecular viewer...'}
        </p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen bg-primary">
        <div className="text-center p-8 bg-secondary rounded-lg shadow-xl max-w-md">
          <h2 className="text-2xl font-heading text-text-primary mb-4">Access Required</h2>
          <p className="text-text-secondary font-body mb-6">
            Please log in to access the Protein Structure Generator
          </p>
          <button
            className="w-full bg-accent hover:bg-accent/90 text-primary font-heading font-bold py-3 px-6 rounded-lg transition-all duration-200"
            onClick={() => (window.location.href = '/login')}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-primary p-4 md:p-8">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <header className="text-center mb-8">
          <h1 className="text-3xl md:text-4xl font-heading font-bold text-accent mb-2">
            Protein Structure Generator
          </h1>
          <p className="text-xs md:text-sm text-accent-secondary font-label font-semibold">
            (Powered by Gemini and MolMIM Nvidia Model)
          </p>
        </header>

        {/* Error Message */}
        {error && (
          <div className="bg-error/20 border border-error text-text-primary px-4 py-3 rounded-lg mb-6 flex justify-between items-center">
            <p className="font-body">{error}</p>
            <button
              className="text-text-primary hover:text-accent font-bold ml-4"
              onClick={() => setError(null)}
            >
              ✕
            </button>
          </div>
        )}

        {/* Main Content Grid */}
        <div className="grid grid-cols-1 lg:grid-cols-1 gap-6">
          {/* Left Column - Generate Form */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <h2 className="text-xl font-heading font-semibold text-accent mb-4">
                Generate New Structure
              </h2>
              <StructureForm onSubmit={handleSubmit} loading={loading} />
            </div>
          </div>

          {/* Middle Column - Structure Details */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <h2 className="text-xl font-heading font-semibold text-accent mb-4">
                Structure Details
              </h2>
              {selectedStructure ? (
                <StructureDetails structure={selectedStructure} rdkitLoaded={rdkitLoaded} />
              ) : (
                <div className="flex flex-col items-center justify-center h-96">
                  <svg
                    className="w-16 h-16 text-text-secondary mb-4"
                    fill="none"
                    stroke="currentColor"
                    viewBox="0 0 24 24"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth={1.5}
                      d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z"
                    />
                  </svg>
                  <p className="text-text-secondary font-body mb-2">
                    Select a structure to view details
                  </p>
                  <p className="text-text-secondary text-sm font-body">or generate a new one</p>
                </div>
              )}
            </div>
          </div>

          {/* Right Column - Saved Structures */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <div className="flex justify-between items-center mb-4">
                <h2 className="text-xl font-heading font-semibold text-accent">
                  Saved Structures
                </h2>
                {structures.length > 0 && (
                  <span className="bg-accent-secondary text-primary text-xs font-bold px-2 py-1 rounded-full">
                    {structures.length}
                  </span>
                )}
              </div>
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
    </div>
  );
};

export default ProteinStructureApp;