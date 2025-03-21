import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Spinner } from '../../components/Spinner.jsx';
import StructureForm from '../../components/StructureForm.jsx';
import StructureList from '../../components/StructureList.jsx';
import StructureDetails from '../../components/StructureDetails.jsx';
import { toast } from 'react-hot-toast';

// Configure axios with base URL
const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000';
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
});

const ProteinStructureApp = () => {
  const [structures, setStructures] = useState([]);
  const [selectedStructure, setSelectedStructure] = useState(null);
  const [loading, setLoading] = useState(false);
  const [rdkitLoaded, setRdkitLoaded] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    const loadRDKit = async () => {
      try {
        if (window.RDKit) {
          setRdkitLoaded(true);
          console.log('RDKit already loaded:', window.RDKit.version());
          return;
        }
        console.log('Loading RDKit script...');
        await new Promise((resolve, reject) => {
          const script = document.createElement('script');
          script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
          script.async = true;
          script.onload = () => {
            console.log('RDKit script loaded');
            resolve();
          };
          script.onerror = () => {
            console.error('Failed to load RDKit script');
            reject(new Error('Failed to load RDKit script'));
          };
          document.head.appendChild(script);
        });
        console.log('Initializing RDKit module...');
        if (!window.initRDKitModule) {
          throw new Error('initRDKitModule not defined');
        }
        const rdkitModule = await window.initRDKitModule();
        if (!rdkitModule) {
          throw new Error('RDKit module initialization failed');
        }
        window.RDKit = rdkitModule;
        console.log('RDKit module initialized:', window.RDKit.version());
        setRdkitLoaded(true);
      } catch (err) {
        console.error('RDKit initialization failed:', err);
        setError('Failed to load molecular visualization library: ' + err.message);
      }
    };

    const initializeApp = async () => {
      console.log('Initializing app...');
      await loadRDKit();
    };

    initializeApp();
  }, []);

  const handleSubmit = async (formData) => {
    try {
      setLoading(true);
      setError(null);
      console.log('Submitting form data:', formData);
      const response = await axiosInstance.post('/api/protein/postproteinstructure', formData);
      console.log('Structure posted:', response.data);
      setStructures((prev) => [response.data.structure, ...prev]);
      setSelectedStructure(response.data.structure);
      toast.success('Structure generated successfully!');
    } catch (err) {
      console.error('Post structure error:', err.response?.data || err.message);
      const errorMessage = err.response?.data?.message || 'Failed to generate structure';
      const errorDetails = err.response?.data?.error ? `: ${JSON.stringify(err.response.data.error)}` : '';
      setError(`${errorMessage}${errorDetails}`);
    } finally {
      setLoading(false);
    }
  };

  const selectStructure = (structure) => {
    setSelectedStructure(structure);
  };

  if (!rdkitLoaded && !error) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="text-center">
          <Spinner />
          <p className="mt-4 text-gray-600">Initializing molecular viewer...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="container mx-auto px-4 py-8">
      <h1 className="text-3xl font-bold text-gray-800 mb-8">Protein Structure Generator</h1>

      {error && (
        <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
          <p>{error}</p>
          <button
            className="text-red-700 font-bold underline ml-2"
            onClick={() => setError(null)}
          >
            Dismiss
          </button>
        </div>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        <div className="lg:col-span-1">
          <div className="bg-white p-6 rounded-lg shadow-md mb-8">
            <h2 className="text-xl font-semibold mb-4">Generate New Structure</h2>
            <StructureForm onSubmit={handleSubmit} loading={loading} />
          </div>

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

        <div className="lg:col-span-2">
          <div className="bg-white p-6 rounded-lg shadow-md">
            {selectedStructure ? (
              <StructureDetails
                structure={selectedStructure}
                rdkitLoaded={rdkitLoaded}
              />
            ) : (
              <div className="flex flex-col items-center justify-center h-96">
                <p className="text-gray-500 mb-4">Select a structure to view details</p>
                <p className="text-gray-400 text-sm">or generate a new one</p>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureApp;