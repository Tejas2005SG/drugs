import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Spinner } from '../components/Spinner.jsx'; // Adjust path
import { toast } from 'react-hot-toast';
import { useAuthStore } from '../Store/auth.store.js'; // Adjust path

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000';
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const ProteinStructureEvolution = () => {
  const [smilesFirst, setSmilesFirst] = useState('');
  const [smilesSecond, setSmilesSecond] = useState('');
  const [newSmiles, setNewSmiles] = useState('');
  const [newIupacName, setNewIupacName] = useState('');
  const [moleculeId, setMoleculeId] = useState('');
  const [fetchMoleculeId, setFetchMoleculeId] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError('Authentication failed. Please log in.');
        return;
      }
      // Fetch the latest generated molecule after authentication
      await fetchLatestMolecule();
    };
    initializeApp();
  }, [checkAuth]);

  // Function to fetch the latest generated molecule
  const fetchLatestMolecule = async () => {
    if (!user?._id) return;

    setLoading(true);
    try {
      console.log(`Fetching latest molecule for user ${user._id}`);
      // Fetch all molecules for the user and sort by creation date
      const response = await axiosInstance.get(`/api/protein/generatednewmolecule/latest/${user._id}`);
      const { newSmiles, newIupacName, id } = response.data;
      setNewSmiles(newSmiles);
      setNewIupacName(newIupacName);
      setMoleculeId(id);
      setFetchMoleculeId(id); // Pre-fill the fetch input with the latest ID
    } catch (err) {
      console.error('Error fetching latest molecule:', err.response?.data || err.message);
      // Ignore if no molecule exists yet (404), don’t set error
      if (err.response?.status !== 404) {
        setError(err.response?.data?.message || 'Failed to fetch latest molecule');
      }
    } finally {
      setLoading(false);
    }
  };

  const handleGenerate = async (e) => {
    e.preventDefault();
    setError(null);
    setNewSmiles('');
    setNewIupacName('');
    setMoleculeId('');
    setLoading(true);

    if (!user?._id) {
      setError('Please log in to generate a new molecule');
      setLoading(false);
      return;
    }

    try {
      console.log(`Sending POST to /api/protein/generatenewmolecule/${user._id}`);
      const response = await axiosInstance.post(`/api/protein/generatenewmolecule/${user._id}`, {
        smilesoffirst: smilesFirst,
        smilesofsecond: smilesSecond,
      });

      const { newSmiles, newIupacName, id } = response.data;
      setNewSmiles(newSmiles);
      setNewIupacName(newIupacName);
      setMoleculeId(id);
      setFetchMoleculeId(id); // Update fetch input with new ID
      toast.success('New molecule generated successfully!');
    } catch (err) {
      console.error('Frontend error:', err.response?.data || err.message);
      setError(err.response?.data?.message || 'Failed to generate molecule');
    } finally {
      setLoading(false);
    }
  };

  const handleFetch = async () => {
    if (!fetchMoleculeId) {
      setError('Please enter a molecule ID to fetch');
      return;
    }

    setError(null);
    setNewSmiles('');
    setNewIupacName('');
    setLoading(true);

    try {
      console.log(`Sending GET to /api/protein/generatednewmolecule/${fetchMoleculeId}`);
      const response = await axiosInstance.get(`/api/protein/generatednewmolecule/${fetchMoleculeId}`);
      const { newSmiles, newIupacName } = response.data;
      setNewSmiles(newSmiles);
      setNewIupacName(newIupacName);
      setMoleculeId(fetchMoleculeId);
      toast.success('Molecule fetched successfully!');
    } catch (err) {
      console.error('Frontend error:', err.response?.data || err.message);
      setError(err.response?.data?.message || 'Failed to fetch molecule');
    } finally {
      setLoading(false);
    }
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center h-screen">
        <Spinner />
        <p className="mt-4 text-gray-600">Verifying authentication...</p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="text-center">
          <p className="text-gray-600">Please log in to access Protein Structure Evolution</p>
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
      <h1 className="text-3xl font-bold text-gray-800 mb-8 text-center">
        Protein Structure Evolution
      </h1>

      {error && (
        <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-6">
          <p>{error}</p>
          <button
            className="text-red-700 underline ml-2"
            onClick={() => setError(null)}
          >
            Dismiss
          </button>
        </div>
      )}

      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        <div className="bg-white p-6 rounded-lg shadow-md">
          <h2 className="text-xl font-semibold text-gray-800 mb-4">
            Generate New Molecule
          </h2>
          <form onSubmit={handleGenerate} className="space-y-4">
            <div>
              <label className="block text-sm font-medium text-gray-700">
                First SMILES String
              </label>
              <input
                type="text"
                value={smilesFirst}
                onChange={(e) => setSmilesFirst(e.target.value)}
                placeholder="e.g., CCO"
                disabled={loading}
                className="mt-1 w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700">
                Second SMILES String
              </label>
              <input
                type="text"
                value={smilesSecond}
                onChange={(e) => setSmilesSecond(e.target.value)}
                placeholder="e.g., c1ccccc1"
                disabled={loading}
                className="mt-1 w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100"
              />
            </div>
            <button
              type="submit"
              disabled={loading}
              className="w-full py-2 px-4 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-colors"
            >
              {loading ? 'Generating...' : 'Generate Molecule'}
            </button>
          </form>
        </div>

        <div className="bg-white p-6 rounded-lg shadow-md">
          <h2 className="text-xl font-semibold text-gray-800 mb-4">
            Fetch Generated Molecule
          </h2>
          <div className="space-y-4">
            <div>
              <label className="block text-sm font-medium text-gray-700">
                Molecule ID
              </label>
              <input
                type="text"
                value={fetchMoleculeId}
                onChange={(e) => setFetchMoleculeId(e.target.value)}
                placeholder="Enter molecule ID"
                disabled={loading}
                className="mt-1 w-full p-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100"
              />
            </div>
            <button
              onClick={handleFetch}
              disabled={loading}
              className="w-full py-2 px-4 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-colors"
            >
              {loading ? 'Fetching...' : 'Fetch Molecule'}
            </button>
          </div>
        </div>
      </div>

      {newSmiles && newIupacName && (
        <div className="mt-8 bg-white p-6 rounded-lg shadow-md">
          <h2 className="text-xl font-semibold text-gray-800 mb-4">
            Generated Molecule
          </h2>
          <p className="text-gray-700">
            <strong>SMILES:</strong> {newSmiles}
          </p>
          <p className="mt-2 text-gray-700">
            <strong>IUPAC Name:</strong> {newIupacName}
          </p>
          {moleculeId && (
            <p className="mt-2 text-gray-700">
              <strong>Molecule ID:</strong> {moleculeId}
            </p>
          )}
        </div>
      )}

      {loading && (
        <div className="flex items-center justify-center mt-6">
          <Spinner />
        </div>
      )}
    </div>
  );
};

export default ProteinStructureEvolution;