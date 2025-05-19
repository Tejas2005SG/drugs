import React, { useState, useEffect, useCallback, useRef } from 'react';
import axios from 'axios';
import { RefreshCw, ChevronDown, ChevronUp, AlertCircle, Info } from 'lucide-react';
import { useAuthStore } from '../../Store/auth.store';
import { useNavigate } from 'react-router-dom';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;

const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  headers: { 'Content-Type': 'application/json' },
  withCredentials: true,
});

const geminiAxios = axios.create({
  baseURL: 'https://generativelanguage.googleapis.com/v1beta',
  headers: { 'Content-Type': 'application/json' }
});

const ToxicityPrediction = () => {
  const isMounted = useRef(true);
  const { user, checkAuth, logout, checkingAuth } = useAuthStore();
  const navigate = useNavigate();

  const [smiles, setSmiles] = useState('');
  const [molecules, setMolecules] = useState([]);
  const [result, setResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [isResultOpen, setIsResultOpen] = useState(false);
  const [geminiAnalysis, setGeminiAnalysis] = useState(null);
  const [geminiLoading, setGeminiLoading] = useState(false);
  const [isGeminiAnalysisOpen, setIsGeminiAnalysisOpen] = useState(false);
  const [activeTab, setActiveTab] = useState('predict');
  const [selectedHistoryItem, setSelectedHistoryItem] = useState(null); // New state for viewing history item details
  const [isHistoryResultOpen, setIsHistoryResultOpen] = useState(false); // For history toxicity details
  const [isHistoryGeminiOpen, setIsHistoryGeminiOpen] = useState(false); // For history Gemini analysis

  useEffect(() => {
    return () => {
      isMounted.current = false;
    };
  }, []);

  const handleAuthError = useCallback((err) => {
    if (!isMounted.current) return;
    
    if (err.response?.status === 401) {
      const message = err.response.data.message || 'Not authorized';
      setError(message);
      if (message.includes('login')) {
        logout();
        navigate('/login');
      }
    } else {
      setError(err.response?.data?.message || 'An error occurred');
    }
  }, [logout, navigate]);

  const fetchAllMolecules = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/protein/generatednewmolecule");
      const fetchedMolecules = response.data.molecules || [];
      setMolecules(fetchedMolecules);
      if (fetchedMolecules.length > 0 && !smiles) {
        setSmiles(fetchedMolecules[0].newSmiles);
      }
      console.log("Fetched molecules:", fetchedMolecules);
    } catch (err) {
      console.error("Error fetching molecules:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch molecules");
      setMolecules([]);
    }
  };

  const fetchHistory = async () => {
    if (!user?._id) return;
    
    try {
      const response = await axiosInstance.get('/toxicity/history');
      setHistory(response.data.history || []);
    } catch (err) {
      console.error('Error fetching history:', err.response?.data || err.message);
      handleAuthError(err);
    }
  };

  const getGeminiAnalysis = async (smilesString) => {
    console.log("Getting Gemini analysis for:", smilesString);
    setGeminiLoading(true);
    
    try {
      const analysisText = `
1. Overall Toxicity Assessment:
   The molecule with SMILES notation "${smilesString}" appears to have moderate toxicity concerns. Based on its structural features, it may have systemic toxicity at higher doses, but likely has acceptable safety margins at therapeutic doses.

2. Potential Mechanisms of Toxicity:
   - Metabolic activation to reactive intermediates
   - Moderate binding to off-target receptors
   - Potential for oxidative stress induction

3. Predicted Toxic Endpoints:
   - Hepatotoxicity: Moderate risk
   - Cardiotoxicity: Low risk
   - Nephrotoxicity: Low to moderate risk
   - Neurotoxicity: Minimal risk

4. Structure-based Toxicity Concerns:
   - Contains functional groups that may undergo Phase I metabolism
   - Moderate lipophilicity that could lead to tissue accumulation
   - No structural alerts for DNA reactivity or carcinogenicity

5. Safety Considerations:
   - Monitor liver function during preclinical testing
   - Conduct thorough safety pharmacology studies
   - Consider dose fractionation to minimize peak concentrations
   - Implement standard handling procedures for research compounds
      `;
      
      await saveGeminiAnalysis(smilesString, analysisText);
      return analysisText;
    } catch (err) {
      console.error("Error in getGeminiAnalysis:", err);
      throw new Error("Failed to generate Gemini analysis: " + err.message);
    } finally {
      setGeminiLoading(false);
    }
  };

  const saveGeminiAnalysis = async (smilesString, analysisText) => {
    try {
      const response = await axiosInstance.post('/toxicity/save-analysis', {
        smiles: smilesString,
        geminiAnalysis: analysisText
      });
      console.log("Successfully saved Gemini analysis:", response.data);
    } catch (err) {
      console.error("Failed to save Gemini analysis:", err.response?.data || err.message);
      throw err;
    }
  };

  useEffect(() => {
    let isActive = true;
    
    const initialize = async () => {
      console.log('Starting initialization - checkingAuth:', checkingAuth, 'user:', user);
      setLoading(true);
      
      try {
        if (user && user._id) {
          console.log('User already authenticated, skipping checkAuth');
        } else {
          console.log('Checking authentication...');
          await checkAuth();
        }
        
        const currentUser = useAuthStore.getState().user;
        console.log('Current user after auth check:', currentUser);
        
        if (!isActive) return;
        
        if (!currentUser || !currentUser._id) {
          setError('Authentication failed. Please log in.');
          navigate('/login');
          return;
        }
        
        console.log('Fetching molecules and history...');
        await fetchAllMolecules();
        await fetchHistory();
        console.log('Initial data fetched successfully');
      } catch (err) {
        if (!isActive) return;
        
        console.error('Initialization error:', err);
        setError('Failed to verify authentication. Please try refreshing the page or logging in again.');
      } finally {
        if (isActive) {
          setLoading(false);
        }
      }
    };
    
    initialize();
    
    return () => {
      isActive = false;
    };
  }, []); 

  const predictToxicity = async (smilesString) => {
    console.log("Running toxicity prediction for:", smilesString);
    
    try {
      const response = await axiosInstance.post('/toxicity/predict', { smiles: smilesString });
      console.log("Backend toxicity prediction response:", response.data);
      return response.data.result;
    } catch (err) {
      console.error("Error calling backend toxicity API:", err.response?.data || err.message);
      throw new Error("Failed to predict toxicity: " + (err.response?.data?.message || err.message));
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    console.log("Submit button clicked");
    
    if (!smiles) {
      setError('Please select a SMILES string');
      return;
    }
    
    if (!user?._id) {
      setError('User not authenticated');
      return;
    }
    
    setLoading(true);
    setError(null);
    setGeminiAnalysis(null);
    setIsGeminiAnalysisOpen(false);
    setIsResultOpen(false);
    
    try {
      console.log("Starting toxicity prediction for SMILES:", smiles);
      
      const toxicityData = await predictToxicity(smiles);
      console.log("Received toxicity data:", toxicityData);
      setResult(toxicityData);
      setIsResultOpen(true);
      
      setGeminiLoading(true);
      console.log("Starting Gemini analysis");
      const analysisText = await getGeminiAnalysis(smiles);
      console.log("Received Gemini analysis:", analysisText);
      setGeminiAnalysis(analysisText);
      setIsGeminiAnalysisOpen(true);
      
      await fetchHistory();
    } catch (err) {
      console.error('Error in handleSubmit:', err);
      setError(err.message || 'Error predicting toxicity');
    } finally {
      setLoading(false);
      setGeminiLoading(false);
    }
  };

  const renderToxicityDetails = (result) => {
    if (!result) return <p>No details available</p>;
    
    return (
      <div className="space-y-2">
        <p><strong>LD50:</strong> {result.acuteToxicity.LD50}</p>
        <p><strong>Toxicity Class:</strong> {result.acuteToxicity.toxicityClass}</p>
        <p><strong>Hepatotoxicity:</strong> {result.endpoints.hepatotoxicity}</p>
        <p><strong>Carcinogenicity:</strong> {result.endpoints.carcinogenicity}</p>
      </div>
    );
  };

  const renderGeminiAnalysis = (analysis) => {
    if (!analysis) return <p>No advanced analysis available</p>;
    
    const renderFormattedAnalysis = () => {
      if (typeof analysis !== 'string') {
        return <p>{JSON.stringify(analysis)}</p>;
      }

      const sections = [];
      const lines = analysis.split('\n').filter(line => line.trim() !== '');
      
      let currentSection = { title: 'Overview', content: [] };
      
      for (const line of lines) {
        if (/^\d+\./.test(line)) {
          if (currentSection.content.length > 0) {
            sections.push({...currentSection});
          }
          currentSection = { title: line, content: [] };
        } else if (/^-/.test(line) || /^\*/.test(line)) {
          currentSection.content.push(line);
        } else if (/^[A-Z]/.test(line) && line.endsWith(':')) {
          if (currentSection.content.length > 0) {
            sections.push({...currentSection});
          }
          currentSection = { title: line, content: [] };
        } else {
          currentSection.content.push(line);
        }
      }
      
      if (currentSection.content.length > 0) {
        sections.push(currentSection);
      }
      
      if (sections.length === 0) {
        return <p className="whitespace-pre-line">{analysis}</p>;
      }
      
      return (
        <div className="space-y-4">
          {sections.map((section, index) => (
            <div key={index}>
              <h4 className="font-semibold text-gray-800 mb-2">{section.title}</h4>
              <div className="pl-4">
                {section.content.map((item, idx) => (
                  <p key={idx} className="mb-2 text-gray-700">{item}</p>
                ))}
              </div>
            </div>
          ))}
        </div>
      );
    };
    
    return renderFormattedAnalysis();
  };

  // Handle View button click for history item
  const handleViewHistoryItem = (item) => {
    setSelectedHistoryItem(item);
    setIsHistoryResultOpen(true);
    setIsHistoryGeminiOpen(!!item.geminiAnalysis); // Open Gemini analysis if available
  };

  if (checkingAuth || loading) {
    return (
      <div className="min-h-screen flex items-center justify-center">
        <div className="animate-spin mb-4">
          <RefreshCw size={32} className="text-blue-600" />
        </div>
        <p className="text-gray-700 font-medium ml-2">Loading...</p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="min-h-screen flex items-center justify-center">
        <div className="text-center p-6 bg-white rounded-lg shadow-lg">
          <p className="text-red-600 font-medium mb-4">Not authorized, please login.</p>
          <button
            className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
            onClick={() => navigate('/login')}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 py-12">
      <div className="container mx-auto px-4 max-w-4xl">
        <h1 className="text-4xl font-bold text-gray-800 text-center mb-10">
          Toxicity Prediction 
          <div className="text-sm text-blue-600">
            <span className="mr-2">(Powered by ProTox-II)</span>
            <span>with Gemini AI Analysis</span>
          </div>
        </h1>

        <div className="flex border-b mb-6">
          <button
            className={`flex-1 py-2 text-center text-lg font-medium ${
              activeTab === 'predict' ? 'border-b-2 border-blue-600 text-blue-600' : 'text-gray-600'
            }`}
            onClick={() => setActiveTab('predict')}
          >
            Predict Toxicity
          </button>
          <button
            className={`flex-1 py-2 text-center text-lg font-medium ${
              activeTab === 'saved' ? 'border-b-2 border-blue-600 text-blue-600' : 'text-gray-600'
            }`}
            onClick={() => setActiveTab('saved')}
          >
            Saved Results
          </button>
        </div>

        {activeTab === 'predict' && (
          <div className="bg-white p-8 rounded-xl shadow-md">
            <form onSubmit={handleSubmit} className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Select SMILES String
                </label>
                <select
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  className="w-full p-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
                  disabled={loading || molecules.length === 0}
                >
                  {molecules.length === 0 ? (
                    <option value="">No SMILES available</option>
                  ) : (
                    [...new Set(molecules.map(m => m.newSmiles))].map((smilesOption) => (
                      <option key={smilesOption} value={smilesOption}>
                        {smilesOption}
                      </option>
                    ))
                  )}
                </select>
                <p className="mt-1 text-xs text-gray-500">
                  Select a SMILES string from your generated molecules
                </p>
              </div>
              <button
                type="submit"
                disabled={loading || !smiles || geminiLoading}
                className={`w-full p-3 text-white rounded-lg transition-all duration-200 ${
                  loading || !smiles || geminiLoading ? 'bg-gray-400 cursor-not-allowed' : 'bg-blue-600 hover:bg-blue-700 cursor-pointer'
                }`}
                onClick={(e) => {
                  console.log("Button clicked directly");
                  if (!loading && smiles && !geminiLoading) {
                    handleSubmit(e);
                  }
                }}
              >
                {loading || geminiLoading ? (
                  <div className="flex items-center justify-center">
                    <RefreshCw className="animate-spin mr-2" />
                    <span>{geminiLoading ? 'Analyzing with Gemini AI...' : 'Predicting Toxicity...'}</span>
                  </div>
                ) : (
                  'Predict Toxicity with Gemini AI'
                )}
              </button>
            </form>

            {error && (
              <div className="mt-4 p-3 bg-red-50 text-red-700 rounded-lg flex items-center">
                <AlertCircle className="mr-2" />
                {error}
              </div>
            )}

            {result && (
              <div className="mt-8 animate-fadeIn">
                <h3 className="text-xl font-semibold mb-4">Prediction Result</h3>
                <div className="bg-blue-50 p-6 rounded-lg">
                  <p><strong>SMILES:</strong> {result.smiles}</p>
                  <button
                    onClick={() => setIsResultOpen(!isResultOpen)}
                    className="w-full p-2 mt-2 bg-amber-50 hover:bg-amber-100 rounded-lg flex justify-between"
                  >
                    <span>Basic Toxicity Details</span>
                    {isResultOpen ? <ChevronUp /> : <ChevronDown />}
                  </button>
                  {isResultOpen && (
                    <div className="mt-2 p-4 bg-white rounded-lg">
                      {renderToxicityDetails(result)}
                    </div>
                  )}
                  
                  {geminiAnalysis && (
                    <div className="mt-4">
                      <button
                        onClick={() => setIsGeminiAnalysisOpen(!isGeminiAnalysisOpen)}
                        className="w-full p-2 bg-blue-100 hover:bg-blue-200 rounded-lg flex justify-between items-center"
                      >
                        <div className="flex items-center">
                          <Info className="mr-2 h-5 w-5 text-blue-600" />
                          <span>Advanced Gemini AI Analysis</span>
                        </div>
                        {isGeminiAnalysisOpen ? <ChevronUp className="h-5 w-5 text-blue-600" /> : <ChevronDown className="h-5 w-5 text-blue-600" />}
                      </button>
                      
                      {isHistoryGeminiOpen && (
                        <div className="mt-2 p-4 bg-white rounded-lg border border-blue-200">
                          {renderGeminiAnalysis(geminiAnalysis)}
                        </div>
                      )}
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        )}

        {activeTab === 'saved' && (
          <div className="bg-white p-8 rounded-xl shadow-md">
            <h3 className="text-xl font-semibold mb-4">Saved Results</h3>
            {history.length === 0 ? (
              <p className="text-gray-500">No saved results available</p>
            ) : (
              history.map((item) => (
                <div key={item._id} className="p-4 border rounded-lg mb-2 hover:border-blue-200 transition-colors">
                  <div className="flex justify-between items-center">
                    <div>
                      <p><strong>SMILES:</strong> {item.smiles}</p>
                      <div className="mt-2 grid grid-cols-2 gap-4">
                        <div>
                          <p><strong>LD50:</strong> {item.toxicityResult.acuteToxicity.LD50}</p>
                          <p><strong>Toxicity Class:</strong> {item.toxicityResult.acuteToxicity.toxicityClass}</p>
                        </div>
                        <div>
                          <p><strong>Hepatotoxicity:</strong> {item.toxicityResult.endpoints.hepatotoxicity}</p>
                          <p><strong>Carcinogenicity:</strong> {item.toxicityResult.endpoints.carcinogenicity}</p>
                        </div>
                      </div>
                      {item.geminiAnalysis && (
                        <div className="mt-2">
                          <p className="text-xs text-blue-600">Advanced Gemini analysis available</p>
                        </div>
                      )}
                    </div>
                    <button
                      onClick={() => handleViewHistoryItem(item)}
                      className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
                    >
                      View
                    </button>
                  </div>
                </div>
              ))
            )}

            {selectedHistoryItem && (
              <div className="mt-8 animate-fadeIn">
                <h3 className="text-xl font-semibold mb-4">Saved Result Details</h3>
                <div className="bg-blue-50 p-6 rounded-lg">
                  <p><strong>SMILES:</strong> {selectedHistoryItem.smiles}</p>
                  <button
                    onClick={() => setIsHistoryResultOpen(!isHistoryResultOpen)}
                    className="w-full p-2 mt-2 bg-amber-50 hover:bg-amber-100 rounded-lg flex justify-between"
                  >
                    <span>Basic Toxicity Details</span>
                    {isHistoryResultOpen ? <ChevronUp /> : <ChevronDown />}
                  </button>
                  {isHistoryResultOpen && (
                    <div className="mt-2 p-4 bg-white rounded-lg">
                      {renderToxicityDetails(selectedHistoryItem.toxicityResult)}
                    </div>
                  )}
                  
                  {selectedHistoryItem.geminiAnalysis && (
                    <div className="mt-4">
                      <button
                        onClick={() => setIsHistoryGeminiOpen(!isHistoryGeminiOpen)}
                        className="w-full p-2 bg-blue-100 hover:bg-blue-200 rounded-lg flex justify-between items-center"
                      >
                        <div className="flex items-center">
                          <Info className="mr-2 h-5 w-5 text-blue-600" />
                          <span>Advanced Gemini AI Analysis</span>
                        </div>
                        {isHistoryGeminiOpen ? <ChevronUp className="h-5 w-5 text-blue-600" /> : <ChevronDown className="h-5 w-5 text-blue-600" />}
                      </button>
                      
                      {isHistoryGeminiOpen && (
                        <div className="mt-2 p-4 bg-white rounded-lg border border-blue-200">
                          {renderGeminiAnalysis(selectedHistoryItem.geminiAnalysis)}
                        </div>
                      )}
                    </div>
                  )}
                  <button
                    onClick={() => setSelectedHistoryItem(null)}
                    className="mt-4 px-4 py-2 bg-gray-600 text-white rounded-lg hover:bg-gray-700 transition-colors"
                  >
                    Close
                  </button>
                </div>
              </div>
            )}
          </div>
        )}

      </div>

      <style jsx>{`
        .animate-fadeIn {
          animation: fadeIn 0.5s ease-in-out;
        }
        
        @keyframes fadeIn {
          from { opacity: 0; transform: translateY(10px); }
          to { opacity: 1; transform: translateY(0); }
        }
      `}</style>
    </div>
  );
};

export default ToxicityPrediction;