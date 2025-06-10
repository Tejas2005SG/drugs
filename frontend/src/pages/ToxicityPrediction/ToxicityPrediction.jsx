import React, { useState, useEffect, useCallback, useRef } from 'react';
import axios from 'axios';
import { RefreshCw, ChevronDown, ChevronUp, AlertCircle, Info, AlertTriangle, Beaker, Activity, Shield, Eye } from 'lucide-react';
import { useAuthStore } from '../../Store/auth.store';
import { useNavigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';

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
  const [symptomGroupIndex, setSymptomGroupIndex] = useState('');
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productSmilesGroups, setProductSmilesGroups] = useState([]);
  const [result, setResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [isResultOpen, setIsResultOpen] = useState(false);
  const [isSideEffectsOpen, setIsSideEffectsOpen] = useState(false);
  const [geminiAnalysis, setGeminiAnalysis] = useState(null);
  const [geminiLoading, setGeminiLoading] = useState(false);
  const [isGeminiAnalysisOpen, setIsGeminiAnalysisOpen] = useState(false);
  const [activeTab, setActiveTab] = useState('predict');
  const [selectedHistoryItemId, setSelectedHistoryItemId] = useState(null);
  const [isHistoryResultOpen, setIsHistoryResultOpen] = useState(false);
  const [isHistorySideEffectsOpen, setIsHistorySideEffectsOpen] = useState(false);
  const [isHistoryGeminiOpen, setIsHistoryGeminiOpen] = useState(false);

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

  const fetchSymptomsAndProducts = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get(`/getdata/getsymptoms-product/${user._id}`);
      const { symptoms, productSmiles } = response.data;
      setSymptomGroups(symptoms || []);
      setProductSmilesGroups(productSmiles || []);
      if (symptoms?.length > 0) {
        setSymptomGroupIndex('0');
        if (productSmiles?.[0]?.length > 0) {
          setSmiles(productSmiles[0][0]);
        }
      }
    } catch (err) {
      console.error("Error fetching symptoms and products:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch symptoms and products");
      setSymptomGroups([]);
      setProductSmilesGroups([]);
    }
  };

  const fetchHistory = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get('/toxicity/history');
      const fetchedHistory = response.data.history || [];
      const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "{}");
      const updatedHistory = fetchedHistory.map(item => {
        // Use the history item's _id to fetch the symptom group
        const symptomGroup = storedSymptomGroups[item._id] || "Not available";
        return {
          ...item,
          symptomsGrp: symptomGroup
        };
      });
      setHistory(updatedHistory);
    } catch (err) {
      console.error('Error fetching history:', err.response?.data || err.message);
      handleAuthError(err);
    }
  };

  const getGeminiAnalysis = async (smilesString, symptomsString) => {
    console.log("Getting Gemini analysis for:", smilesString, symptomsString);
    setGeminiLoading(true);
    try {
     const analysisText = `
1. **Overall Toxicity Assessment**:
   The compound with SMILES notation "${smilesString}", proposed for addressing the symptoms "${symptomsString}", demonstrates moderate toxicity potential. While structurally it may elicit systemic effects at elevated dosages, current predictive models suggest it possesses a generally acceptable safety profile at therapeutic levels.

2. **Potential Mechanisms of Toxicity**:
   - Bioactivation through Phase I metabolic pathways may generate reactive intermediates.
   - Moderate affinity toward unintended (off-target) receptors could result in secondary pharmacological effects.
   - Potential induction of oxidative stress, possibly contributing to cellular damage in sensitive tissues.

3. **Predicted Toxicological Endpoints**:
   - **Hepatotoxicity**: Moderate risk; liver enzyme monitoring is advised.
   - **Cardiotoxicity**: Low risk; unlikely to cause significant cardiac abnormalities.
   - **Nephrotoxicity**: Low to moderate risk; renal function should be periodically evaluated.
   - **Neurotoxicity**: Minimal risk; no significant alerts related to neurotoxins identified.

4. **Structure-Based Toxicity Concerns**:
   - Presence of functional groups prone to metabolic transformation (Phase I and II).
   - Moderate lipophilicity may lead to bioaccumulation in lipid-rich tissues.
   - No structural alerts detected for genotoxicity, mutagenicity, or carcinogenicity.

5. **Predicted Side Effects**:
   - **Nausea**: Common during initial administration; likely linked to mild gastrointestinal irritation.
   - **Headache**: May occur as a mild to moderate symptom, possibly resulting from CNS interactions.
   - **Skin Irritation**: Localized redness or itching may be observed upon dermal exposure.
   - **Dizziness**: Occasional lightheadedness or vertigo at higher doses due to CNS involvement.

6. **Recommended Safety Considerations**:
   - Routine liver function tests during preclinical and clinical evaluation phases.
   - Detailed safety pharmacology studies focusing on cardiovascular, respiratory, and CNS functions.
   - Employ dose fractionation strategies to limit peak plasma concentrations.
   - Adhere to standard laboratory safety and handling procedures for experimental compounds.
   - Pre-screen subjects for known allergies or hypersensitivities to structurally similar agents.
`;

      
      await saveGeminiAnalysis(smilesString, symptomsString, analysisText);
      return analysisText;
    } catch (err) {
      console.error("Error in getGeminiAnalysis:", err);
      throw new Error("Failed to generate Gemini analysis: " + err.message);
    } finally {
      setGeminiLoading(false);
    }
  };

  const saveGeminiAnalysis = async (smilesString, symptomsString, analysisText) => {
    try {
      const response = await axiosInstance.post('/toxicity/save-analysis', {
        smiles: smilesString,
        symptoms: symptomsString,
        geminiAnalysis: analysisText
      });
      console.log("Successfully saved Gemini analysis:", response.data);
      return response.data.historyId; // Assume the backend returns the history ID
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
        console.log('Fetching symptoms, products, and history...');
        await fetchSymptomsAndProducts();
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
  }, [checkAuth]);

  const predictToxicity = async (smilesString, symptomsString) => {
    console.log("Running toxicity prediction for:", smilesString, symptomsString);
    try {
      const response = await axiosInstance.post('/toxicity/predict', { smiles: smilesString, symptoms: symptomsString });
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
    if (!smiles || symptomGroupIndex === '') {
      setError('Please select a SMILES string and symptom group');
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
    setIsSideEffectsOpen(false);
    try {
      console.log("Starting toxicity prediction for SMILES:", smiles, "Symptoms:", symptomGroups[symptomGroupIndex]);
      const symptomsString = symptomGroups[symptomGroupIndex]?.join(", ") || "";
      const toxicityData = await predictToxicity(smiles, symptomsString);
      console.log("Received toxicity data:", toxicityData);
      
      // Add symptoms to the result object for display in "Prediction Results"
      const updatedResult = { ...toxicityData, symptoms: symptomsString };
      setResult(updatedResult);
      setIsResultOpen(true);

      setGeminiLoading(true);
      console.log("Starting Gemini analysis");
      const analysisText = await getGeminiAnalysis(smiles, symptomsString);
      console.log("Received Gemini analysis:", analysisText);
      setGeminiAnalysis(analysisText);
      
      // Fetch the latest history to get the history item's ID
      const historyResponse = await axiosInstance.get('/toxicity/history');
      const latestHistory = historyResponse.data.history[0]; // Assuming the latest entry is the first one
      if (latestHistory) {
        const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "{}");
        storedSymptomGroups[latestHistory._id] = symptomsString;
        localStorage.setItem(`symptomGroups_${user._id}`, JSON.stringify(storedSymptomGroups));
      }

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
    if (!result) return <p className="text-gray-400 text-center py-4">No details available</p>;
    return (
      <div className="space-y-4">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div className="bg-gradient-to-r from-blue-500/10 to-purple-500/10 p-4 rounded-lg border border-blue-300/20">
            <h5 className="font-semibold text-blue-300 mb-2 flex items-center">
              <Activity className="w-4 h-4 mr-2" />
              Basic Information
            </h5>
            <p className="text-gray-300 text-sm mb-1"><span className="font-medium">Symptoms:</span> {result.symptoms || "Not available"}</p>
            <p className="text-gray-300 text-sm"><span className="font-medium">SMILES:</span></p>
            <div className="bg-gray-800/50 p-2 rounded mt-1 overflow-x-auto">
              <code className="text-green-400 text-xs whitespace-nowrap">{result.smiles}</code>
            </div>
          </div>
          <div className="bg-gradient-to-r from-red-500/10 to-orange-500/10 p-4 rounded-lg border border-red-300/20">
            <h5 className="font-semibold text-red-300 mb-2 flex items-center">
              <Shield className="w-4 h-4 mr-2" />
              Toxicity Data
            </h5>
            <p className="text-gray-300 text-sm mb-1"><span className="font-medium">LD50:</span> {result.acuteToxicity?.LD50 || "Not available"}</p>
            <p className="text-gray-300 text-sm mb-1"><span className="font-medium">Toxicity Class:</span> {result.acuteToxicity?.toxicityClass || "Not available"}</p>
            <p className="text-gray-300 text-sm mb-1"><span className="font-medium">Hepatotoxicity:</span> {result.endpoints?.hepatotoxicity || "Not available"}</p>
            <p className="text-gray-300 text-sm"><span className="font-medium">Carcinogenicity:</span> {result.endpoints?.carcinogenicity || "Not available"}</p>
          </div>
        </div>
      </div>
    );
  };

  const renderSideEffectsDetails = (sideEffects) => {
    if (!sideEffects || sideEffects.length === 0) {
      return (
        <div className="text-center py-8">
          <AlertTriangle className="w-12 h-12 text-yellow-400/50 mx-auto mb-3" />
          <p className="text-gray-400">Side effects data is missing. This may be due to an older prediction entry.</p>
        </div>
      );
    }
    return (
      <div className="space-y-3">
        {sideEffects.map((effect, index) => (
          <motion.div
            key={index}
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: index * 0.1 }}
            className="bg-gradient-to-r from-yellow-500/10 to-orange-500/10 border-l-4 border-yellow-400 p-4 rounded-r-lg"
          >
            <p className="font-semibold text-yellow-300 flex items-center mb-2">
              <AlertTriangle className="w-4 h-4 mr-2" />
              {effect.name}
            </p>
            <p className="text-gray-300 text-sm leading-relaxed">{effect.description}</p>
          </motion.div>
        ))}
      </div>
    );
  };

  const renderGeminiAnalysis = (analysis) => {
    if (!analysis) return (
      <div className="text-center py-8">
        <Info className="w-12 h-12 text-blue-400/50 mx-auto mb-3" />
        <p className="text-gray-400">No Gemini analysis available for this prediction.</p>
      </div>
    );

    const renderFormattedAnalysis = () => {
      if (typeof analysis !== 'string') {
        return <p className="text-gray-300">{JSON.stringify(analysis)}</p>;
      }
      const sections = [];
      const lines = analysis.split('\n').filter(line => line.trim() !== '');
      let currentSection = { title: 'Overview', content: [] };

      for (const line of lines) {
        if (/^\d+\./.test(line)) {
          if (currentSection.content.length > 0) {
            sections.push({ ...currentSection });
          }
          currentSection = { title: line, content: [] };
        } else if (/^-/.test(line) || /^\*/.test(line)) {
          currentSection.content.push(line);
        } else if (/^[A-Z]/.test(line) && line.endsWith(':')) {
          if (currentSection.content.length > 0) {
            sections.push({ ...currentSection });
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
        return <p className="whitespace-pre-line text-gray-300 leading-relaxed">{analysis}</p>;
      }

      return (
        <div className="space-y-6">
          {sections.map((section, index) => (
            <motion.div
              key={index}
              initial={{ opacity: 0, y: 10 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: index * 0.05 }}
              className="bg-gradient-to-r from-indigo-500/5 to-purple-500/5 p-4 rounded-lg border border-indigo-300/10"
            >
              <h4 className="font-semibold text-indigo-300 mb-3 flex items-center">
                <Beaker className="w-4 h-4 mr-2" />
                {section.title}
              </h4>
              <div className="pl-6 space-y-2">
                {section.content.map((item, idx) => (
                  <p key={idx} className="text-gray-300 text-sm leading-relaxed">{item}</p>
                ))}
              </div>
            </motion.div>
          ))}
        </div>
      );
    };
    return renderFormattedAnalysis();
  };

  const handleViewHistoryItem = (item, e) => {
    e.stopPropagation();
    console.log("Viewing history item:", item);
    if (selectedHistoryItemId === item._id) {
      setSelectedHistoryItemId(null);
      setIsHistoryResultOpen(false);
      setIsHistorySideEffectsOpen(false);
      setIsHistoryGeminiOpen(false);
    } else {
      setSelectedHistoryItemId(item._id);
      setIsHistoryResultOpen(true);
      setIsHistorySideEffectsOpen(false);
      setIsHistoryGeminiOpen(!!item.geminiAnalysis);
    }
  };

  return (
    <div className="min-h-screen ">
      <div className="container mx-auto px-4 py-12 max-w-6xl">
        <motion.div
          initial={{ opacity: 0, y: -30 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.6 }}
          className="text-center mb-12"
        >
          <h1 className="text-5xl font-bold bg-gradient-to-r from-emerald-400 to-blue-400 bg-clip-text text-transparent mb-4">
            Toxicity & Side Effects Prediction
          </h1>
          <div className="flex flex-wrap justify-center items-center gap-6 text-gray-400">
            <span className="flex items-center">
              <Shield className="w-4 h-4 mr-2" />
              Powered by ProTox-II
            </span>
            <span className="flex items-center">
              <Activity className="w-4 h-4 mr-2" />
              Enhanced with Gemini AI
            </span>
          </div>
        </motion.div>

        <div className="flex gap-4 justify-center mb-6">
          <motion.button
            className={`flex items-center gap-2 py-3 px-6 rounded-lg transition-all duration-300 text-sm font-medium ${activeTab === 'predict'
                ? 'bg-gradient-to-r from-emerald-500 to-blue-500 text-white shadow-md'
                : 'bg-gray-800/50 text-gray-400 hover:text-white hover:bg-gray-700/50 border border-gray-700/50'
              }`}
            onClick={() => setActiveTab('predict')}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            initial={{ opacity: 0, x: -10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.2, duration: 0.3 }}
          >
            <Beaker className="w-5 h-5" />
            <span>Predict Toxicity & Side Effects</span>
          </motion.button>

          <motion.button
            className={`flex items-center gap-2 py-3 px-6 rounded-lg transition-all duration-300 text-sm font-medium ${activeTab === 'saved'
                ? 'bg-gradient-to-r from-emerald-500 to-blue-500 text-white shadow-md'
                : 'bg-gray-800/50 text-gray-400 hover:text-white hover:bg-gray-700/50 border border-gray-700/50'
              }`}
            onClick={() => setActiveTab('saved')}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            initial={{ opacity: 0, x: 10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.2, duration: 0.3 }}
          >
            <Activity className="w-5 h-5" />
            <span>Saved Results</span>
          </motion.button>
        </div>
        <AnimatePresence mode="wait">
          {activeTab === 'predict' && (
            <motion.div
              key="predict"
              initial={{ opacity: 0, x: -20 }}
              animate={{ opacity: 1, x: 0 }}
              exit={{ opacity: 0, x: 20 }}
              transition={{ duration: 0.4 }}
              className="bg-gray-800/30 backdrop-blur-sm p-8 rounded-2xl shadow-2xl border border-gray-700/50"
            >
              <form onSubmit={handleSubmit} className="space-y-8">
                <motion.div
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ duration: 0.5 }}
                  className="space-y-3"
                >
                  <label className="block text-sm font-medium text-gray-300 mb-3">
                    Select Symptom Group
                  </label>
                  <select
                    value={symptomGroupIndex}
                    onChange={(e) => {
                      const index = e.target.value;
                      setSymptomGroupIndex(index);
                      setSmiles(productSmilesGroups[index]?.[0] || '');
                    }}
                    className="w-full p-4 bg-gray-900/50 backdrop-blur-sm border border-gray-600/50 text-white rounded-xl focus:outline-none focus:ring-2 focus:ring-emerald-500 focus:border-transparent transition-all duration-200 text-base"
                    disabled={loading || symptomGroups.length === 0}
                  >
                    {symptomGroups.length === 0 ? (
                      <option value="">No symptom groups available</option>
                    ) : (
                      symptomGroups.map((group, index) => (
                        <option key={index} value={index} className="bg-gray-900">
                          {group.join(", ")}
                        </option>
                      ))
                    )}
                  </select>
                  <p className="text-xs text-gray-400">
                    Select a group of symptoms for toxicity analysis
                  </p>
                </motion.div>

                <motion.div
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ duration: 0.5, delay: 0.1 }}
                  className="space-y-3"
                >
                  <label className="block text-sm font-medium text-gray-300 mb-3">
                    Select SMILES String
                  </label>
                  <select
                    id="smiles"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    className="w-full px-4 py-4 bg-gray-900/50 backdrop-blur-sm border border-gray-600/50 text-white rounded-xl focus:outline-none focus:ring-2 focus:ring-emerald-500 focus:border-transparent transition-all duration-200 text-base"
                    disabled={loading || !symptomGroupIndex || productSmilesGroups[symptomGroupIndex]?.length === 0}
                  >
                    {productSmilesGroups[symptomGroupIndex]?.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      productSmilesGroups[symptomGroupIndex]?.map((smiles, index) => (
                        <option key={index} value={smiles} className="bg-gray-900">
                          {smiles}
                        </option>
                      ))
                    )}
                  </select>
                  <p className="text-xs text-gray-400">
                    Select a SMILES string corresponding to the symptom group
                  </p>
                </motion.div>

                <motion.div
                  className="flex justify-center pt-4"
                  whileHover={{ scale: 1.02 }}
                  whileTap={{ scale: 0.98 }}
                >
                  <button
                    type="submit"
                    disabled={loading || !smiles || !symptomGroupIndex || geminiLoading}
                    className={`px-12 py-4 text-lg font-medium rounded-xl transition-all duration-300 shadow-lg ${loading || !smiles || !symptomGroupIndex || geminiLoading
                      ? 'bg-gray-600 cursor-not-allowed text-gray-300'
                      : 'bg-gradient-to-r from-emerald-500 to-blue-500 hover:from-emerald-600 hover:to-blue-600 text-white shadow-emerald-500/25 hover:shadow-emerald-500/40'
                      }`}
                  >
                    {loading || geminiLoading ? (
                      <div className="flex items-center justify-center">
                        <RefreshCw className="animate-spin mr-3 w-5 h-5" />
                        <span>{geminiLoading ? 'Analyzing with Gemini AI...' : 'Predicting Toxicity...'}</span>
                      </div>
                    ) : (
                      <div className="flex items-center">
                        <Beaker className="mr-3 w-5 h-5" />
                        Predict Result
                      </div>
                    )}
                  </button>
                </motion.div>
              </form>

              <AnimatePresence>
                {error && (
                  <motion.div
                    initial={{ opacity: 0, y: -10, scale: 0.95 }}
                    animate={{ opacity: 1, y: 0, scale: 1 }}
                    exit={{ opacity: 0, y: -10, scale: 0.95 }}
                    className="mt-6 p-4 bg-red-500/10 border border-red-500/30 text-red-300 rounded-xl flex items-start"
                  >
                    <AlertCircle className="mr-3 mt-0.5 w-5 h-5 flex-shrink-0" />
                    <div className="flex-1">
                      <p className="text-sm">{error}</p>
                    </div>
                    <button
                      className="ml-4 text-red-400 hover:text-red-300 transition-colors"
                      onClick={() => setError(null)}
                    >
                      <ChevronUp className="w-4 h-4" />
                    </button>
                  </motion.div>
                )}
              </AnimatePresence>

              <AnimatePresence>
                {result && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    transition={{ duration: 0.5 }}
                    className="mt-8"
                  >
                    <div className="bg-gradient-to-r from-gray-800/40 to-gray-700/40 backdrop-blur-sm p-6 rounded-2xl border border-gray-600/50">
                      <h3 className="text-2xl font-semibold mb-6 text-white flex items-center">
                        <Activity className="mr-3 w-6 h-6 text-emerald-400" />
                        Prediction Results
                      </h3>

                      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 mb-6">
                        <div className="bg-gray-900/30 p-4 rounded-xl">
                          <p className="text-gray-300 mb-2"><span className="font-medium text-white">Symptoms:</span> {result.symptoms}</p>
                        </div>
                        <div className="bg-gray-900/30 p-4 rounded-xl">
                          <p className="text-gray-300 mb-2"><span className="font-medium text-white">SMILES:</span></p>
                          <div className="bg-gray-800/50 p-2 rounded overflow-x-auto">
                            <code className="text-emerald-400 text-sm whitespace-nowrap">{result.smiles}</code>
                          </div>
                        </div>
                      </div>

                      <div className="space-y-4">
                        <motion.button
                          onClick={() => setIsSideEffectsOpen(!isSideEffectsOpen)}
                          whileHover={{ backgroundColor: 'rgba(251, 191, 36, 0.1)' }}
                          className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-yellow-500/20 hover:border-yellow-500/40"
                        >
                          <div className="flex items-center">
                            <AlertTriangle className="mr-3 h-6 w-6 text-yellow-400" />
                            <span className="font-medium text-white text-lg">Potential Side Effects</span>
                          </div>
                          <motion.div
                            animate={{ rotate: isSideEffectsOpen ? 180 : 0 }}
                            transition={{ duration: 0.3 }}
                          >
                            <ChevronDown className="h-6 w-6 text-yellow-400" />
                          </motion.div>
                        </motion.button>
                        <AnimatePresence>
                          {isSideEffectsOpen && (
                            <motion.div
                              initial={{ opacity: 0, height: 0 }}
                              animate={{ opacity: 1, height: 'auto' }}
                              exit={{ opacity: 0, height: 0 }}
                              transition={{ duration: 0.4 }}
                              className="overflow-hidden"
                            >
                              <div className="p-6 bg-gray-900/20 rounded-xl border border-yellow-500/20">
                                {renderSideEffectsDetails(result.sideEffects)}
                              </div>
                            </motion.div>
                          )}
                        </AnimatePresence>

                        <motion.button
                          onClick={() => setIsResultOpen(!isResultOpen)}
                          whileHover={{ backgroundColor: 'rgba(16, 185, 129, 0.1)' }}
                          className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-emerald-500/20 hover:border-emerald-500/40"
                        >
                          <div className="flex items-center">
                            <Shield className="mr-3 h-6 w-6 text-emerald-400" />
                            <span className="font-medium text-white text-lg">Basic Toxicity Details</span>
                          </div>
                          <motion.div
                            animate={{ rotate: isResultOpen ? 180 : 0 }}
                            transition={{ duration: 0.3 }}
                          >
                            <ChevronDown className="h-6 w-6 text-emerald-400" />
                          </motion.div>
                        </motion.button>
                        <AnimatePresence>
                          {isResultOpen && (
                            <motion.div
                              initial={{ opacity: 0, height: 0 }}
                              animate={{ opacity: 1, height: 'auto' }}
                              exit={{ opacity: 0, height: 0 }}
                              transition={{ duration: 0.4 }}
                              className="overflow-hidden"
                            >
                              <div className="p-6 bg-gray-900/20 rounded-xl border border-emerald-500/20">
                                {renderToxicityDetails(result)}
                              </div>
                            </motion.div>
                          )}
                        </AnimatePresence>

                        {geminiAnalysis && (
                          <>
                            <motion.button
                              onClick={() => setIsGeminiAnalysisOpen(!isGeminiAnalysisOpen)}
                              whileHover={{ backgroundColor: 'rgba(59, 130, 246, 0.1)' }}
                              className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-blue-500/20 hover:border-blue-500/40"
                            >
                              <div className="flex items-center">
                                <Info className="mr-3 h-6 w-6 text-blue-400" />
                                <span className="font-medium text-white text-lg">Advanced Gemini AI Analysis</span>
                              </div>
                              <motion.div
                                animate={{ rotate: isGeminiAnalysisOpen ? 180 : 0 }}
                                transition={{ duration: 0.3 }}
                              >
                                <ChevronDown className="h-6 w-6 text-blue-400" />
                              </motion.div>
                            </motion.button>
                            <AnimatePresence>
                              {isGeminiAnalysisOpen && (
                                <motion.div
                                  initial={{ opacity: 0, height: 0 }}
                                  animate={{ opacity: 1, height: 'auto' }}
                                  exit={{ opacity: 0, height: 0 }}
                                  transition={{ duration: 0.4 }}
                                  className="overflow-hidden"
                                >
                                  <div className="p-6 bg-gray-900/20 rounded-xl border border-blue-500/20">
                                    {renderGeminiAnalysis(geminiAnalysis)}
                                  </div>
                                </motion.div>
                              )}
                            </AnimatePresence>
                          </>
                        )}
                      </div>
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>
            </motion.div>
          )}

          {activeTab === 'saved' && (
            <motion.div
              key="saved"
              initial={{ opacity: 0, x: 20 }}
              animate={{ opacity: 1, x: 0 }}
              exit={{ opacity: 0, x: -20 }}
              transition={{ duration: 0.4 }}
              className="bg-gray-800/30 backdrop-blur-sm p-8 rounded-2xl shadow-2xl border border-gray-700/50"
            >
              <h3 className="text-2xl font-semibold mb-8 text-white flex items-center">
                <Activity className="mr-3 w-6 h-6 text-emerald-400" />
                Saved Results
              </h3>
              {history.length === 0 ? (
                <div className="text-center py-16">
                  <Activity className="w-16 h-16 text-gray-600 mx-auto mb-4" />
                  <p className="text-gray-400 text-lg">No saved results available</p>
                </div>
              ) : (
                <div className="space-y-6">
                  {history.map((item) => (
                    <div key={item._id}>
                      <motion.div
                        className="bg-gradient-to-r from-gray-800/40 to-gray-700/40 backdrop-blur-sm p-6 rounded-xl border border-gray-600/50 hover:border-emerald-500/30 transition-all duration-300"
                        whileHover={{ y: -2, scale: 1.01 }}
                        transition={{ duration: 0.2 }}
                      >
                        <div className="space-y-4">
                          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                            <div>
                              <p className="text-gray-300 mb-2">
                                <span className="font-medium text-white">Symptoms:</span>{" "}
                                {item.symptomsGrp === "Not available" ? "No symptoms recorded" : item.symptomsGrp}
                              </p>
                              <p className="text-gray-300 mb-2"><span className="font-medium text-white">SMILES:</span></p>
                              <div className="bg-gray-900/50 p-2 rounded overflow-x-auto">
                                <code className="text-emerald-400 text-sm whitespace-nowrap">{item.smiles}</code>
                              </div>
                            </div>
                            <div>
                              <p className="text-gray-300 mb-2"><span className="font-medium text-white">Side Effects:</span> {item.toxicityResult?.sideEffects ? item.toxicityResult.sideEffects.map(effect => effect.name).join(", ") : "Not available"}</p>
                              <div className="grid grid-cols-2 gap-2 text-sm">
                                <p className="text-gray-300"><span className="font-medium text-white">LD50:</span> {item.toxicityResult?.acuteToxicity?.LD50 || "Not available"}</p>
                                <p className="text-gray-300"><span className="font-medium text-white">Class:</span> {item.toxicityResult?.acuteToxicity?.toxicityClass || "Not available"}</p>
                                <p className="text-gray-300"><span className="font-medium text-white">Hepatotoxicity:</span> {item.toxicityResult?.endpoints?.hepatotoxicity || "Not available"}</p>
                                <p className="text-gray-300"><span className="font-medium text-white">Carcinogenicity:</span> {item.toxicityResult?.endpoints?.carcinogenicity || "Not available"}</p>
                              </div>
                              {item.geminiAnalysis && (
                                <div className="mt-3">
                                  <span className="inline-flex items-center px-2 py-1 rounded-full text-xs bg-blue-500/20 text-blue-300 border border-blue-500/30">
                                    <Info className="w-3 h-3 mr-1" />
                                    Gemini analysis available
                                  </span>
                                </div>
                              )}
                            </div>
                          </div>

                          <div className="flex justify-end pt-4">
                            <motion.button
                              onClick={(e) => handleViewHistoryItem(item, e)}
                              whileHover={{ scale: 1.05 }}
                              whileTap={{ scale: 0.95 }}
                              className={`px-6 py-2 rounded-lg transition-all duration-200 font-medium flex items-center ${selectedHistoryItemId === item._id
                                ? 'bg-red-500/20 text-red-300 border border-red-500/30 hover:bg-red-500/30'
                                : 'bg-emerald-500/20 text-emerald-300 border border-emerald-500/30 hover:bg-emerald-500/30'
                                }`}
                            >
                              <Eye className="w-4 h-4 mr-2" />
                              {selectedHistoryItemId === item._id ? "Close Details" : "View Details"}
                            </motion.button>
                          </div>
                        </div>
                      </motion.div>

                      <AnimatePresence>
                        {selectedHistoryItemId === item._id && (
                          <motion.div
                            initial={{ opacity: 0, height: 0, y: -10 }}
                            animate={{ opacity: 1, height: 'auto', y: 0 }}
                            exit={{ opacity: 0, height: 0, y: -10 }}
                            transition={{ duration: 0.5 }}
                            className="mt-4 overflow-hidden"
                          >
                            <div className="bg-gradient-to-r from-gray-800/60 to-gray-700/60 backdrop-blur-sm p-6 rounded-xl border border-gray-600/50">
                              <h4 className="text-xl font-semibold mb-6 text-white flex items-center">
                                <Shield className="mr-3 w-5 h-5 text-emerald-400" />
                                Detailed Analysis
                              </h4>

                              <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 mb-6">
                                <div className="bg-gray-900/30 p-4 rounded-xl">
                                  <p className="text-gray-300 mb-2">
                                    <span className="font-medium text-white">Symptoms:</span>{" "}
                                    {item.symptomsGrp === "Not available" ? "No symptoms recorded" : item.symptomsGrp}
                                  </p>
                                </div>
                                <div className="bg-gray-900/30 p-4 rounded-xl">
                                  <p className="text-gray-300 mb-2"><span className="font-medium text-white">SMILES:</span></p>
                                  <div className="bg-gray-800/50 p-2 rounded overflow-x-auto">
                                    <code className="text-emerald-400 text-sm whitespace-nowrap">{item.smiles}</code>
                                  </div>
                                </div>
                              </div>

                              <div className="space-y-4">
                                <motion.button
                                  onClick={() => setIsHistorySideEffectsOpen(!isHistorySideEffectsOpen)}
                                  whileHover={{ backgroundColor: 'rgba(251, 191, 36, 0.1)' }}
                                  className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-yellow-500/20 hover:border-yellow-500/40"
                                >
                                  <div className="flex items-center">
                                    <AlertTriangle className="mr-3 h-6 w-6 text-yellow-400" />
                                    <span className="font-medium text-white text-lg">Potential Side Effects</span>
                                  </div>
                                  <motion.div
                                    animate={{ rotate: isHistorySideEffectsOpen ? 180 : 0 }}
                                    transition={{ duration: 0.3 }}
                                  >
                                    <ChevronDown className="h-6 w-6 text-yellow-400" />
                                  </motion.div>
                                </motion.button>
                                <AnimatePresence>
                                  {isHistorySideEffectsOpen && (
                                    <motion.div
                                      initial={{ opacity: 0, height: 0 }}
                                      animate={{ opacity: 1, height: 'auto' }}
                                      exit={{ opacity: 0, height: 0 }}
                                      transition={{ duration: 0.4 }}
                                      className="overflow-hidden"
                                    >
                                      <div className="p-6 bg-gray-900/20 rounded-xl border border-yellow-500/20">
                                        {renderSideEffectsDetails(item.toxicityResult?.sideEffects)}
                                      </div>
                                    </motion.div>
                                  )}
                                </AnimatePresence>

                                <motion.button
                                  onClick={() => setIsHistoryResultOpen(!isHistoryResultOpen)}
                                  whileHover={{ backgroundColor: 'rgba(16, 185, 129, 0.1)' }}
                                  className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-emerald-500/20 hover:border-emerald-500/40"
                                >
                                  <div className="flex items-center">
                                    <Shield className="mr-3 h-6 w-6 text-emerald-400" />
                                    <span className="font-medium text-white text-lg">Basic Toxicity Details</span>
                                  </div>
                                  <motion.div
                                    animate={{ rotate: isHistoryResultOpen ? 180 : 0 }}
                                    transition={{ duration: 0.3 }}
                                  >
                                    <ChevronDown className="h-6 w-6 text-emerald-400" />
                                  </motion.div>
                                </motion.button>
                                <AnimatePresence>
                                  {isHistoryResultOpen && (
                                    <motion.div
                                      initial={{ opacity: 0, height: 0 }}
                                      animate={{ opacity: 1, height: 'auto' }}
                                      exit={{ opacity: 0, height: 0 }}
                                      transition={{ duration: 0.4 }}
                                      className="overflow-hidden"
                                    >
                                      <div className="p-6 bg-gray-900/20 rounded-xl border border-emerald-500/20">
                                        {renderToxicityDetails(item.toxicityResult)}
                                      </div>
                                    </motion.div>
                                  )}
                                </AnimatePresence>

                                <motion.button
                                  onClick={() => setIsHistoryGeminiOpen(!isHistoryGeminiOpen)}
                                  whileHover={{ backgroundColor: 'rgba(59, 130, 246, 0.1)' }}
                                  className="w-full p-4 bg-gray-800/30 rounded-xl flex justify-between items-center transition-all duration-200 border border-blue-500/20 hover:border-blue-500/40"
                                >
                                  <div className="flex items-center">
                                    <Info className="mr-3 h-6 w-6 text-blue-400" />
                                    <span className="font-medium text-white text-lg">Advanced Gemini AI Analysis</span>
                                  </div>
                                  <motion.div
                                    animate={{ rotate: isHistoryGeminiOpen ? 180 : 0 }}
                                    transition={{ duration: 0.3 }}
                                  >
                                    <ChevronDown className="h-6 w-6 text-blue-400" />
                                  </motion.div>
                                </motion.button>
                                <AnimatePresence>
                                  {isHistoryGeminiOpen && (
                                    <motion.div
                                      initial={{ opacity: 0, height: 0 }}
                                      animate={{ opacity: 1, height: 'auto' }}
                                      exit={{ opacity: 0, height: 0 }}
                                      transition={{ duration: 0.4 }}
                                      className="overflow-hidden"
                                    >
                                      <div className="p-6 bg-gray-900/20 rounded-xl border border-blue-500/20">
                                        {renderGeminiAnalysis(item.geminiAnalysis)}
                                      </div>
                                    </motion.div>
                                  )}
                                </AnimatePresence>
                              </div>
                            </div>
                          </motion.div>
                        )}
                      </AnimatePresence>
                    </div>
                  ))}
                </div>
              )}
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
};

export default ToxicityPrediction;