import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import * as NGL from 'ngl';

const AlphaFoldExplorer = () => {
  const [uniprotId, setUniprotId] = useState('');
  const [predictionResult, setPredictionResult] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [summary, setSummary] = useState(null);
  const [annotations, setAnnotations] = useState(null);
  const [previousJobs, setPreviousJobs] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [activeTab, setActiveTab] = useState('predict');
  const [structureRendered, setStructureRendered] = useState(false);
  const [nglComponents, setNglComponents] = useState(null);
  const stageRef = useRef(null);

  const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
  const axiosInstance = axios.create({
    baseURL: import.meta.mode === "development" ? API_BASE_URL : '/api',
    withCredentials: true,
  });

  // Handle prediction submission (unchanged)
  const handlePredictionSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{5,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (6-10 alphanumeric characters)');
      return;
    }
    setLoading(true);
    setError('');
    setPredictionResult(null);
    setJobStatus(null);
    setStructureRendered(false);
    try {
      const response = await axiosInstance.post(`${API_BASE_URL}/alphafold/predict`, { uniprot_id: uniprotId }, {
        timeout: 30000,
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPredictionResult(response.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Prediction submission failed');
      setLoading(false);
    }
  };

  // Handle UniProt data fetch (unchanged)
  const handleUniProtSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{6,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (6-10 alphanumeric characters)');
      return;
    }
    setLoading(true);
    setError('');
    setSummary(null);
    setAnnotations(null);

    try {
      const [summaryRes, annotationsRes] = await Promise.all([
        axiosInstance.get(`${API_BASE_URL}/alphafold/uniprot/summary/${uniprotId}`, { timeout: 30000 }),
        axiosInstance.get(`${API_BASE_URL}/alphafold/uniprot/annotations/${uniprotId}`, { timeout: 30000 }),
      ]);
      setSummary(summaryRes.data);
      setAnnotations(annotationsRes.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Failed to fetch UniProt data');
    } finally {
      setLoading(false);
    }
  };

  // Fetch previous jobs (unchanged)
  const fetchPreviousJobs = async () => {
    try {
      const response = await axiosInstance.get(`${API_BASE_URL}/alphafold/previous-jobs`, {
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPreviousJobs(response.data);
    } catch (err) {
      setError('Failed to fetch previous jobs');
    }
  };

  // Poll job status (unchanged)
  useEffect(() => {
    let intervalId;
    if (predictionResult && predictionResult.jobId && !jobStatus?.pdbUrl && !structureRendered) {
      setLoading(true);
      const pollStatus = async () => {
        try {
          const response = await axiosInstance.get(`${API_BASE_URL}/alphafold/status/${predictionResult.jobId}`, {
            headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
          });
          setJobStatus(response.data);
          if (response.data.status === 'completed' || response.data.status === 'failed') {
            clearInterval(intervalId);
            fetchPreviousJobs();
          }
        } catch (err) {
          setError(err.response?.data?.error || 'Failed to fetch job status');
          clearInterval(intervalId);
          setLoading(false);
        }
      };

      pollStatus();
      intervalId = setInterval(pollStatus, 5000);
    }

    return () => {
      if (intervalId) clearInterval(intervalId);
    };
  }, [predictionResult, structureRendered]);

  // Initialize NGL stage and load structure (unchanged)
  useEffect(() => {
    if (jobStatus?.pdbUrl && !structureRendered) {
      if (!stageRef.current) {
        stageRef.current = new NGL.Stage('viewport', {
          backgroundColor: 'black',
          quality: 'high',
          antialias: true,
          clipNear: 0,
          clipFar: 100,
          clipDist: 10
        });
      } else {
        stageRef.current.removeAllComponents();
      }

      const handleResize = () => {
        stageRef.current.handleResize();
      };

      window.addEventListener('resize', handleResize);

      setLoading(true);
      setError(null);

      stageRef.current.loadFile(jobStatus.pdbUrl, { ext: 'pdb' }).then((structure) => {
        structure.addRepresentation('cartoon', {
          colorScheme: 'residueindex',
          colorScale: 'rainbow',
          opacity: 0.85,
          side: 'front'
        });

        structure.addRepresentation('cartoon', {
          colorScheme: 'sstruc',
          opacity: 0.7,
          visible: false,
          name: 'secondary_structure'
        });

        structure.addRepresentation('ball+stick', {
          sele: 'hetero and not water',
          colorScheme: 'element',
          radius: 0.3,
          multipleBond: 'symmetric'
        });

        structure.addRepresentation('surface', {
          sele: 'protein',
          colorScheme: 'electrostatic',
          opacity: 0.5,
          visible: false,
          name: 'surface_view'
        });

        structure.addRepresentation('licorice', {
          sele: 'CYS or HIS or ASP or GLU or LYS or ARG',
          colorValue: 'red',
          radius: 0.25,
          name: 'active_site'
        });

        structure.autoView(500);

        setNglComponents({
          main: structure,
          secondary: stageRef.current.getRepresentationsByName('secondary_structure')[0],
          surface: stageRef.current.getRepresentationsByName('surface_view')[0],
          activeSite: stageRef.current.getRepresentationsByName('active_site')[0]
        });

        setStructureRendered(true);
        setLoading(false);
      }).catch((err) => {
        console.error('Structure loading failed:', err);
        setError(`Failed to load structure: ${err.message}`);
        setLoading(false);
        setStructureRendered(true);
      });

      return () => {
        window.removeEventListener('resize', handleResize);
      };
    }
  }, [jobStatus?.pdbUrl, structureRendered]);

  // Reset structure and prediction data when switching tabs (unchanged)
  useEffect(() => {
    setPredictionResult(null);
    setJobStatus(null);
    setStructureRendered(false);
    setNglComponents(null);
    setError('');
    setLoading(false);

    if (stageRef.current) {
      stageRef.current.removeAllComponents();
    }
  }, [activeTab]);

  // Initial fetch of previous jobs (unchanged)
  useEffect(() => {
    fetchPreviousJobs();
  }, []);

  // Cleanup NGL stage onhomme unmount (unchanged)
  useEffect(() => {
    return () => {
      if (stageRef.current) {
        stageRef.current.dispose();
        stageRef.current = null;
      }
    };
  }, []);

  const loadPreviousStructure = (pdbUrl) => {
    setJobStatus({ pdbUrl });
    setStructureRendered(false);
  };

  const toggleRepresentation = (repName) => {
    if (nglComponents && nglComponents[repName]) {
      nglComponents[repName].toggleVisibility();
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-b from-blue-50 to-white flex flex-col">
      <div className="container mx-auto px-4 py-6 sm:py-8 flex-grow">
        <h1 className="text-2xl sm:text-3xl md:text-4xl font-bold text-center mb-6 sm:mb-8 text-blue-800">
          AlphaFold Explorer
        </h1>

        {/* Tab Navigation */}
        <div className="flex flex-col sm:flex-row justify-center mb-6 sm:mb-8 gap-2 sm:gap-0">
          <div className="bg-white rounded-lg shadow-lg p-1 flex flex-row justify-center gap-2 sm:gap-0 flex-wrap">
            <button
              onClick={() => setActiveTab('predict')}
              className={`px-4 py-2 sm:px-6 sm:py-3 rounded-lg font-semibold text-sm sm:text-base w-full sm:w-auto ${
                activeTab === 'predict' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'
              }`}
            >
              Predict Structure
            </button>
            <button
              onClick={() => setActiveTab('uniprot')}
              className={`px-4 py-2 sm:px-6 sm:py-3 rounded-lg font-semibold text-sm sm:text-base w-full sm:w-auto ${
                activeTab === 'uniprot' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'
              }`}
            >
              UniProt Info
            </button>
            <button
              onClick={() => setActiveTab('previous')}
              className={`px-4 py-2 sm:px-6 sm:py-3 rounded-lg font-semibold text-sm sm:text-base w-full sm:w-auto ${
                activeTab === 'previous' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'
              }`}
            >
              Previous Jobs
            </button>
          </div>
        </div>

        <div className="w-full max-w-4xl mx-auto">
          <div className="bg-white rounded-xl shadow-lg p-4 sm:p-6 md:p-8">
            {activeTab !== 'previous' && (
              <form onSubmit={activeTab === 'predict' ? handlePredictionSubmit : handleUniProtSubmit}>
                <div className="mb-4 sm:mb-6">
                  <label className="block text-gray-700 font-semibold mb-2 text-sm sm:text-base" htmlFor="uniprotId">
                    UniProt ID
                  </label>
                  <input
                    type="text"
                    id="uniprotId"
                    value={uniprotId}
                    onChange={(e) => setUniprotId(e.target.value.toUpperCase())}
                    className="w-full p-2 sm:p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 text-sm sm:text-base"
                    placeholder="Enter UniProt ID (e.g., P00520)"
                    required
                  />
                  <p className="mt-1 text-xs sm:text-sm text-gray-500">
                    Example IDs: P00520 (BRCA2), P69905 (HBA1), P0DTD1 (SARS-CoV-2 Spike)
                  </p>
                </div>
                <button
                  type="submit"
                  disabled={loading}
                  className="w-full bg-blue-600 text-white px-4 py-2 sm:px-6 sm:py-3 rounded-lg font-semibold text-sm sm:text-base hover:bg-blue-700 disabled:bg-blue-400"
                >
                  {loading ? (
                    <div className="flex items-center justify-center">
                      <svg className="animate-spin -ml-1 mr-2 sm:mr-3 h-4 w-4 sm:h-5 sm:w-5 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                      </svg>
                      Processing...
                    </div>
                  ) : activeTab === 'predict' ? 'Submit Prediction' : 'Get UniProt Info'}
                </button>
              </form>
            )}

            {activeTab === 'previous' && (
              <div>
                <h2 className="text-lg sm:text-xl font-semibold text-blue-800 mb-4">Previous Predictions</h2>
                {previousJobs.length === 0 ? (
                  <p className="text-gray-600 text-sm sm:text-base">No previous jobs found</p>
                ) : (
                  <div className="space-y-4">
                    {previousJobs.map(job => (
                      <div
                        key={job.jobId}
                        className="bg-blue-50 p-4 rounded-lg hover:bg-blue-100 transition-colors cursor-pointer"
                        onClick={() => loadPreviousStructure(job.pdbUrl)}
                      >
                        <p className="text-sm sm:text-base"><span className="font-medium">UniProt ID:</span> {job.uniprotId}</p>
                        <p className="text-sm sm:text-base">
                          <span className="font-medium">Status:</span>
                          <span className={`ml-2 font-semibold ${
                            job.status === 'completed' ? 'text-green-600' :
                            job.status === 'failed' ? 'text-red-600' : 'text-blue-600'
                          }`}>
                            {job.status}
                          </span>
                        </p>
                        <p className="text-sm sm:text-base"><span className="font-medium">Created:</span> {new Date(job.createdAt).toLocaleString()}</p>
                        {job.completedAt && (
                          <p className="text-sm sm:text-base"><span className="font-medium">Completed:</span> {new Date(job.completedAt).toLocaleString()}</p>
                        )}
                        {job.error && (
                          <p className="text-red-600 mt-2 text-sm sm:text-base">Error: {job.error}</p>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}

            {error && (
              <div className="mt-4 sm:mt-6 p-4 bg-red-50 border border-red-200 text-red-700 rounded-lg text-sm sm:text-base">
                {error}
              </div>
            )}

            {activeTab === 'predict' && (predictionResult || jobStatus) && (
              <div className="mt-6 sm:mt-8 bg-blue-50 rounded-lg p-4 sm:p-6 border border-blue-100">
                <h2 className="text-lg sm:text-xl font-semibold text-blue-800 mb-4">Prediction Status</h2>
                {predictionResult && (
                  <div className="space-y-2 text-gray-700 text-sm sm:text-base">
                    <p><span className="font-medium">UniProt ID:</span> {predictionResult.uniprotId}</p>
                    <p>
                      <span className="font-medium">Status:</span>
                      <span className={`ml-2 font-semibold ${
                        jobStatus?.status === 'completed' ? 'text-green-600' :
                        jobStatus?.status === 'failed' ? 'text-red-600' : 'text-blue-600'
                      }`}>
                        {jobStatus?.status || predictionResult.status}
                      </span>
                    </p>
                    <p><span className="font-medium">Submitted:</span> {new Date(predictionResult.createdAt).toLocaleString()}</p>
                  </div>
                )}
                {jobStatus?.pdbUrl && (
                  <div className="mt-4 sm:mt-6">
                    <div className="flex flex-wrap gap-2 mb-3">
                      <button
                        onClick={() => toggleRepresentation('secondary')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200 text-sm"
                      >
                        Toggle Secondary
                      </button>
                      <button
                        onClick={() => toggleRepresentation('surface')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200 text-sm"
                      >
                        Toggle Surface
                      </button>
                      <button
                        onClick={() => toggleRepresentation('activeSite')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200 text-sm"
                      >
                        Toggle Active Sites
                      </button>
                      <button
                        onClick={() => stageRef.current.autoView()}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200 text-sm"
                      >
                        Reset View
                      </button>
                    </div>
                    <div
                      id="viewport"
                      className="w-full h-[300px] sm:h-[400px] md:h-[500px] rounded-lg border border-gray-200"
                    ></div>
                  </div>
                )}
                {jobStatus?.error && (
                  <p className="text-red-600 mt-2 text-sm sm:text-base">Error: {jobStatus.error}</p>
                )}
              </div>
            )}

            {activeTab === 'uniprot' && summary && (
              <div className="mt-6 sm:mt-8 bg-blue-50 rounded-lg p-4 sm:p-6 border border-blue-100">
                <h2 className="text-lg sm:text-xl font-semibold text-blue-800 mb-4">Protein Information</h2>
                <div className="flex flex-col gap-4 text-gray-700 text-sm sm:text-base">
                  <div>
                    <h3 className="font-medium text-blue-700 mb-2">Basic Information</h3>
                    <p><span className="font-medium">Protein Name:</span> {summary.protein?.recommendedName?.fullName?.value || 'Not available'}</p>
                    <p><span className="font-medium">Gene:</span> {summary.gene?.[0]?.name?.value || 'Not available'}</p>
                    <p><span className="font-medium">Organism:</span> {summary.organism?.scientificName || 'Not available'}</p>
                    <p><span className="font-medium">Length:</span> {summary.sequence?.length ? `${summary.sequence.length} amino acids` : 'Not available'}</p>
                  </div>
                  <div className="bg-blue-50 rounded-lg p-4 border border-blue-100">
                    <h3 className="text-base sm:text-lg font-semibold text-blue-800 mb-3 border-b border-blue-200 pb-2">
                      Functional Summary
                    </h3>
                    <div className="mb-4">
                      <h4 className="font-medium text-blue-700 flex items-center mb-2 text-sm sm:text-base">
                        <svg className="w-4 h-4 mr-2 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
                        </svg>
                        Protein Function
                      </h4>
                      <div className="bg-white p-3 rounded-md shadow-sm text-sm sm:text-base">
                        {summary.comments?.find(c => c.type === 'FUNCTION')?.text?.[0]?.value || (
                          <span className="text-gray-500">Not specified</span>
                        )}
                      </div>
                    </div>
                    {annotations?.features?.length > 0 && (
                      <div className="mb-4">
                        <h4 className="font-medium text-blue-700 flex items-center mb-2 text-sm sm:text-base">
                          <svg className="w-4 h-4 mr-2 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                          </svg>
                          Key Features
                        </h4>
                        <div className="bg-white p-3 rounded-md shadow-sm">
                          <ul className="space-y-2">
                            {annotations.features.slice(0, 5).map((feat, i) => (
                              <li key={i} className="flex items-start text-sm sm:text-base">
                                <span className="inline-block bg-blue-100 text-blue-800 text-xs px-2 py-1 rounded-full mr-2 mt-1">
                                  {feat.type}
                                </span>
                                <span className="flex-1">{feat.description || 'No description available'}</span>
                              </li>
                            ))}
                          </ul>
                        </div>
                      </div>
                    )}
                    <div className="grid grid-cols-1 gap-4">
                      <div>
                        <h4 className="font-medium text-blue-700 mb-1 text-sm sm:text-base">Subcellular Location</h4>
                        <div className="bg-white p-2 rounded text-sm sm:text-base">
                          {summary.comments?.find(c => c.type === 'SUBCELLULAR_LOCATION')?.text?.[0]?.value || 'Unknown'}
                        </div>
                      </div>
                      <div>
                        <h4 className="font-medium text-blue-700 mb-1 text-sm sm:text-base">Catalytic Activity</h4>
                        <div className="bg-white p-2 rounded text-sm sm:text-base">
                          {summary.comments?.find(c => c.type === 'CATALYTIC_ACTIVITY')?.text?.[0]?.value || 'Not enzymatic'}
                        </div>
                      </div>
                    </div>
                  </div>
                </div>
                {summary?.sequence?.sequence && typeof summary.sequence.sequence === 'string' && (
                  <div className="mt-4">
                    <h3 className="font-medium text-blue-700 mb-2 text-sm sm:text-base">Sequence Features</h3>
                    <div className="bg-white p-3 rounded border border-gray-200 overflow-x-auto">
                      <pre className="text-xs font-mono whitespace-pre-wrap break-all">
                        {summary.sequence.sequence.match(/.{1,10}/g)?.map((chunk, i) => (
                          <div key={i}>
                            {(i * 10 + 1).toString().padStart(6, ' ')} {chunk.split('').join(' ')}
                          </div>
                        )) || 'No sequence data available'}
                      </pre>
                    </div>
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default AlphaFoldExplorer;