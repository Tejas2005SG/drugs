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

  const API_BASE_URL = 'http://localhost:5000/api/alphafold';

  // Handle prediction submission
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
      const response = await axios.post(`${API_BASE_URL}/predict`, { uniprot_id: uniprotId }, {
        timeout: 30000,
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPredictionResult(response.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Prediction submission failed');
      setLoading(false);
    }
  };

  // Handle UniProt data fetch
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
        axios.get(`${API_BASE_URL}/uniprot/summary/${uniprotId}`, { timeout: 30000 }),
        axios.get(`${API_BASE_URL}/uniprot/annotations/${uniprotId}`, { timeout: 30000 }),
      ]);
      setSummary(summaryRes.data);
      setAnnotations(annotationsRes.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Failed to fetch UniProt data');
    } finally {
      setLoading(false);
    }
  };

  // Fetch previous jobs
  const fetchPreviousJobs = async () => {
    try {
      const response = await axios.get(`${API_BASE_URL}/previous-jobs`, {
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPreviousJobs(response.data);
    } catch (err) {
      setError('Failed to fetch previous jobs');
    }
  };

  // Poll job status
  useEffect(() => {
    let intervalId;
    if (predictionResult && predictionResult.jobId && !jobStatus?.pdbUrl && !structureRendered) {
      setLoading(true);
      const pollStatus = async () => {
        try {
          const response = await axios.get(`${API_BASE_URL}/status/${predictionResult.jobId}`, {
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

  // Initialize NGL stage and load structure
  useEffect(() => {
    if (jobStatus?.pdbUrl && !structureRendered) {
      // Initialize stage if it doesn't exist
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
        // Clear existing components if re-rendering
        stageRef.current.removeAllComponents();
      }

      const handleResize = () => {
        stageRef.current.handleResize();
      };

      window.addEventListener('resize', handleResize);

      setLoading(true);
      setError(null);

      stageRef.current.loadFile(jobStatus.pdbUrl, { ext: 'pdb' }).then((structure) => {
        // Main cartoon representation
        structure.addRepresentation('cartoon', {
          colorScheme: 'residueindex',
          colorScale: 'rainbow',
          opacity: 0.85,
          side: 'front'
        });
       

        // Secondary structure coloring (hidden by default)
        structure.addRepresentation('cartoon', {
          colorScheme: 'sstruc',
          opacity: 0.7,
          visible: false,
          name: 'secondary_structure'
        });

        // Ligands and hetero atoms
        structure.addRepresentation('ball+stick', {
          sele: 'hetero and not water',
          colorScheme: 'element',
          radius: 0.3,
          multipleBond: 'symmetric'
        });

        // Surface representation (hidden by default)
        structure.addRepresentation('surface', {
          sele: 'protein',
          colorScheme: 'electrostatic',
          opacity: 0.5,
          visible: false,
          name: 'surface_view'
        });

        // Active site residues
        structure.addRepresentation('licorice', {
          sele: 'CYS or HIS or ASP or GLU or LYS or ARG',
          colorValue: 'red',
          radius: 0.25,
          name: 'active_site'
        });

        structure.autoView(500);

        // Store component references for UI control
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

      // Cleanup function
      return () => {
        window.removeEventListener('resize', handleResize);
      };
    }
  }, [jobStatus?.pdbUrl, structureRendered]);

  // Initial fetch of previous jobs
  useEffect(() => {
    fetchPreviousJobs();
  }, []);

  // Cleanup NGL stage on unmount
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

  // Toggle visualization features
  const toggleRepresentation = (repName) => {
    if (nglComponents && nglComponents[repName]) {
      nglComponents[repName].toggleVisibility();
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-b from-blue-50 to-white">
      <div className="container mx-auto px-4 py-8">
        <h1 className="text-4xl font-bold text-center mb-8 text-blue-800">
          AlphaFold Explorer
        </h1>

        <div className="flex justify-center mb-8">
          <div className="bg-white rounded-lg shadow-lg p-1">
            <button onClick={() => setActiveTab('predict')} className={`px-6 py-3 rounded-lg font-semibold ${activeTab === 'predict' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'}`}>
              Predict Structure
            </button>
            <button onClick={() => setActiveTab('uniprot')} className={`px-6 py-3 rounded-lg font-semibold ml-2 ${activeTab === 'uniprot' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'}`}>
              UniProt Info
            </button>
            <button onClick={() => setActiveTab('previous')} className={`px-6 py-3 rounded-lg font-semibold ml-2 ${activeTab === 'previous' ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-700'}`}>
              Previous Jobs
            </button>
          </div>
        </div>

        <div className="max-w-3xl mx-auto">
          <div className="bg-white rounded-xl shadow-lg p-8">
            {activeTab !== 'previous' && (
              <form onSubmit={activeTab === 'predict' ? handlePredictionSubmit : handleUniProtSubmit}>
                <div className="mb-6">
                  <label className="block text-gray-700 font-semibold mb-2" htmlFor="uniprotId">
                    UniProt ID
                  </label>
                  <input
                    type="text"
                    id="uniprotId"
                    value={uniprotId}
                    onChange={(e) => setUniprotId(e.target.value.toUpperCase())}
                    className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
                    placeholder="Enter UniProt ID (e.g., P00520)"
                    required
                  />
                  <p className="mt-1 text-sm text-gray-500">
                    Example IDs: P00520 (BRCA2), P69905 (HBA1), P0DTD1 (SARS-CoV-2 Spike)
                  </p>
                </div>
                <button
                  type="submit"
                  disabled={loading}
                  className="w-full bg-blue-600 text-white px-6 py-3 rounded-lg font-semibold hover:bg-blue-700 disabled:bg-blue-400"
                >
                  {loading ? (
                    <>
                      <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-white inline" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                      </svg>
                      Processing...
                    </>
                  ) : activeTab === 'predict' ? 'Submit Prediction' : 'Get UniProt Info'}
                </button>
              </form>
            )}

            {activeTab === 'previous' && (
              <div>
                <h2 className="text-xl font-semibold text-blue-800 mb-4">Previous Predictions</h2>
                {previousJobs.length === 0 ? (
                  <p className="text-gray-600">No previous jobs found</p>
                ) : (
                  <div className="space-y-4">
                    {previousJobs.map(job => (
                      <div key={job.jobId} className="bg-blue-50 p-4 rounded-lg hover:bg-blue-100 transition-colors cursor-pointer" onClick={() => loadPreviousStructure(job.pdbUrl)}>
                        <p><span className="font-medium">UniProt ID:</span> {job.uniprotId}</p>
                        <p><span className="font-medium">Status:</span> 
                          <span className={`ml-2 font-semibold ${
                            job.status === 'completed' ? 'text-green-600' :
                            job.status === 'failed' ? 'text-red-600' : 'text-blue-600'
                          }`}>
                            {job.status}
                          </span>
                        </p>
                        <p><span className="font-medium">Created:</span> {new Date(job.createdAt).toLocaleString()}</p>
                        {job.completedAt && (
                          <p><span className="font-medium">Completed:</span> {new Date(job.completedAt).toLocaleString()}</p>
                        )}
                        {job.error && (
                          <p className="text-red-600 mt-2">Error: {job.error}</p>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}

            {error && (
              <div className="mt-6 p-4 bg-red-50 border border-red-200 text-red-700 rounded-lg">
                {error}
              </div>
            )}

            {activeTab === 'predict' && (predictionResult || jobStatus) && (
              <div className="mt-8 bg-blue-50 rounded-lg p-6 border border-blue-100">
                <h2 className="text-xl font-semibold text-blue-800 mb-4">Prediction Status</h2>
                {predictionResult && (
                  <div className="space-y-2 text-gray-700">
                    <p><span className="font-medium">UniProt ID:</span> {predictionResult.uniprotId}</p>
                    <p><span className="font-medium">Status:</span> 
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
                  <div className="mt-6">
                    <div className="flex space-x-2 mb-3">
                      <button 
                        onClick={() => toggleRepresentation('secondary')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200"
                      >
                        Toggle Secondary Structure
                      </button>
                      <button 
                        onClick={() => toggleRepresentation('surface')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200"
                      >
                        Toggle Surface View
                      </button>
                      <button 
                        onClick={() => toggleRepresentation('activeSite')}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200"
                      >
                        Toggle Active Sites
                      </button>
                      <button 
                        onClick={() => stageRef.current.autoView()}
                        className="px-3 py-1 bg-blue-100 text-blue-800 rounded hover:bg-blue-200"
                      >
                        Reset View
                      </button>
                    </div>
                    <div id="viewport" className="w-full h-[500px] rounded-lg border border-gray-200"></div>
                  </div>
                )}
                {jobStatus?.error && (
                  <p className="text-red-600 mt-2">Error: {jobStatus.error}</p>
                )}
              </div>
            )}

            {activeTab === 'uniprot' && summary && (
              <div className="mt-8 bg-blue-50 rounded-lg p-6 border border-blue-100">
                <h2 className="text-xl font-semibold text-blue-800 mb-4">Protein Information</h2>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-gray-700">
                  <div>
                    <h3 className="font-medium text-blue-700 mb-2">Basic Information</h3>
                    <p><span className="font-medium">Name:</span> {summary.proteinName}</p>
                    <p><span className="font-medium">Gene:</span> {summary.geneName}</p>
                    <p><span className="font-medium">Organism:</span> {summary.organism}</p>
                    <p><span className="font-medium">Length:</span> {summary.sequence.length} amino acids</p>
                  </div>
                  <div>
                    <h3 className="font-medium text-blue-700 mb-2">Functional Data</h3>
                    <p><span className="font-medium">Function:</span> {summary.function || 'Not available'}</p>
                    {annotations?.features?.length > 0 && (
                      <div className="mt-2">
                        <p className="font-medium">Key Features:</p>
                        <ul className="list-disc pl-5">
                          {annotations.features.slice(0, 3).map((feat, i) => (
                            <li key={i}>{feat.type}: {feat.description}</li>
                          ))}
                        </ul>
                      </div>
                    )}
                  </div>
                </div>
                {annotations?.features?.length > 0 && (
                  <div className="mt-4">
                    <h3 className="font-medium text-blue-700 mb-2">Sequence Features</h3>
                    <div className="bg-white p-3 rounded border border-gray-200 overflow-x-auto">
                      <pre className="text-xs font-mono">
                        {summary.sequence.match(/.{1,10}/g).map((chunk, i) => (
                          <div key={i}>
                            {(i * 10 + 1).toString().padStart(6, ' ')} {chunk.split('').join(' ')}
                          </div>
                        ))}
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