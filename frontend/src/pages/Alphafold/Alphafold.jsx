import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import * as NGL from 'ngl';

const AlphaFoldExplorer = () => {
  const [uniprotId, setUniprotId] = useState('');
  const [predictionResult, setPredictionResult] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [summary, setSummary] = useState(null);
  const [annotations, setAnnotations] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [activeTab, setActiveTab] = useState('predict');
  const [structureRendered, setStructureRendered] = useState(false);
  const stageRef = useRef(null);

  const API_BASE_URL = 'http://localhost:5000/api/alphafold';
  const STATIC_BASE_URL = 'http://localhost:5000';

  const handlePredictionSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{6,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (6-10 alphanumeric characters)');
      return;
    }
    setLoading(true);
    setError('');
    setPredictionResult(null);
    setJobStatus(null);
    setSummary(null);
    setAnnotations(null);
    setStructureRendered(false);
    try {
      const response = await axios.post(`${API_BASE_URL}/predict`, { uniprot_id: uniprotId }, { timeout: 30000 });
      console.log('Prediction response:', response.data);
      setPredictionResult(response.data);
    } catch (err) {
      console.error('Prediction error:', err);
      setError(err.response?.data?.error || err.message || 'Prediction submission failed');
      setLoading(false);
    }
  };

  const handleUniProtSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{6,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (6-10 alphanumeric characters)');
      return;
    }
    setLoading(true);
    setError('');
    setPredictionResult(null);
    setJobStatus(null);
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
      console.error('UniProt fetch error:', err);
      setError(err.response?.data?.error || err.message || 'Failed to fetch UniProt data');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    let intervalId;
    if (predictionResult && predictionResult.jobId && !jobStatus?.pdb_url && !structureRendered) {
      setLoading(true);
      const pollStatus = async () => {
        try {
          const response = await axios.get(`${API_BASE_URL}/status/${predictionResult.jobId}`);
          console.log('Job status:', response.data);
          setJobStatus(response.data);
          if (response.data.status === 'completed' || response.data.status === 'failed') {
            clearInterval(intervalId);
          }
        } catch (err) {
          console.error('Status poll error:', err);
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

  useEffect(() => {
    if (jobStatus?.pdb_url && !stageRef.current) {
      stageRef.current = new NGL.Stage('viewport', {
        backgroundColor: 'white',
        quality: 'high',
        cameraType: 'perspective',
        cameraFov: 40,
      });
  
      if (stageRef.current.viewer.renderer.useLegacyLights !== undefined) {
        stageRef.current.viewer.renderer.useLegacyLights = false;
      }
  
      const handleResize = () => stageRef.current.handleResize();
      window.addEventListener('resize', handleResize);
  
      const pdbUrl = `${STATIC_BASE_URL}${jobStatus.pdb_url}`;
      console.log('Loading PDB from:', pdbUrl);
  
      stageRef.current.loadFile(pdbUrl, { ext: 'pdb' }).then((structure) => {
        structure.addRepresentation('cartoon', {
          color: 'sstruc', // Color by secondary structure
          opacity: 0.9,
          smoothSheet: true,
          quality: 'high'
        });
  
        structure.addRepresentation('ball+stick', {
          sele: 'hetero and not water',
          color: 'element',
          opacity: 0.8
        });
  
        structure.addRepresentation('spacefill', {
          sele: 'water',
          color: 'element',
          opacity: 0.5,
          scale: 0.25
        });
  
        stageRef.current.setParameters({
          clipNear: 0,
          clipFar: 100,
          clipDist: 10,
          fogNear: 50,
          fogFar: 100,
        });
  
        stageRef.current.setParameters({
          lightColor: 0xffffff, // White light
          lightIntensity: 1.0,
          ambientColor: 0xffffff, // White ambient light
          ambientIntensity: 0.4,
        });
  
        structure.autoView(1000);
        setStructureRendered(true);
        setLoading(false);
      }).catch((err) => {
        console.error('Error loading PDB:', err);
        setError('Failed to load protein structure. The PDB file may be invalid or inaccessible.');
        setLoading(false);
        setStructureRendered(true);
      });
    }
  
    return () => {
      if (stageRef.current) {
        stageRef.current.dispose();
        stageRef.current = null;
      }
      window.removeEventListener('resize', () => {});
    };
  }, [jobStatus]);

  return (
    <div className="min-h-screen bg-gradient-to-b from-blue-50 to-white">
      <div className="container mx-auto px-4 py-8">
        <h1 className="text-4xl font-bold text-center mb-8 text-blue-800">
          AlphaFold Explorer
        </h1>
        
        <div className="flex justify-center mb-8">
          <div className="bg-white rounded-lg shadow-lg p-1">
            <button
              onClick={() => setActiveTab('predict')}
              className={`px-6 py-3 rounded-lg font-semibold transition-all duration-200 ${
                activeTab === 'predict' 
                  ? 'bg-blue-600 text-white shadow-md' 
                  : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
              }`}
            >
              Predict Structure
            </button>
            <button
              onClick={() => setActiveTab('uniprot')}
              className={`px-6 py-3 rounded-lg font-semibold ml-2 transition-all duration-200 ${
                activeTab === 'uniprot' 
                  ? 'bg-blue-600 text-white shadow-md' 
                  : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
              }`}
            >
              UniProt Info
            </button>
          </div>
        </div>

        <div className="max-w-3xl mx-auto">
          <div className="bg-white rounded-xl shadow-lg p-8">
            <form onSubmit={activeTab === 'predict' ? handlePredictionSubmit : handleUniProtSubmit}>
              <div className="mb-6">
                <label className="block text-gray-700 font-semibold mb-2" htmlFor="uniprotId">
                  UniProt ID
                </label>
                <input
                  type="text"
                  id="uniprotId"
                  value={uniprotId}
                  onChange={(e) => setUniprotId(e.target.value)}
                  className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                  placeholder="Enter UniProt ID (e.g., P00520)"
                  required
                />
              </div>
              <button
                type="submit"
                disabled={loading}
                className="w-full bg-blue-600 text-white px-6 py-3 rounded-lg font-semibold hover:bg-blue-700 disabled:bg-blue-400 transition-all duration-200 transform hover:scale-[1.02]"
              >
                {loading
                  ? 'Processing...'
                  : activeTab === 'predict'
                  ? 'Submit Prediction'
                  : 'Get UniProt Info'}
              </button>
            </form>

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
                    <p><span className="font-medium">Job ID:</span> {predictionResult.jobId}</p>
                    <p><span className="font-medium">Initial Status:</span> {predictionResult.status}</p>
                    <p><span className="font-medium">UniProt ID:</span> {predictionResult.uniprotId}</p>
                  </div>
                )}
                {jobStatus && (
                  <div className="mt-4">
                    <p className="font-medium text-gray-700">
                      Current Status: <span className="text-blue-600">{jobStatus.status}</span>
                    </p>
                    {jobStatus.pdb_url && (
                      <div className="mt-6">
                        <h3 className="text-lg font-semibold text-blue-800 mb-3">Protein Structure</h3>
                        <div 
                          id="viewport" 
                          className="w-full h-[500px] rounded-lg overflow-hidden border border-gray-200 shadow-inner"
                        ></div>
                      </div>
                    )}
                    {jobStatus.error && (
                      <p className="text-red-600 mt-2">Error: {jobStatus.error}</p>
                    )}
                  </div>
                )}
              </div>
            )}

            {activeTab === 'uniprot' && summary && (
              <div className="mt-8 bg-blue-50 rounded-lg p-6 border border-blue-100">
                <h2 className="text-xl font-semibold text-blue-800 mb-4">UniProt Protein Summary</h2>
                <div className="space-y-3 text-gray-700">
                  <p><span className="font-medium">Accession:</span> {summary.accession}</p>
                  <p><span className="font-medium">ID:</span> {summary.id}</p>
                  <p><span className="font-medium">Protein Name:</span> {summary.protein.recommendedName.fullName.value}</p>
                  <p><span className="font-medium">Alternative Names:</span> {summary.protein.alternativeName.map(n => n.fullName.value).join(', ')}</p>
                  <p><span className="font-medium">Organism:</span> {summary.organism.names.find(n => n.type === 'scientific')?.value} ({summary.organism.names.find(n => n.type === 'common')?.value})</p>
                  <p><span className="font-medium">Function:</span> {summary.comments.find(c => c.type === 'FUNCTION')?.text[0].value}</p>
                  <p><span className="font-medium">Gene:</span> {summary.gene[0].name.value} {summary.gene[0].synonyms ? `(Synonyms: ${summary.gene[0].synonyms.map(s => s.value).join(', ')})` : ''}</p>
                  <p><span className="font-medium">Protein Existence:</span> {summary.proteinExistence}</p>
                  <p><span className="font-medium">Entry Info:</span> Type: {summary.info.type}, Created: {summary.info.created}, Modified: {summary.info.modified}, Version: {summary.info.version}</p>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
      
      <footer className="mt-12 pb-8 text-center text-gray-600">
        <p className="text-sm">Powered by AlphaFold & UniProt API</p>
      </footer>
    </div>
  );
};

export default AlphaFoldExplorer;