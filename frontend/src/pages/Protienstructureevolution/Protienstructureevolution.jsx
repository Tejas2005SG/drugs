import React, { useState, useEffect } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import { jsPDF } from "jspdf";
import { FiCopy, FiDownload, FiInfo, FiChevronDown, FiChevronUp, FiRotateCw } from "react-icons/fi";
import { Dna, LogIn } from "lucide-react";

const ProteinStructureEvolution = () => {
  const [formData, setFormData] = useState({ 
    smilesoffirst: "", 
    smilesofsecond: "", 
    newmoleculetitle: "" 
  });
  const [molecules, setMolecules] = useState([]);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [realTimeOutput, setRealTimeOutput] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [expandedInfoId, setExpandedInfoId] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError("Please log in to access this feature");
        return;
      }
      await fetchAllMolecules();
    };
    initializeApp();
  }, [checkAuth]);

  const fetchAllMolecules = async () => {
    if (!user?._id) return;
    setLoading(true);
    try {
      const { data } = await axios.get("/api/protein/generatednewmolecule");
      setMolecules(data.molecules || []);
      if (data.molecules?.length > 0) {
        setSelectedMolecule(data.molecules[0]);
      }
    } catch (err) {
      setError("Failed to fetch molecules. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
  };

  const handleGenerate = async (e) => {
    e.preventDefault();
    if (!formData.smilesoffirst || !formData.smilesofsecond || !formData.newmoleculetitle) {
      setError("All fields are required");
      return;
    }
    if (!user?._id) {
      setError("Please log in to generate molecules");
      return;
    }

    setLoading(true);
    setError(null);
    setRealTimeOutput("");

    try {
      const response = await fetch(`${API_BASE_URL}/api/protein/generatenewmolecule/${user._id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(formData),
      });

      if (!response.ok) throw new Error("Generation failed");

      const data = await response.json();
      setMolecules(prev => [data.molecule, ...prev]);
      setSelectedMolecule(data.molecule);
      toast.success("Molecule generated successfully!");
      setFormData({ smilesoffirst: "", smilesofsecond: "", newmoleculetitle: "" });
    } catch (err) {
      setError(err.message || "Generation failed");
    } finally {
      setLoading(false);
    }
  };

  const handleCopySmiles = (smiles) => {
    navigator.clipboard.writeText(smiles)
      .then(() => toast.success("SMILES copied to clipboard!"))
      .catch(() => toast.error("Failed to copy SMILES"));
  };

  const toggleInfo = (id) => {
    setExpandedInfoId(prev => prev === id ? null : id);
  };

  const exportToPDF = (molecule) => {
    const doc = new jsPDF();
    doc.text(`Molecule Report - ${molecule.newmoleculetitle}`, 10, 10);
    doc.text(`SMILES: ${molecule.newSmiles || "N/A"}`, 10, 20);
    doc.text(`Generated on: ${new Date(molecule.created).toLocaleString()}`, 10, 30);
    doc.save(`molecule-${molecule._id}.pdf`);
  };

  if (checkingAuth) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gray-50">
        <FiRotateCw className="animate-spin h-12 w-12 text-blue-600 mb-4" />
        <p className="text-lg text-gray-700">Verifying authentication...</p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gray-50">
        <div className="bg-white p-8 rounded-xl shadow-lg max-w-md w-full text-center">
          <Dna className="h-12 w-12 text-blue-600 mx-auto mb-4" />
          <h2 className="text-2xl font-bold text-gray-800 mb-4">Access Required</h2>
          <p className="text-gray-600 mb-6">
            Please log in to access the Molecule Generator
          </p>
          <button
            onClick={() => window.location.href = '/login'}
            className="w-full px-4 py-3 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition flex items-center justify-center"
          >
            <LogIn className="mr-2" size={18} />
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 py-8 px-4">
      <div className="max-w-5xl mx-auto">
        <div className="text-center mb-10">
          <h1 className="text-3xl font-bold text-gray-800 mb-2 flex items-center justify-center">
            <Dna className="mr-2 text-blue-600" size={28} />
            New Drug Molecule Generation
          </h1>
          <p className="text-xs text-blue-700 font-semibold">(Powered by Gemini)</p>
        </div>

        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-6 flex justify-between items-center">
            <span>{error}</span>
            <button onClick={() => setError(null)} className="text-red-700 hover:text-red-900">
              ×
            </button>
          </div>
        )}

        <div className="space-y-6">
          {/* Generate New Molecule */}
          <div className="bg-white p-6 rounded-xl shadow-md">
            <h2 className="text-xl font-semibold text-gray-800 mb-4">Generate New Molecule</h2>
            <form onSubmit={handleGenerate} className="space-y-4">
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                {["newmoleculetitle", "smilesoffirst", "smilesofsecond"].map((field) => (
                  <div key={field}>
                    <label className="block text-sm font-medium text-gray-700 mb-1">
                      {field === "newmoleculetitle" ? "Molecule Title" : 
                       field.replace(/([A-Z])/g, " $1").trim()}
                    </label>
                    <input
                      type="text"
                      name={field}
                      value={formData[field]}
                      onChange={handleInputChange}
                      placeholder={
                        field === "newmoleculetitle" ? "Enter molecule title" : 
                        `Enter ${field.replace(/([A-Z])/g, " $1").trim().toLowerCase()}`
                      }
                      disabled={loading}
                      required
                      className="w-full p-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-400"
                    />
                  </div>
                ))}
              </div>
              <button
                type="submit"
                disabled={loading}
                className={`w-full py-2 px-4 rounded-lg text-white font-medium ${
                  loading ? 'bg-blue-400' : 'bg-blue-600 hover:bg-blue-700'
                } transition flex items-center justify-center`}
              >
                {loading ? (
                  <>
                    <FiRotateCw className="animate-spin mr-2" />
                    Generating...
                  </>
                ) : 'Generate Molecule'}
              </button>
            </form>
          </div>

          {/* Real-Time Output */}
          {realTimeOutput && (
            <div className="bg-white p-6 rounded-xl shadow-md">
              <h3 className="text-lg font-semibold text-gray-800 mb-3">Live Output</h3>
              <div className="p-3 bg-gray-50 border border-gray-200 rounded-lg max-h-48 overflow-y-auto">
                <pre className="text-sm text-gray-700 whitespace-pre-wrap">{realTimeOutput}</pre>
              </div>
            </div>
          )}

          {/* Molecule Details */}
          {selectedMolecule && (
            <div className="bg-white p-6 rounded-xl shadow-md">
              <div className="flex justify-between items-center mb-4">
                <h2 className="text-xl font-semibold text-gray-800">Molecule Details</h2>
                <div className="flex space-x-2">
                  <button
                    onClick={() => exportToPDF(selectedMolecule)}
                    className="flex items-center px-3 py-1 bg-green-600 text-white rounded-lg hover:bg-green-700"
                  >
                    <FiDownload className="mr-1" />
                    PDF
                  </button>
                  <button
                    onClick={() => toggleInfo(selectedMolecule.id)}
                    className="flex items-center px-3 py-1 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
                  >
                    <FiInfo className="mr-1" />
                    {expandedInfoId === selectedMolecule.id ? 'Hide' : 'Info'}
                  </button>
                </div>
              </div>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {[
                  { label: "Title", value: selectedMolecule.newmoleculetitle },
                  { label: "SMILES", value: selectedMolecule.newSmiles || "N/A", copy: true },
                  { label: "IUPAC Name", value: selectedMolecule.newIupacName || "N/A" },
                  { label: "Created", value: new Date(selectedMolecule.created).toLocaleString() },
                ].map(({ label, value, copy }) => (
                  <div key={label} className="bg-gray-50 p-3 rounded-lg">
                    <p className="text-sm font-medium text-gray-500">{label}</p>
                    <div className="flex items-center mt-1">
                      <p className="text-gray-800 break-all">{value}</p>
                      {copy && (
                        <button 
                          onClick={() => handleCopySmiles(value)}
                          className="ml-2 text-gray-500 hover:text-blue-600"
                        >
                          <FiCopy />
                        </button>
                      )}
                    </div>
                  </div>
                ))}
              </div>

              {expandedInfoId === selectedMolecule.id && (
                <div className="mt-4 p-4 bg-gray-50 border border-gray-200 rounded-lg">
                  <h4 className="font-medium text-gray-700 mb-2">Additional Information</h4>
                  <div className="prose prose-sm max-w-none">
                    {selectedMolecule.information || "No additional information available"}
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Your Molecules */}
          <div className="bg-white p-6 rounded-xl shadow-md">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-semibold text-gray-800">Your Molecules</h2>
              <button 
                onClick={fetchAllMolecules}
                disabled={loading}
                className="text-blue-600 hover:text-blue-800"
              >
                <FiRotateCw className={`${loading ? 'animate-spin' : ''}`} />
              </button>
            </div>

            {loading && molecules.length === 0 ? (
              <div className="flex justify-center py-8">
                <FiRotateCw className="animate-spin text-blue-600 text-2xl" />
              </div>
            ) : molecules.length === 0 ? (
              <p className="text-gray-500 text-center py-8">No molecules generated yet</p>
            ) : (
              <div className="space-y-3">
                {molecules.map((molecule) => (
                  <div 
                    key={molecule.id} 
                    className={`p-4 border rounded-lg cursor-pointer transition ${
                      selectedMolecule?.id === molecule.id 
                        ? 'border-blue-500 bg-blue-50' 
                        : 'border-gray-200 hover:border-blue-300'
                    }`}
                    onClick={() => setSelectedMolecule(molecule)}
                  >
                    <div className="flex justify-between items-start">
                      <div>
                        <h3 className="font-medium text-gray-800">{molecule.newmoleculetitle}</h3>
                        <p className="text-sm text-gray-500 mt-1">
                          {new Date(molecule.created).toLocaleString()}
                        </p>
                      </div>
                      <button
                        onClick={(e) => {
                          e.stopPropagation();
                          toggleInfo(molecule.id);
                        }}
                        className="text-blue-600 hover:text-blue-800"
                      >
                        {expandedInfoId === molecule.id ? <FiChevronUp /> : <FiChevronDown />}
                      </button>
                    </div>

                    {expandedInfoId === molecule.id && (
                      <div className="mt-3 pt-3 border-t border-gray-200">
                        <div className="flex items-center text-sm text-gray-600 mb-2">
                          <span className="font-medium mr-2">SMILES:</span>
                          <span className="break-all">{molecule.newSmiles || "N/A"}</span>
                          <button 
                            onClick={(e) => {
                              e.stopPropagation();
                              handleCopySmiles(molecule.newSmiles);
                            }}
                            className="ml-2 text-gray-500 hover:text-blue-600"
                          >
                            <FiCopy size={14} />
                          </button>
                        </div>
                        <div className="prose prose-sm max-w-none text-gray-600">
                          {molecule.information || "No additional information available"}
                        </div>
                      </div>
                    )}
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureEvolution;