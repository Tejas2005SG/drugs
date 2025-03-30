import React, { useState, useEffect } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import { jsPDF } from "jspdf";
import { FiCopy } from "react-icons/fi";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode==="development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

const ProteinStructureEvolution = () => {
  const [formData, setFormData] = useState({ smilesoffirst: "", smilesofsecond: "", newmoleculetitle: "" });
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
        setError("Authentication failed. Please log in.");
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
      const { data } = await axiosInstance.get("/protein/generatednewmolecule");
      const fetchedMolecules = data.molecules || [];
      setMolecules(fetchedMolecules);
    } catch (err) {
      console.error("Fetch error:", err);
      setMolecules([]);
      setError(err.response?.data?.message || "Failed to fetch molecules.");
    } finally {
      setLoading(false);
    }
  };

  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setFormData((prev) => ({ ...prev, [name]: value }));
  };

  const handleGenerate = async (e) => {
    e.preventDefault();
    const { smilesoffirst, smilesofsecond, newmoleculetitle } = formData;
    if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
      setError("All fields are required.");
      return;
    }
    if (!user?._id) {
      setError("Please log in to generate a molecule.");
      return;
    }

    setLoading(true);
    setError(null);
    setRealTimeOutput("");

    try {
      console.log("Generating molecule for user ID:", user._id);
      const response = await fetch(`${API_BASE_URL}/protein/generatenewmolecule/${user._id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json", "Accept": "text/event-stream" },
        credentials: "include",
        body: JSON.stringify(formData),
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        const errorMessage = errorData.message || `Server error: ${response.status}`;
        throw new Error(errorMessage);
      }

      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let fullResponse = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) {
          setLoading(false);
          toast.success("Molecule generated successfully!");
          const { data } = await axiosInstance.get("/protein/generatednewmolecule");
          const fetchedMolecules = data.molecules || [];
          setMolecules(fetchedMolecules);
          if (fetchedMolecules.length > 0) {
            setSelectedMolecule(fetchedMolecules[fetchedMolecules.length - 1]);
          }
          setFormData({ smilesoffirst: "", smilesofsecond: "", newmoleculetitle: "" });
          break;
        }
        const chunk = decoder.decode(value).split("\n").map(line => line.replace(/^data:\s*/, "").trim()).filter(line => line && line !== "[DONE]").join(" ");
        fullResponse += chunk + " ";
        setRealTimeOutput(fullResponse);
      }
    } catch (err) {
      console.error("Generate error:", err);
      setError(err.message || "Failed to generate molecule.");
      setLoading(false);
    }
  };

  const handleCopySmiles = (smiles) => {
    if (!smiles || smiles === "Not available") {
      toast.error("No SMILES available to copy.");
      return;
    }
    navigator.clipboard.writeText(smiles)
      .then(() => toast.success("SMILES copied!"))
      .catch((err) => toast.error("Copy failed: " + err.message));
  };

  const toggleInfo = (id) => {
    if (expandedInfoId === id) {
      setExpandedInfoId(null);
    } else {
      setExpandedInfoId(id);
    }
  };

  const exportToPDF = (molecule) => {
    const doc = new jsPDF();
    const margin = 20;
    const maxWidth = doc.internal.pageSize.width - 2 * margin;
    let y = margin;

    const addText = (text, size, bold = false) => {
      doc.setFontSize(size);
      doc.setFont(undefined, bold ? "bold" : "normal");
      const lines = doc.splitTextToSize(text, maxWidth);
      lines.forEach(line => {
        if (y + size / 2 > doc.internal.pageSize.height - margin) {
          doc.addPage();
          y = margin;
        }
        doc.text(line, margin, y);
        y += size / 2 + 2;
      });
    };

    addText("Molecule Details", 16, true);
    addText(`Title: ${molecule.newmoleculetitle}`, 12);
    addText(`SMILES: ${molecule.newSmiles || "Not available"}`, 12);
    addText(`IUPAC Name: ${molecule.newIupacName || "Not available"}`, 12);
    addText(`Created: ${new Date(molecule.created).toLocaleString()}`, 12);
    addText("Information:", 12, true);
    addText(molecule.information || "No information available", 10);

    doc.save(`${molecule.newmoleculetitle}_details.pdf`);
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-100">
        <div className="text-lg font-semibold text-gray-700">Verifying Authentication...</div>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gray-100">
        <div className="bg-white p-6 rounded-lg border border-gray-300 text-center">
          <h2 className="text-xl font-bold text-gray-800 mb-4">Access Denied</h2>
          <p className="text-gray-600 mb-4">Please log in to explore Protein Structure Evolution.</p>
          <a
            href="/login"
            className="inline-block px-4 py-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600"
          >
            Login Now
          </a>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-100 py-8 px-4">
      <div className="max-w-5xl mx-auto"> {/* Changed from max-w-3xl to max-w-5xl */}
        <h1 className="text-3xl font-bold text-gray-800 mb-8 text-center">
          New Drug Molecule Generation
          <p className="text-xs p-1 text-blue-700 font-semibold">(Powered by Gemini)</p>

        </h1>

        {error && (
          <div className="mb-6 p-4 bg-red-100 border border-red-400 text-red-700 rounded-lg flex justify-between items-center">
            <span>{error}</span>
            <button onClick={() => setError(null)} className="text-red-900 hover:text-red-700">Dismiss</button>
          </div>
        )}

        <div className="space-y-6">
          {/* Generate New Molecule */}
          <div className="bg-white p-6 rounded-lg border border-gray-300">
            <h2 className="text-xl font-semibold text-gray-800 mb-4">Generate New Molecule</h2>
            <form onSubmit={handleGenerate} className="space-y-4">
              {["newmoleculetitle", "smilesoffirst", "smilesofsecond"].map((field) => (
                <div key={field}>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    {field === "newmoleculetitle" ? "New Molecule Title" : field.replace(/([A-Z])/g, " $1").trim()}
                  </label>
                  <input
                    type="text"
                    name={field}
                    value={formData[field]}
                    onChange={handleInputChange}
                    placeholder={`Enter ${field === "newmoleculetitle" ? "new molecule title" : field.replace(/([A-Z])/g, " $1").trim().toLowerCase()}`}
                    disabled={loading}
                    required
                    className="w-full p-2 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-400"
                  />
                </div>
              ))}
              <button
                type="submit"
                disabled={loading}
                className="w-full py-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600 disabled:bg-gray-400 disabled:cursor-not-allowed"
              >
                {loading ? (
                  <span className="flex items-center justify-center">
                    <svg className="animate-spin h-5 w-5 mr-2 text-white" viewBox="0 0 24 24">
                      <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                      <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v8h8a8 8 0 01-16 0z" />
                    </svg>
                    Generating...
                  </span>
                ) : (
                  "Generate Molecule"
                )}
              </button>
            </form>
          </div>

          {/* Real-Time Output */}
          {realTimeOutput && (
            <div className="bg-white p-6 rounded-lg border border-gray-300">
              <h3 className="text-lg font-semibold text-gray-800 mb-3">Live Output</h3>
              <div className="max-h-48 overflow-y-auto text-gray-700 text-sm p-3 border border-gray-200 rounded-lg">
                {realTimeOutput}
              </div>
            </div>
          )}

          {/* Molecule Details */}
          {selectedMolecule && (
            <div className="bg-white p-6 rounded-lg border border-gray-300">
              <h2 className="text-lg font-semibold text-gray-800 mb-3">Molecule Details</h2>
              <div className="space-y-3 text-sm">
                {[
                  { label: "Title", value: selectedMolecule.newmoleculetitle },
                  { label: "SMILES", value: selectedMolecule.newSmiles || "N/A", copy: true },
                  { label: "IUPAC Name", value: selectedMolecule.newIupacName || "N/A" },
                  { label: "Conversion", value: selectedMolecule.conversionDetails || "N/A" },
                  { label: "Diseases", value: selectedMolecule.potentialDiseases || "N/A" },
                  { label: "Created", value: new Date(selectedMolecule.created).toLocaleString() },
                ].map(({ label, value, copy }) => (
                  <div key={label} className="flex flex-col sm:flex-row sm:items-center">
                    <p className="font-medium text-gray-700 sm:w-28">{label}:</p>
                    <div className="flex items-center mt-1 sm:mt-0">
                      <p className="flex-1 text-gray-600 break-all">{value}</p>
                      {copy && (
                        <FiCopy
                          onClick={() => handleCopySmiles(value)}
                          className="ml-2 text-gray-500 hover:text-blue-600 cursor-pointer"
                        />
                      )}
                    </div>
                  </div>
                ))}
                <div className="flex space-x-2 mt-4">
                  <button
                    onClick={() => toggleInfo(selectedMolecule.id)}
                    className="px-3 py-1 bg-blue-500 text-white rounded-lg hover:bg-blue-600"
                  >
                    {expandedInfoId === selectedMolecule.id ? "Hide Info" : "Show Info"}
                  </button>
                  <button
                    onClick={() => exportToPDF(selectedMolecule)}
                    className="px-3 py-1 bg-teal-500 text-white rounded-lg hover:bg-teal-600"
                  >
                    Export PDF
                  </button>
                </div>
                {expandedInfoId === selectedMolecule.id && (
                  <div className="mt-3 p-3 bg-gray-50 border border-gray-200 rounded-lg text-sm text-gray-700 break-all">
                    {selectedMolecule.information}
                  </div>
                )}
              </div>
            </div>
          )}

          {/* Your Molecules */}
          <div className="bg-white p-6 rounded-lg border border-gray-300">
            <h2 className="text-lg font-semibold text-gray-800 mb-3">Your Molecules</h2>
            {loading && (
              <p className="text-gray-500 flex items-center">
                <svg className="animate-spin h-5 w-5 mr-2 text-gray-500" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v8h8a8 8 0 01-16 0z" />
                </svg>
                Loading...
              </p>
            )}
            {!loading && molecules.length === 0 && (
              <p className="text-gray-500 italic">No molecules generated yet.</p>
            )}
            {!loading && molecules.length > 0 && (
              <ul className="space-y-4 max-h-80 overflow-y-auto">
                {molecules.map((molecule) => (
                  <li key={molecule.id} className="border-b border-gray-200 pb-3 last:border-b-0">
                    <div className="flex flex-col sm:flex-row sm:justify-between sm:items-start">
                      <div
                        onClick={() => setSelectedMolecule(molecule)}
                        className="flex-1 cursor-pointer p-2 rounded-lg hover:bg-gray-100"
                      >
                        <p className="font-medium text-gray-800">{molecule.newmoleculetitle}</p>
                        <div className="flex items-center space-x-2 mt-1">
                          <p className="text-gray-600 text-sm break-all">SMILES: {molecule.newSmiles || "N/A"}</p>
                          <FiCopy
                            onClick={(e) => { e.stopPropagation(); handleCopySmiles(molecule.newSmiles); }}
                            className="text-gray-500 hover:text-blue-600"
                          />
                        </div>
                        <p className="text-gray-500 text-xs mt-1">{new Date(molecule.created).toLocaleString()}</p>
                      </div>
                      <div className="flex space-x-2 mt-2 sm:mt-0">
                        <button
                          onClick={() => toggleInfo(molecule.id)}
                          className="px-3 py-1 bg-blue-500 text-white rounded-lg hover:bg-blue-600"
                        >
                          {expandedInfoId === molecule.id ? "Hide" : "Info"}
                        </button>
                        <button
                          onClick={() => exportToPDF(molecule)}
                          className="px-3 py-1 bg-teal-500 text-white rounded-lg hover:bg-teal-600"
                        >
                          PDF
                        </button>
                      </div>
                    </div>
                    {expandedInfoId === molecule.id && (
                      <div className="mt-2 p-3 bg-gray-50 border border-gray-200 rounded-lg text-sm text-gray-700 break-all">
                        {molecule.information}
                      </div>
                    )}
                  </li>
                ))}
              </ul>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureEvolution;
