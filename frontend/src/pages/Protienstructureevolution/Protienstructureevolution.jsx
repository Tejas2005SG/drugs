import React, { useState, useEffect } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js"; // Adjust path
import { jsPDF } from "jspdf";

// Import a copy icon from react-icons (you'll need to install react-icons if not already installed)
// Run: npm install react-icons
import { FiCopy } from "react-icons/fi";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000";
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const ProteinStructureEvolution = () => {
  const [smilesFirst, setSmilesFirst] = useState("");
  const [smilesSecond, setSmilesSecond] = useState("");
  const [newMoleculeTitle, setNewMoleculeTitle] = useState("");
  const [molecules, setMolecules] = useState([]);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [realTimeOutput, setRealTimeOutput] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [expandedInfo, setExpandedInfo] = useState({});
  const [selectedInfoExpanded, setSelectedInfoExpanded] = useState(false); // For Selected Molecule Details

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
      console.log(`Fetching all molecules for user ${user._id}`);
      const response = await axiosInstance.get("/api/protein/generatednewmolecule");
      const fetchedMolecules = response.data.molecules || [];
      setMolecules(fetchedMolecules);
      if (fetchedMolecules.length > 0 && !selectedMolecule) {
        setSelectedMolecule(fetchedMolecules[0]);
      }
    } catch (err) {
      console.error("Error fetching molecules:", err.response?.data || err.message);
      if (err.response?.status === 404) {
        setMolecules([]);
      } else {
        setError(err.response?.data?.message || "Failed to fetch molecules");
      }
    } finally {
      setLoading(false);
    }
  };

  const handleGenerate = async (e) => {
    e.preventDefault();
    setError(null);
    setLoading(true);
    setRealTimeOutput("");

    if (!user?._id) {
      setError("Please log in to generate a new molecule");
      setLoading(false);
      return;
    }

    try {
      console.log(`Streaming POST to /api/protein/generatenewmolecule/${user._id}`);
      const response = await fetch(`${API_BASE_URL}/api/protein/generatenewmolecule/${user._id}`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "Accept": "text/event-stream",
        },
        credentials: "include",
        body: JSON.stringify({
          smilesoffirst: smilesFirst,
          smilesofsecond: smilesSecond,
          newmoleculetitle: newMoleculeTitle,
        }),
      });

      if (!response.ok) {
        throw new Error(`Server responded with ${response.status}: ${response.statusText}`);
      }

      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let fullResponse = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) {
          setLoading(false);
          toast.success("New molecule generated successfully!");
          fetchAllMolecules();
          break;
        }
        const chunk = decoder.decode(value, { stream: true });
        const cleanedChunk = chunk
          .split("\n")
          .map(line => line.replace(/^data:\s*/, "").trim())
          .filter(line => line && line !== "[DONE]")
          .join(" ");
        fullResponse += cleanedChunk + " ";
        setRealTimeOutput(fullResponse);
      }

      setSmilesFirst("");
      setSmilesSecond("");
      setNewMoleculeTitle("");
    } catch (err) {
      console.error("Frontend error:", err);
      setError(err.message || "Failed to generate molecule");
      setLoading(false);
    }
  };

  // Function to handle copying SMILES to clipboard
  const handleCopySmiles = (smiles) => {
    if (!smiles || smiles === "Not available") {
      toast.error("No SMILES string available to copy");
      return;
    }
    navigator.clipboard.writeText(smiles)
      .then(() => {
        toast.success("SMILES string copied to clipboard!");
      })
      .catch((err) => {
        console.error("Failed to copy SMILES:", err);
        toast.error("Failed to copy SMILES string");
      });
  };

  const toggleInfo = (moleculeId) => {
    setExpandedInfo((prev) => ({
      ...prev,
      [moleculeId]: !prev[moleculeId],
    }));
  };

  const toggleSelectedInfo = () => {
    setSelectedInfoExpanded((prev) => !prev);
  };

  const exportToPDF = (molecule) => {
    const doc = new jsPDF();
    
    // Define page dimensions and margins
    const pageHeight = doc.internal.pageSize.height;
    const pageWidth = doc.internal.pageSize.width;
    const margin = 20;
    const maxWidth = pageWidth - 2 * margin;
    let yPosition = margin; // Starting Y position

    // Set font sizes
    const titleFontSize = 16;
    const labelFontSize = 12;
    const textFontSize = 10;
    const lineHeight = 7; // Approximate line height for font size 10

    // Helper function to add a new page if needed
    const checkAndAddPage = () => {
      if (yPosition + lineHeight > pageHeight - margin) {
        doc.addPage();
        yPosition = margin; // Reset Y position for the new page
      }
    };

    // Add Title
    doc.setFontSize(titleFontSize);
    doc.text("Molecule Details", margin, yPosition);
    yPosition += lineHeight * 2; // Extra space after title

    // Add Molecule Details
    doc.setFontSize(labelFontSize);
    const details = [
      `Title: ${molecule.newmoleculetitle}`,
      `SMILES: ${molecule.newSmiles || "Not available"}`,
      `IUPAC Name: ${molecule.newIupacName || "Not available"}`,
      `Created: ${new Date(molecule.created).toLocaleString()}`,
    ];

    details.forEach((line) => {
      checkAndAddPage();
      doc.text(line, margin, yPosition);
      yPosition += lineHeight;
    });

    // Add Information Section
    yPosition += lineHeight; // Extra space before Information
    checkAndAddPage();
    doc.text("Information:", margin, yPosition);
    yPosition += lineHeight;

    // Switch to smaller font for the information content
    doc.setFontSize(textFontSize);
    const infoText = molecule.information || "No information available";
    const splitInfo = doc.splitTextToSize(infoText, maxWidth); // Split text to fit within page width

    // Add each line of the information, adding new pages as needed
    splitInfo.forEach((line) => {
      checkAndAddPage();
      doc.text(line, margin, yPosition);
      yPosition += lineHeight;
    });

    // Save the PDF
    doc.save(`${molecule.newmoleculetitle}_details.pdf`);
  };

  const handleMoleculeClick = (molecule) => {
    setSelectedMolecule(molecule);
    setSelectedInfoExpanded(false); // Reset the expanded state when a new molecule is selected
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <p className="mt-4 text-gray-600">Verifying authentication...</p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="text-center">
          <p className="text-gray-600">Please log in to access Protein Structure Evolution</p>
          <button
            className="mt-4 px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-500 transition-colors"
            onClick={() => (window.location.href = "/login")}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 py-12">
      <div className="container mx-auto px-4">
        <h1 className="text-4xl font-bold text-blue-600 mb-10 text-center">
          Protein Structure Evolution
        </h1>

        {error && (
          <div className="bg-red-50 border border-red-300 text-red-700 px-4 py 🙂-3 rounded-lg mb-6 max-w-2xl mx-auto">
            <p>{error}</p>
            <button
              className="text-red-700 underline ml-2 hover:text-red-900"
              onClick={() => setError(null)}
            >
              Dismiss
            </button>
          </div>
        )}

        {/* Generate New Molecule Section - Full Width */}
        <div className="bg-white p-6 rounded-xl shadow-lg border border-gray-100 mb-8">
          <h2 className="text-2xl font-semibold text-blue-600 mb-6">
            Generate New Molecule
          </h2>
          <form onSubmit={handleGenerate} className="space-y-5">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Molecule Title
              </label>
              <input
                type="text"
                value={newMoleculeTitle}
                onChange={(e) => setNewMoleculeTitle(e.target.value)}
                placeholder="e.g., My New Molecule"
                disabled={loading}
                required
                className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100 transition-all"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                First SMILES String
              </label>
              <input
                type="text"
                value={smilesFirst}
                onChange={(e) => setSmilesFirst(e.target.value)}
                placeholder="e.g., CCO"
                disabled={loading}
                required
                className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100 transition-all"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Second SMILES String
              </label>
              <input
                type="text"
                value={smilesSecond}
                onChange={(e) => setSmilesSecond(e.target.value)}
                placeholder="e.g., c1ccccc1"
                disabled={loading}
                required
                className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:bg-gray-100 transition-all"
              />
            </div>
            <button
              type="submit"
              disabled={loading}
              className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-500 disabled:bg-gray-400 disabled:cursor-not-allowed transition-colors"
            >
              {loading ? "Generating..." : "Generate Molecule"}
            </button>
          </form>
        </div>

        {/* Real-Time Output Section - Full Width, Increased Height */}
        {realTimeOutput && (
          <div className="bg-blue-50 p-6 rounded-xl border border-blue-200 mb-8">
            <h3 className="text-lg font-semibold text-blue-600 mb-4">Real-Time Output</h3>
            <div className="max-h-64 overflow-y-auto">
              <p className="text-gray-700 leading-relaxed break-words">{realTimeOutput}</p>
            </div>
          </div>
        )}

        {/* Your Generated Molecules and Selected Molecule Details - Side by Side */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
          {/* Your Generated Molecules */}
          <div className="bg-white p-6 rounded-xl shadow-lg border border-gray-100">
            <h2 className="text-2xl font-semibold text-blue-600 mb-6">
              Your Generated Molecules
            </h2>
            {molecules.length === 0 && !loading && (
              <p className="text-gray-600">No molecules generated yet.</p>
            )}
            {molecules.length > 0 && (
              <ul className="space-y-4 max-h-96 overflow-y-auto">
                {molecules.map((molecule) => (
                  <li
                    key={molecule.id}
                    className="border-b border-gray-200 pb-3 rounded-lg p-2"
                  >
                    <div className="flex justify-between items-center">
                      <div
                        className="cursor-pointer hover:bg-blue-50 transition-colors flex-1 p-2 rounded-lg"
                        onClick={() => handleMoleculeClick(molecule)}
                      >
                        <p className="text-gray-800 font-medium">
                          <strong>Title:</strong> {molecule.newmoleculetitle}
                        </p>
                        <div className="flex items-center space-x-2">
                          <p className="text-gray-700 truncate">
                            <strong>SMILES:</strong> {molecule.newSmiles || "Not available"}
                          </p>
                          <button
                            onClick={(e) => {
                              e.stopPropagation(); // Prevent triggering the molecule click
                              handleCopySmiles(molecule.newSmiles);
                            }}
                            className="text-gray-500 hover:text-blue-600 transition-colors"
                            title="Copy SMILES"
                          >
                            <FiCopy size={16} />
                          </button>
                        </div>
                        <p className="text-gray-500 text-sm">
                          <strong>Created:</strong>{" "}
                          {new Date(molecule.created).toLocaleString()}
                        </p>
                      </div>
                      <div className="flex space-x-2">
                        <button
                          onClick={() => toggleInfo(molecule.id)}
                          className="px-3 py-1 bg-blue-500 text-white rounded-lg hover:bg-blue-400 transition-colors text-sm"
                        >
                          {expandedInfo[molecule.id] ? "Hide Info" : "Information"}
                        </button>
                        <button
                          onClick={() => exportToPDF(molecule)}
                          className="px-3 py-1 bg-blue-600 text-white rounded-lg hover:bg-blue-500 transition-colors text-sm"
                        >
                          Export PDF
                        </button>
                      </div>
                    </div>
                    {expandedInfo[molecule.id] && (
                      <div className="mt-2 p-3 bg-blue-50 rounded-lg">
                        <p className="text-gray-700 leading-relaxed break-words">
                          {expandedInfo[molecule.id]
                            ? molecule.information
                            : molecule.information.length > 100
                            ? `${molecule.information.substring(0, 100)}...`
                            : molecule.information}
                        </p>
                        {molecule.information.length > 100 && (
                          <button
                            onClick={() => toggleInfo(molecule.id)}
                            className="text-blue-600 underline hover:text-blue-500 mt-1 text-sm"
                          >
                            {expandedInfo[molecule.id] ? "Read Less" : "Read More"}
                          </button>
                        )}
                      </div>
                    )}
                  </li>
                ))}
              </ul>
            )}
          </div>

          {/* Selected Molecule Details */}
          {selectedMolecule && (
            <div className="bg-white p-6 rounded-xl shadow-lg border border-gray-100">
              <h2 className="text-2xl font-semibold text-blue-600 mb-6">
                Selected Molecule Details
              </h2>
              <div className="space-y-3">
                <p className="text-gray-800">
                  <strong className="text-blue-600">Title:</strong> {selectedMolecule.newmoleculetitle}
                </p>
                <div className="flex items-center space-x-2">
                  <p className="text-gray-800">
                    <strong className="text-blue-600">SMILES:</strong>{" "}
                    {selectedMolecule.newSmiles || "Not available"}
                  </p>
                  <button
                    onClick={() => handleCopySmiles(selectedMolecule.newSmiles)}
                    className="text-gray-500 hover:text-blue-600 transition-colors"
                    title="Copy SMILES"
                  >
                    <FiCopy size={16} />
                  </button>
                </div>
                <p className="text-gray-800">
                  <strong className="text-blue-600">IUPAC Name:</strong>{" "}
                  {selectedMolecule.newIupacName || "Not available"}
                </p>
                <p className="text-gray-800">
                  <strong className="text-blue-600">Conversion Details:</strong>{" "}
                  {selectedMolecule.conversionDetails || "Not available"}
                </p>
                <p className="text-gray-800">
                  <strong className="text-blue-600">Potential Diseases:</strong>{" "}
                  {selectedMolecule.potentialDiseases || "Not available"}
                </p>
                <p className="text-gray-500 text-sm">
                  <strong className="text-blue-600">Created:</strong>{" "}
                  {new Date(selectedMolecule.created).toLocaleString()}
                </p>
                {/* Add Information and Export PDF Buttons */}
                <div className="flex space-x-2 mt-4">
                  <button
                    onClick={toggleSelectedInfo}
                    className="px-3 py-1 bg-blue-500 text-white rounded-lg hover:bg-blue-400 transition-colors text-sm"
                  >
                    {selectedInfoExpanded ? "Hide Info" : "Information"}
                  </button>
                  <button
                    onClick={() => exportToPDF(selectedMolecule)}
                    className="px-3 py-1 bg-blue-600 text-white rounded-lg hover:bg-blue-500 transition-colors text-sm"
                  >
                    Export PDF
                  </button>
                </div>
                {selectedInfoExpanded && (
                  <div className="mt-2 p-3 bg-blue-50 rounded-lg">
                    <p className="text-gray-700 leading-relaxed break-words">
                      {selectedInfoExpanded
                        ? selectedMolecule.information
                        : selectedMolecule.information.length > 100
                        ? `${selectedMolecule.information.substring(0, 100)}...`
                        : selectedMolecule.information}
                    </p>
                    {selectedMolecule.information.length > 100 && (
                      <button
                        onClick={toggleSelectedInfo}
                        className="text-blue-600 underline hover:text-blue-500 mt-1 text-sm"
                      >
                        {selectedInfoExpanded ? "Read Less" : "Read More"}
                      </button>
                    )}
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureEvolution;