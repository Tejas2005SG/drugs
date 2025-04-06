import React, { useState, useEffect } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";

// Utility to truncate long SMILES strings
const truncateSmiles = (smiles, maxLength = 20) => {
  if (smiles.length <= maxLength) return smiles;
  return `${smiles.substring(0, maxLength)}...`;
};

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const AINamingSuggestion = () => {
  const [activeTab, setActiveTab] = useState("generate");
  const [molecules, setMolecules] = useState([]);
  const [selectedTitle, setSelectedTitle] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [suggestedNames, setSuggestedNames] = useState([]);
  const [savedNames, setSavedNames] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [fallbackMessage, setFallbackMessage] = useState(null);

  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const initialize = async () => {
      await checkAuth();
      if (!user) {
        setError("Authentication failed. Please log in.");
        return;
      }
      await fetchAllMolecules();
      await fetchSavedNames();
    };
    initialize();
  }, [checkAuth]);

  const fetchAllMolecules = async () => {
    if (!user?._id) return;

    setLoading(true);
    try {
      const response = await axiosInstance.get("/protein/generatednewmolecule");
      const fetchedMolecules = response.data.molecules || [];
      setMolecules(fetchedMolecules);
      if (fetchedMolecules.length > 0 && !selectedTitle && !selectedSmiles) {
        setSelectedTitle(fetchedMolecules[0].newmoleculetitle);
        setSelectedSmiles(fetchedMolecules[0].newSmiles);
      }
    } catch (err) {
      console.error("Error fetching molecules:", err);
      setError(err.response?.data?.message || "Failed to fetch molecules");
    } finally {
      setLoading(false);
    }
  };

  const fetchSavedNames = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/protein/saved-drug-names");
      setSavedNames(response.data.drugNames || []);
    } catch (err) {
      console.error("Error fetching saved names:", err);
      setError(err.response?.data?.message || "Failed to fetch saved names");
    }
  };

  const checkIfNameExists = async (title, smiles) => {
    try {
      const response = await axiosInstance.get("/protein/check-saved-drug-name", {
        params: { moleculeTitle: title, smiles },
      });
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved name:", err);
      return false;
    }
  };

  const handleGenerateName = async () => {
    if (!selectedTitle || !selectedSmiles) {
      toast.error("Please select both a title and SMILES string");
      return;
    }

    const nameExists = await checkIfNameExists(selectedTitle, selectedSmiles);
    if (nameExists) {
      const savedName = savedNames.find((n) => n.moleculeTitle === selectedTitle && n.smiles === selectedSmiles);
      if (savedName?.status === "accepted") {
        toast("An accepted drug name already exists for this molecule.", { type: "info" });
        setActiveTab("saved");
        return;
      }
    }

    setLoading(true);
    setError(null);
    setSuggestedNames([]);
    setFallbackMessage(null);

    try {
      const response = await axiosInstance.post(`/protein/generate-drug-name/${user._id}`, {
        moleculeTitle: selectedTitle,
        smiles: selectedSmiles,
      });

      if (response.status === 409) {
        toast("An accepted drug name already exists. Redirecting to Saved Names.", { type: "info" });
        setActiveTab("saved");
        return;
      }

      setSuggestedNames(response.data.allCandidates);
      setFallbackMessage(response.data.fallback);
      toast.success("Drug names generated successfully!");
    } catch (err) {
      console.error("Error generating drug name:", err);
      setError(err.response?.data?.message || "Failed to generate drug name");
      toast.error("Failed to generate drug name");
    } finally {
      setLoading(false);
    }
  };

  const handleAcceptName = async (candidate) => {
    const confirm = window.confirm(
      `By accepting "${candidate.name}", this name will become the final title for this molecule across the database and cannot be changed later. Proceed?`
    );
    if (!confirm) return;

    setLoading(true);
    try {
      const response = await axiosInstance.post(`/protein/accept-drug-name/${user._id}`, {
        moleculeTitle: selectedTitle,
        smiles: selectedSmiles,
        selectedName: candidate.name,
        rationale: candidate.rationale,
        compliance: candidate.compliance,
      });

      toast.success("Drug name accepted and molecule title updated!");
      setSuggestedNames([]);
      setSelectedTitle(candidate.name); // Update local state
      await fetchAllMolecules(); // Refresh molecule list
      await fetchSavedNames(); // Refresh saved names
    } catch (err) {
      console.error("Error accepting drug name:", err);
      setError(err.response?.data?.message || "Failed to accept drug name");
      toast.error("Failed to accept drug name");
    } finally {
      setLoading(false);
    }
  };

  const handleRejectName = () => {
    setSuggestedNames([]);
    toast.success("Generated names rejected. You can generate new ones.");
  };

  const handleTabChange = async (tab) => {
    if (activeTab === "generate" && suggestedNames.length > 0) {
      // Save as pending if switching tabs without accepting/rejecting
      try {
        await axiosInstance.post(`/protein/save-pending-drug-name/${user._id}`, {
          moleculeTitle: selectedTitle,
          smiles: selectedSmiles,
          candidates: suggestedNames,
        });
        toast("Generated names saved as pending.", { type: "info" });
      } catch (err) {
        console.error("Error saving pending drug name:", err);
        toast.error("Failed to save pending drug name");
      }
    }
    setActiveTab(tab);
    setSuggestedNames([]);
    setError(null);
    setFallbackMessage(null);
  };

  const getComplianceText = (compliance) => {
    if (typeof compliance === "string") return compliance;
    if (compliance && typeof compliance === "object") {
      return compliance.status || JSON.stringify(compliance);
    }
    return "Unknown";
  };

  const isComplianceFail = (compliance) => {
    const text = getComplianceText(compliance).toLowerCase();
    return text.includes("fail");
  };

  if (checkingAuth) {
    return <div className="text-center py-10 text-gray-600">Verifying authentication...</div>;
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-100">
        <div className="text-center p-6 bg-white rounded-lg shadow-lg max-w-md w-full">
          <p className="text-gray-600 mb-4">Please log in to access AI Naming Suggestion</p>
          <button
            className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors w-full sm:w-auto"
            onClick={() => (window.location.href = "/login")}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-100 py-6 px-4 sm:py-8 sm:px-6 lg:py-12 lg:px-8">
      <div className="max-w-5xl mx-auto">
        <h1 className="text-3xl sm:text-4xl font-bold text-blue-700 mb-6 sm:mb-10 text-center">
          AI Drug Naming Suggestion
          <p className="text-xs sm:text-sm p-1 text-blue-700 font-semibold">(Powered by Gemini)</p>
        </h1>

        <div className="flex flex-col sm:flex-row justify-center mb-4 sm:mb-8 space-y-2 sm:space-y-0 sm:space-x-4">
          {["generate", "saved"].map((tab) => (
            <button
              key={tab}
              className={`px-4 sm:px-6 py-2 rounded-lg font-medium text-sm sm:text-base transition-all duration-300 ${
                activeTab === tab
                  ? "bg-blue-600 text-white shadow-md"
                  : "bg-gray-200 text-gray-700 hover:bg-gray-300"
              }`}
              onClick={() => handleTabChange(tab)}
            >
              {tab === "generate" && "Generate Drug Name"}
              {tab === "saved" && "Saved Names"}
            </button>
          ))}
        </div>

        <div className="bg-white p-4 sm:p-6 lg:p-8 rounded-xl shadow-lg border border-gray-200">
          {activeTab === "generate" && (
            <>
              <h2 className="text-xl sm:text-2xl font-semibold text-blue-700 mb-4 sm:mb-6">Generate Drug Name</h2>

              <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 sm:gap-6 mb-4 sm:mb-6">
                <div>
                  <label className="block text-sm sm:text-base font-medium text-gray-700 mb-2">
                    Select Molecule Title
                  </label>
                  <select
                    value={selectedTitle}
                    onChange={(e) => setSelectedTitle(e.target.value)}
                    className="w-full p-2 sm:p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all text-sm sm:text-base"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No titles available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newmoleculetitle))].map((title) => (
                        <option key={title} value={title}>
                          {title}
                        </option>
                      ))
                    )}
                  </select>
                </div>
                <div>
                  <label className="block text-sm sm:text-base font-medium text-gray-700 mb-2">
                    Select SMILES String
                  </label>
                  <select
                    value={selectedSmiles}
                    onChange={(e) => setSelectedSmiles(e.target.value)}
                    className="w-full p-2 sm:p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all text-sm sm:text-base"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newSmiles))].map((smiles) => (
                        <option key={smiles} value={smiles}>
                          {truncateSmiles(smiles)}
                        </option>
                      ))
                    )}
                  </select>
                </div>
              </div>

              <button
                onClick={handleGenerateName}
                disabled={loading || !selectedTitle || !selectedSmiles}
                className="w-full py-2 sm:py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300 text-sm sm:text-base"
              >
                {loading ? "Generating Names..." : "Start Prediction"}
              </button>

              {error && (
                <div className="mt-4 sm:mt-6 bg-red-50 border border-red-300 text-red-700 px-4 py-3 rounded-lg flex justify-between items-center">
                  <p className="text-sm sm:text-base">{error}</p>
                  <button
                    className="text-red-700 underline hover:text-red-900 text-sm sm:text-base"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </div>
              )}

              {suggestedNames.length > 0 && (
                <div className="mt-6 sm:mt-8 bg-blue-50 p-4 sm:p-6 rounded-xl border border-blue-200 animate-fadeIn">
                  <h3 className="text-lg sm:text-xl font-semibold text-blue-700 mb-4">Suggested Drug Names</h3>
                  <div className="space-y-4 sm:space-y-6">
                    {suggestedNames.map((candidate) => (
                      <div
                        key={candidate.rank}
                        className="border-b border-gray-200 pb-2 sm:pb-4 last:border-b-0"
                      >
                        <div className="flex items-center mb-2">
                          <span className="text-sm sm:text-base font-medium text-gray-600 mr-2">
                            Rank {candidate.rank}:
                          </span>
                          <p className="text-lg sm:text-xl font-bold text-blue-600">{candidate.name}</p>
                        </div>
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Structural Rationale</p>
                          <p className="text-gray-700">{candidate.rationale}</p>
                        </div>
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Compliance Status</p>
                          {/* <p
                            className={`text-gray-700 ${
                              isComplianceFail(candidate.compliance) ? "text-red-600" : "text-green-600"
                            }`}
                          >
                            {getComplianceText(candidate.compliance)}
                          </p> */}
                        </div>
                        <div className="mt-2 sm:mt-4 flex space-x-2 sm:space-x-4">
                          <button
                            onClick={() => handleAcceptName(candidate)}
                            className="px-3 sm:px-4 py-1 sm:py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors text-sm sm:text-base"
                            disabled={loading}
                          >
                            Accept
                          </button>
                          <button
                            onClick={handleRejectName}
                            className="px-3 sm:px-4 py-1 sm:py-2 bg-red-600 text-white rounded-lg hover:bg-red-700 transition-colors text-sm sm:text-base"
                            disabled={loading}
                          >
                            Reject
                          </button>
                        </div>
                      </div>
                    ))}
                  </div>
                  {fallbackMessage && (
                    <div className="mt-4 sm:mt-6 bg-yellow-50 p-3 sm:p-4 rounded-lg border border-yellow-200">
                      <p className="text-yellow-700 font-medium text-sm sm:text-base">Fallback Notice:</p>
                      <p className="text-yellow-600 text-sm sm:text-base">{fallbackMessage}</p>
                    </div>
                  )}
                </div>
              )}
            </>
          )}

          {activeTab === "saved" && (
            <>
              <h2 className="text-xl sm:text-2xl font-semibold text-blue-700 mb-4 sm:mb-6">Saved Drug Names</h2>
              {savedNames.length > 0 ? (
                <div className="space-y-4 sm:space-y-6">
                  {savedNames.map((drugName) => (
                    <div
                      key={drugName._id}
                      className="bg-blue-50 p-4 sm:p-6 rounded-xl border border-blue-200 transition-all duration-200 hover:shadow-md"
                    >
                      <h3 className="text-lg sm:text-xl font-semibold text-gray-800 mb-2 sm:mb-4">
                        {drugName.moleculeTitle} (SMILES: {truncateSmiles(drugName.smiles)})
                      </h3>
                      <div className="space-y-2 sm:space-y-3">
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Suggested Name</p>
                          <p className="text-lg sm:text-xl font-bold text-blue-600">{drugName.suggestedName}</p>
                        </div>
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Naming Details</p>
                          <p className="text-gray-700">{drugName.namingDetails}</p>
                        </div>
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Status</p>
                          <p
                            className={`text-gray-700 ${
                              drugName.status === "accepted" ? "text-green-600" : "text-yellow-600"
                            } text-sm sm:text-base`}
                          >
                            {drugName.status.charAt(0).toUpperCase() + drugName.status.slice(1)}
                          </p>
                        </div>
                        <div>
                          <p className="text-sm sm:text-base text-gray-600">Created</p>
                          <p className="text-gray-700 text-sm sm:text-base">
                            {new Date(drugName.createdAt).toLocaleString()}
                          </p>
                        </div>
                        {drugName.status === "pending" && (
                          <div className="mt-2 sm:mt-4 flex space-x-2 sm:space-x-4">
                            <button
                              onClick={() =>
                                handleAcceptName({
                                  name: drugName.suggestedName,
                                  rationale: drugName.namingDetails.split(" | Compliance: ")[0],
                                  compliance: drugName.namingDetails.split(" | Compliance: ")[1],
                                })
                              }
                              className="px-3 sm:px-4 py-1 sm:py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors text-sm sm:text-base"
                              disabled={loading}
                            >
                              Accept
                            </button>
                            <button
                              onClick={async () => {
                                await axiosInstance.delete(`/protein/delete-drug-name/${drugName._id}`);
                                toast.success("Pending name rejected and removed.");
                                await fetchSavedNames();
                              }}
                              className="px-3 sm:px-4 py-1 sm:py-2 bg-red-600 text-white rounded-lg hover:bg-red-700 transition-colors text-sm sm:text-base"
                              disabled={loading}
                            >
                              Reject
                            </button>
                          </div>
                        )}
                      </div>
                    </div>
                  ))}
                </div>
              ) : (
                <p className="text-gray-600 text-center text-sm sm:text-base">No saved drug names found.</p>
              )}
            </>
          )}
        </div>
      </div>

      <style jsx>{`
        .animate-fadeIn {
          animation: fadeIn 0.5s ease-in-out;
        }
        @keyframes fadeIn {
          from {
            opacity: 0;
            transform: translateY(10px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }
      `}</style>
    </div>
  );
};

export default AINamingSuggestion;