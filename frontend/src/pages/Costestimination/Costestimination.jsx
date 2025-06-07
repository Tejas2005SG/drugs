"use client";

import { useState, useEffect, useRef } from "react";
import { useAuthStore } from "../../Store/auth.store.js";
import { postCostEstimation, getCostEstimations } from "../../api/costestimination.jsx";
import axios from "axios";
import { AlertCircle, DollarSign, Clock, Database, X, Info, RefreshCw, LogIn, FileText, ChevronDown, ChevronUp, Download } from "lucide-react";
import { jsPDF } from "jspdf";
import domtoimage from "dom-to-image";
import Loader from "../../components/Loader.jsx"
const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const CostEstimationForm = () => {
  const [smiles, setSmiles] = useState("");
  const [molecules, setMolecules] = useState([]);
  const [result, setResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [historyLoading, setHistoryLoading] = useState(false);
  const [isResultInfoOpen, setIsResultInfoOpen] = useState(false);
  const [openHistoryItems, setOpenHistoryItems] = useState({});

  const { user, checkAuth, checkingAuth } = useAuthStore();
  const infoRef = useRef(null);

  // Fetch all generated molecules
  const fetchAllMolecules = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/protein/generatednewmolecule");
      const fetchedMolecules = response.data.molecules || [];
      setMolecules(fetchedMolecules);
      if (fetchedMolecules.length > 0 && !smiles) {
        setSmiles(fetchedMolecules[0].newSmiles);
      }
    } catch (err) {
      console.error("Error fetching molecules:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch molecules");
      setMolecules([]);
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!user?._id) {
      setError("Please log in to estimate costs");
      return;
    }
    if (!smiles) {
      setError("Please select a SMILES string");
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const data = await postCostEstimation(smiles);
      setResult(data.data);
      await fetchHistory();
    } catch (error) {
      setError("Failed to estimate cost. Please try again.");
      console.error("Error estimating cost:", error);
    } finally {
      setLoading(false);
    }
  };

  const fetchHistory = async () => {
    if (!user?._id) {
      setError("Please log in to view history");
      return;
    }
    setHistoryLoading(true);
    try {
      const data = await getCostEstimations(user._id);
      const validHistory = Array.isArray(data.data) ? data.data.filter((item) => item && typeof item === "object") : [];
      setHistory(validHistory);
    } catch (error) {
      setError("Failed to fetch history.");
      console.error("Error fetching history:", error);
    } finally {
      setHistoryLoading(false);
    }
  };

  const toggleResultInfo = () => {
    setIsResultInfoOpen(!isResultInfoOpen);
  };

  const toggleHistoryItem = (id) => {
    setOpenHistoryItems((prev) => ({
      ...prev,
      [id]: !prev[id],
    }));
  };

  const renderInformation = (info) => {
    if (!info) return <p className="text-text-secondary">No additional information available</p>;

    const lines = info.split("\n").filter((line) => line.trim() !== "");
    let introParagraph = "";
    const sections = [];
    let currentSection = null;

    lines.forEach((line, index) => {
      if (index === 0 && !line.match(/^\d+\./)) {
        introParagraph = line.trim();
      } else if (line.match(/^\d+\.\s/)) {
        if (currentSection) sections.push(currentSection);
        currentSection = { title: line.trim(), bullets: [] };
      } else if (line.trim().startsWith("-")) {
        if (currentSection) currentSection.bullets.push(line.trim().replace("-", "").trim());
      } else if (currentSection) {
        currentSection.bullets.push(line.trim());
      }
    });
    if (currentSection) sections.push(currentSection);

    return (
      <div className="font-body">
        {introParagraph && <p className="mb-2 text-text-primary">{introParagraph}</p>}
        {sections.length > 0 && (
          <ul className="list-disc pl-5 space-y-2">
            {sections.map((section, index) => (
              <li key={index} className="font-semibold text-accent">
                {section.title}
                {section.bullets.length > 0 && (
                  <ul className="list-disc pl-5 mt-1 space-y-1">
                    {section.bullets.map((bullet, bulletIndex) => (
                      <li key={bulletIndex} className="text-text-primary">{bullet}</li>
                    ))}
                  </ul>
                )}
              </li>
            ))}
          </ul>
        )}
      </div>
    );
  };

  const exportToPDF = async () => {
    const pdf = new jsPDF("p", "mm", "a4");
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 10;
    let yPosition = margin;

    pdf.setFontSize(16);
    pdf.text("Drug Cost Estimation Report", margin, yPosition);
    yPosition += 10;

    pdf.setFontSize(12);
    pdf.text("SMILES:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result.smiles || "N/A", margin, yPosition);
    yPosition += 8;

    pdf.setFontSize(12);
    pdf.text("Estimated Cost:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result.estimatedcost || "N/A", margin, yPosition);
    yPosition += 10;

    if (infoRef.current) {
      pdf.setFontSize(12);
      pdf.text("Information:", margin, yPosition);
      yPosition += 8;

      const infoImgData = await domtoimage.toPng(infoRef.current, { quality: 1 });
      const infoImgProps = pdf.getImageProperties(infoImgData);
      const infoImgWidth = pageWidth - 2 * margin;
      let infoImgHeight = (infoImgProps.height * infoImgWidth) / infoImgProps.width;

      let remainingHeight = infoImgHeight;
      let yOffset = 0;

      while (remainingHeight > 0) {
        const spaceLeft = pageHeight - yPosition - margin;
        const heightToRender = Math.min(remainingHeight, spaceLeft);

        const tempCanvas = document.createElement("canvas");
        const tempCtx = tempCanvas.getContext("2d");
        const img = new Image();
        img.src = infoImgData;
        await new Promise((resolve) => {
          img.onload = resolve;
        });
        tempCanvas.width = img.width;
        tempCanvas.height = (heightToRender / infoImgHeight) * img.height;
        tempCtx.drawImage(img, 0, yOffset, img.width, tempCanvas.height, 0, 0, img.width, tempCanvas.height);
        const croppedImgData = tempCanvas.toDataURL("image/png");

        pdf.addImage(croppedImgData, "PNG", margin, yPosition, infoImgWidth, heightToRender);
        yPosition += heightToRender;
        remainingHeight -= heightToRender;
        yOffset += tempCanvas.height;

        if (remainingHeight > 0) {
          pdf.addPage();
          yPosition = margin;
        }
      }
    }

    pdf.save("drug-cost-estimation-report.pdf");
  };

  useEffect(() => {
    const initialize = async () => {
      console.log("Starting initialization - checkingAuth:", checkingAuth, "user:", user);
      setLoading(true);
      try {
        const checkedUser = await checkAuth();
        console.log("After checkAuth - checkingAuth:", useAuthStore.getState().checkingAuth, "user:", useAuthStore.getState().user);
        
        if (!checkedUser || !useAuthStore.getState().user?._id) {
          setError("Authentication failed. Please log in.");
          setLoading(false);
          return;
        }

        await fetchAllMolecules();
        await fetchHistory();
      } catch (err) {
        console.error("Initialization error:", err);
        setError("Failed to verify authentication. Please try refreshing the page or logging in again.");
      } finally {
        setLoading(false);
      }
    };

    initialize();
  }, [checkAuth]);

  // if (checkingAuth || loading) {
  //   return (
  //     <div className="flex flex-col items-center justify-center h-screen bg-primary">
  //       <Loader />
  //     </div>
  //   );
  // }

  if (!user || !user._id) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-primary">
        <div className="text-center max-w-md p-8 bg-secondary rounded-xl shadow-lg border border-secondary">
          <LogIn size={48} className="mx-auto text-accent mb-4" />
          <h2 className="text-2xl font-bold text-text-primary mb-3 font-heading">Authentication Required</h2>
          <p className="text-text-secondary mb-6 font-body">
            Please log in to access the Drug Cost Estimator tool and view your estimation history.
          </p>
          <button
            className="w-full px-4 py-3 bg-accent text-primary font-medium rounded-lg hover:bg-accent/90 transition-colors duration-200 flex items-center justify-center font-heading"
            onClick={() => (window.location.href = "/login")}
          >
            <LogIn size={18} className="mr-2" />
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-primary py-12">
      <div className="container mx-auto px-4 max-w-4xl">
        <div className="mb-10 text-center">
          <h1 className="text-4xl font-bold text-text-primary mb-2 font-heading">
            Drug Cost Estimator
            <p className="text-xs p-1 text-accent-secondary font-semibold font-label">(Powered by Gemini)</p>
          </h1>
          <p className="text-text-secondary max-w-2xl mx-auto font-body">
            Select a SMILES string from your generated molecules to estimate the cost of drug synthesis and production.
          </p>
        </div>

        <div className="space-y-8">
          {/* Form Section */}
          <div className="bg-secondary p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg border border-secondary">
            <div className="flex items-center mb-6">
              <Database className="h-6 w-6 text-accent mr-2" />
              <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimate New Cost</h2>
            </div>

            <form onSubmit={handleSubmit} className="space-y-4">
              <div>
                <label htmlFor="smiles" className="block text-sm font-medium text-text-secondary mb-1 font-body">
                  Select SMILES String
                </label>
                <select
                  id="smiles"
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  className="w-full px-4 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent transition-all duration-200 font-body"
                  disabled={loading || molecules.length === 0}
                >
                  {molecules.length === 0 ? (
                    <option value="" className="bg-primary">No SMILES available</option>
                  ) : (
                    [...new Set(molecules.map((m) => m.newSmiles))].map((smilesOption) => (
                      <option key={smilesOption} value={smilesOption} className="bg-primary">
                        {smilesOption}
                      </option>
                    ))
                  )}
                </select>
                <p className="mt-1 text-xs text-text-secondary font-body">
                  Select a SMILES string from your generated molecules
                </p>
              </div>

              <button
                type="submit"
                disabled={loading || !smiles}
                className={`w-full px-4 py-3 text-primary rounded-lg transition-all duration-200 flex items-center justify-center font-heading ${
                  loading || !smiles ? "bg-gray-400 cursor-not-allowed" : "bg-accent hover:bg-accent/90"
                }`}
              >
                {loading ? (
                  <>
                    <RefreshCw size={20} className="animate-spin mr-2" />
                    Estimating...
                  </>
                ) : (
                  <>
                    <DollarSign size={20} className="mr-2" />
                    Estimate Cost
                  </>
                )}
              </button>
            </form>

            {error && (
              <div className="mt-6 bg-error/10 border border-error text-text-primary px-4 py-3 rounded-lg flex justify-between items-center font-body">
                <p>{error}</p>
                <button
                  className="text-error underline hover:text-error/80"
                  onClick={() => setError(null)}
                >
                  Dismiss
                </button>
              </div>
            )}

            {result && (
              <div className="mt-8 animate-fadeIn">
                <div className="flex items-center mb-4">
                  <DollarSign className="h-5 w-5 text-success mr-2" />
                  <h3 className="text-xl font-semibold text-text-primary font-heading">Estimation Result</h3>
                </div>
                <div className="bg-gradient-to-r from-primary/20 to-secondary/20 p-6 rounded-lg border border-accent/20">
                  <div className="space-y-3">
                    <div>
                      <p className="text-sm text-text-secondary mb-1 font-body">SMILES</p>
                      <p className="font-mono text-sm bg-primary p-2 rounded border border-secondary overflow-x-auto font-code">
                        {result.smiles || "N/A"}
                      </p>
                    </div>

                    <div>
                      <p className="text-sm text-text-secondary mb-1 font-body">Estimated Cost</p>
                      <p className="text-3xl font-bold text-success font-heading">
                        {result.estimatedcost || "N/A"}
                      </p>
                    </div>

                    <div className="space-y-2">
                      <button
                        onClick={toggleResultInfo}
                        className="flex items-center justify-between w-full p-3 bg-accent-secondary/10 hover:bg-accent-secondary/20 rounded-lg transition-colors duration-200 border border-accent-secondary/30 focus:outline-none focus:ring-2 focus:ring-accent-secondary/50 focus:ring-offset-2 focus:ring-offset-secondary font-body"
                      >
                        <div className="flex items-center space-x-2">
                          <Info className="h-5 w-5 text-accent-secondary" />
                          <span className="text-base font-semibold text-text-primary">
                            Detailed Analysis <span className="text-sm text-accent font-label">(powered by Gemini)</span>
                          </span>
                        </div>
                        {isResultInfoOpen ? (
                          <ChevronUp className="h-5 w-5 text-accent-secondary transform transition-transform duration-300" />
                        ) : (
                          <ChevronDown className="h-5 w-5 text-accent-secondary transform transition-transform duration-300" />
                        )}
                      </button>

                      {isResultInfoOpen && (
                        <div className="space-y-4">
                          <div
                            ref={infoRef}
                            className="p-4 bg-secondary rounded-lg shadow-lg border border-secondary transition-all duration-300"
                          >
                            <div className="prose prose-sm text-text-primary max-w-none font-body">
                              {renderInformation(result.information)}
                            </div>
                          </div>

                          <div className="flex justify-end">
                            <button
                              onClick={exportToPDF}
                              className="flex items-center px-4 py-2.5 bg-gradient-to-br from-success to-success/90 text-primary font-medium rounded-md hover:from-success/90 hover:to-success transition-all duration-200 transform hover:scale-[1.02] shadow-sm hover:shadow-md focus:outline-none focus:ring-2 focus:ring-success focus:ring-offset-2 focus:ring-offset-secondary font-heading"
                            >
                              <Download size={18} className="mr-2 transform transition-transform hover:-translate-y-0.5" />
                              Export Report
                            </button>
                          </div>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>

          {/* History Section */}
          <div className="bg-secondary p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg border border-secondary">
            <div className="flex items-center justify-between mb-6">
              <div className="flex items-center">
                <Clock className="h-6 w-6 text-accent mr-2" />
                <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimation History</h2>
              </div>
              <button
                onClick={fetchHistory}
                disabled={historyLoading}
                className="text-accent hover:text-accent/80 transition-colors duration-200"
                title="Refresh history"
              >
                <RefreshCw size={20} className={historyLoading ? "animate-spin" : ""} />
              </button>
            </div>

            {historyLoading ? (
              <div className="flex flex-col items-center justify-center py-12">
                <RefreshCw size={24} className="animate-spin mb-4 text-accent" />
                <p className="text-text-secondary font-body">Loading history...</p>
              </div>
            ) : history.length > 0 ? (
              <div className="space-y-4 max-h-[500px] overflow-y-auto pr-2 custom-scrollbar">
                {history.map((item) => (
                  <div
                    key={item._id}
                    className="border border-secondary p-4 rounded-lg hover:border-accent-secondary transition-all duration-200 hover:shadow-sm"
                  >
                    <div className="flex justify-between items-start mb-2">
                      <div className="flex items-center">
                        <DollarSign size={16} className="text-success mr-1" />
                        <span className="font-bold text-lg text-success font-heading">
                          {item.estimatedcost || "N/A"}
                        </span>
                      </div>
                      <span className="text-xs text-text-secondary bg-primary px-2 py-1 rounded-full font-body">
                        {item.created ? new Date(item.created).toLocaleDateString() : "N/A"}
                      </span>
                    </div>

                    <div className="mt-2">
                      <p className="text-xs text-text-secondary mb-1 font-body">SMILES</p>
                      <p className="font-mono text-sm bg-primary p-2 rounded border border-secondary overflow-x-auto font-code">
                        {item.smiles || "N/A"}
                      </p>
                    </div>

                    {item.information && (
                      <div className="mt-3 space-y-2">
                        <button
                          onClick={() => toggleHistoryItem(item._id)}
                          className="flex items-center justify-between w-full p-2.5 bg-accent-secondary/10 hover:bg-accent-secondary/20 rounded-md transition-colors duration-200 border border-accent-secondary/30 focus:outline-none focus:ring-2 focus:ring-accent-secondary/50 focus:ring-offset-1 focus:ring-offset-secondary font-body"
                        >
                          <div className="flex items-center space-x-2">
                            <Info className="h-4 w-4 text-accent-secondary" />
                            <span className="text-xs font-medium text-text-primary uppercase tracking-wide">
                              Details <span className="text-xs text-accent font-label">(powered by Gemini)</span>
                            </span>
                          </div>
                          {openHistoryItems[item._id] ? (
                            <ChevronUp className="h-4 w-4 text-accent-secondary transform transition-transform duration-300" />
                          ) : (
                            <ChevronDown className="h-4 w-4 text-accent-secondary transform transition-transform duration-300" />
                          )}
                        </button>

                        {openHistoryItems[item._id] && (
                          <div className="space-y-3">
                            <div className="p-3 bg-secondary rounded-md shadow-sm border border-secondary animate-fade-in">
                              <div className="prose prose-sm text-text-primary max-w-none font-body">
                                {renderInformation(item.information)}
                              </div>
                            </div>
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            ) : (
              <div className="flex flex-col items-center justify-center py-12 text-center">
                <Database size={32} className="text-text-secondary mb-4" />
                <p className="text-text-secondary mb-2 font-body">No previous estimations found.</p>
                <p className="text-sm text-text-secondary font-body">
                  Your estimation history will appear here after you submit your first request.
                </p>
              </div>
            )}
          </div>
        </div>
      </div>

      <style jsx>{`
        .animate-fadeIn {
          animation: fadeIn 0.5s ease-in-out;
        }
        
        @keyframes fadeIn {
          from { opacity: 0; transform: translateY(10px); }
          to { opacity: 1; transform: translateY(0); }
        }
        
        .custom-scrollbar::-webkit-scrollbar {
          width: 6px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-track {
          background: #172A45;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb {
          background: #5E81F4;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb:hover {
          background: #00F5D4;
        }
      `}</style>
    </div>
  );
};

export default CostEstimationForm;