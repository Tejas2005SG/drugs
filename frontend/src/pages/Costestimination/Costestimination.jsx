// src/components/CostEstimationForm.jsx
"use client"

import { useState, useEffect, useRef } from "react"
import { useAuthStore } from "../../Store/auth.store.js"
import { postCostEstimation, getCostEstimations } from "../../api/costestimination.jsx"
import { AlertCircle, DollarSign, Clock, Database, X, Info, RefreshCw, LogIn, FileText, ChevronDown, ChevronUp, Download } from "lucide-react"
import { jsPDF } from "jspdf"
import domtoimage from "dom-to-image"
const CostEstimationForm = () => {
  const [smiles, setSmiles] = useState("")
  const [result, setResult] = useState(null)
  const [history, setHistory] = useState([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [historyLoading, setHistoryLoading] = useState(false)
  const [isResultInfoOpen, setIsResultInfoOpen] = useState(false) // State for result dropdown
  const [openHistoryItems, setOpenHistoryItems] = useState({}) // State for history dropdowns

  const { user, checkAuth, checkingAuth } = useAuthStore()

  // Ref to capture the information section for PDF export
  const infoRef = useRef(null)

  const handleSubmit = async (e) => {
    e.preventDefault()
    if (!user?._id) {
      setError("Please log in to estimate costs")
      return
    }
    setLoading(true)
    setError(null)
    try {
      const data = await postCostEstimation(smiles)
      setResult(data.data)
      await fetchHistory()
    } catch (error) {
      setError("Failed to estimate cost. Please try again.")
      console.error("Error estimating cost:", error)
    } finally {
      setLoading(false)
    }
  }

  const fetchHistory = async () => {
    if (!user?._id) {
      setError("Please log in to view history")
      return
    }
    setHistoryLoading(true)
    try {
      const data = await getCostEstimations(user._id)
      const validHistory = Array.isArray(data.data) ? data.data.filter((item) => item && typeof item === "object") : []
      setHistory(validHistory)
    } catch (error) {
      setError("Failed to fetch history.")
      console.error("Error fetching history:", error)
    } finally {
      setHistoryLoading(false)
    }
  }

  const toggleResultInfo = () => {
    setIsResultInfoOpen(!isResultInfoOpen)
  }

  const toggleHistoryItem = (id) => {
    setOpenHistoryItems((prev) => ({
      ...prev,
      [id]: !prev[id],
    }))
  }

  // Function to render the information content
  const renderInformation = (info) => {
    if (!info) return <p style={{ color: '#374151' }}>No additional information available</p>;

    const lines = info.split('\n').filter(line => line.trim() !== '');
    let introParagraph = '';
    const sections = [];
    let currentSection = null;

    lines.forEach((line, index) => {
      if (index === 0 && !line.match(/^\d+\./)) {
        introParagraph = line.trim();
      } else if (line.match(/^\d+\.\s/)) {
        if (currentSection) {
          sections.push(currentSection);
        }
        currentSection = { title: line.trim(), bullets: [] };
      } else if (line.trim().startsWith('-')) {
        if (currentSection) {
          currentSection.bullets.push(line.trim().replace('-', '').trim());
        }
      } else {
        if (currentSection) {
          currentSection.bullets.push(line.trim());
        }
      }
    });

    if (currentSection) {
      sections.push(currentSection);
    }

    return (
      <div>
        {introParagraph && <p className="mb-2" style={{ color: '#374151' }}>{introParagraph}</p>}
        {sections.length > 0 && (
          <ul className="list-disc pl-5 space-y-2">
            {sections.map((section, index) => (
              <li key={index} className="font-semibold" style={{ color: '#1f2937' }}>
                {section.title}
                {section.bullets.length > 0 && (
                  <ul className="list-disc pl-5 mt-1 space-y-1">
                    {section.bullets.map((bullet, bulletIndex) => (
                      <li key={bulletIndex} style={{ color: '#374151' }}>
                        {bullet}
                      </li>
                    ))}
                  </ul>
                )}
              </li>
            ))}
          </ul>
        )}
      </div>
    );
  }

  // Function to export the SMILES, estimated cost, and information as a PDF
  const exportToPDF = async () => {
    const pdf = new jsPDF('p', 'mm', 'a4');
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 10;
    let yPosition = margin;

    // Add title
    pdf.setFontSize(16);
    pdf.text("Drug Cost Estimation Report", margin, yPosition);
    yPosition += 10;

    // Add SMILES
    pdf.setFontSize(12);
    pdf.text("SMILES:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result.smiles || "N/A", margin, yPosition);
    yPosition += 8;

    // Add Estimated Cost
    pdf.setFontSize(12);
    pdf.text("Estimated Cost:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result.estimatedcost || "N/A", margin, yPosition);
    yPosition += 10;

    // Add Information section
    if (infoRef.current) {
      pdf.setFontSize(12);
      pdf.text("Information:", margin, yPosition);
      yPosition += 8;

      // Use dom-to-image to capture the information section
      const infoImgData = await domtoimage.toPng(infoRef.current, { quality: 1 });
      const infoImgProps = pdf.getImageProperties(infoImgData);
      const infoImgWidth = pageWidth - 2 * margin;
      let infoImgHeight = (infoImgProps.height * infoImgWidth) / infoImgProps.width;

      // Handle pagination for the information section
      let remainingHeight = infoImgHeight;
      let yOffset = 0;

      while (remainingHeight > 0) {
        const spaceLeft = pageHeight - yPosition - margin;
        const heightToRender = Math.min(remainingHeight, spaceLeft);

        // Crop the image to fit the current page
        const tempCanvas = document.createElement('canvas');
        const tempCtx = tempCanvas.getContext('2d');
        const img = new Image();
        img.src = infoImgData;
        await new Promise((resolve) => { img.onload = resolve; });
        tempCanvas.width = img.width;
        tempCanvas.height = (heightToRender / infoImgHeight) * img.height;
        tempCtx.drawImage(img, 0, yOffset, img.width, tempCanvas.height, 0, 0, img.width, tempCanvas.height);
        const croppedImgData = tempCanvas.toDataURL('image/png');

        pdf.addImage(croppedImgData, 'PNG', margin, yPosition, infoImgWidth, heightToRender);
        yPosition += heightToRender;
        remainingHeight -= heightToRender;
        yOffset += tempCanvas.height;

        if (remainingHeight > 0) {
          pdf.addPage();
          yPosition = margin;
        }
      }
    }

    // Save the PDF
    pdf.save("drug-cost-estimation-report.pdf");
  };

  useEffect(() => {
    const initialize = async () => {
      try {
        console.log("Starting auth check...")
        await checkAuth()
        console.log("After checkAuth - User:", user, "CheckingAuth:", checkingAuth)
        if (user?._id) {
          console.log("Fetching history for user:", user._id)
          await fetchHistory()
        } else {
          console.log("No user ID found after checkAuth")
        }
      } catch (err) {
        setError("Failed to verify authentication. Please try refreshing the page.")
        console.error("Error during initialization:", err)
      }
    }
    initialize()
  }, [checkAuth])

  if (checkingAuth) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-gray-50">
        <div className="animate-spin mb-4">
          <RefreshCw size={32} className="text-blue-600" />
        </div>
        <p className="text-gray-700 font-medium">Verifying authentication...</p>
      </div>
    )
  }

  if (!user || !user._id) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-gray-50">
        <div className="text-center max-w-md p-8 bg-white rounded-xl shadow-lg">
          <LogIn size={48} className="mx-auto text-blue-600 mb-4" />
          <h2 className="text-2xl font-bold text-gray-800 mb-3">Authentication Required</h2>
          <p className="text-gray-600 mb-6">
            Please log in to access the Drug Cost Estimator tool and view your estimation history.
          </p>
          <button
            className="w-full px-4 py-3 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors duration-200 flex items-center justify-center"
            onClick={() => (window.location.href = "/login")}
          >
            <LogIn size={18} className="mr-2" />
            Go to Login
          </button>
        </div>
      </div>
    )
  }

  return (
    <div className="min-h-screen bg-gray-50 py-12">
      <div className="container mx-auto px-4 max-w-4xl">
        <div className="mb-10 text-center">
          <h1 className="text-4xl font-bold text-gray-800 mb-2">Drug Cost Estimator</h1>
          <p className="text-gray-600 max-w-2xl mx-auto">
            Enter a SMILES string to estimate the cost of drug synthesis and production.
          </p>
        </div>

        <div className="space-y-8">
          {/* Form Section */}
          <div className="bg-white p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg">
            <div className="flex items-center mb-6">
              <Database className="h-6 w-6 text-blue-600 mr-2" />
              <h2 className="text-2xl font-semibold text-gray-800">Estimate New Cost</h2>
            </div>

            <form onSubmit={handleSubmit} className="space-y-4">
              <div>
                <label htmlFor="smiles" className="block text-sm font-medium text-gray-700 mb-1">
                  SMILES String
                </label>
                <div className="relative">
                  <input
                    id="smiles"
                    type="text"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="Enter SMILES string (e.g., CC(=O)OC1=CC=CC=C1C(=O)O)"
                    required
                    disabled={loading}
                    className="w-full px-4 py-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent transition-all duration-200"
                  />
                  <FileText className="absolute right-3 top-3 text-gray-400" size={20} />
                </div>
                <p className="mt-1 text-xs text-gray-500">Enter a valid SMILES notation for the chemical compound</p>
              </div>

              <button
                type="submit"
                disabled={loading}
                className={`w-full px-4 py-3 text-white rounded-lg transition-all duration-200 flex items-center justify-center ${loading ? "bg-gray-400 cursor-not-allowed" : "bg-blue-600 hover:bg-blue-700"
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

            {result && (
              <div className="mt-8 animate-fadeIn">
                <div className="flex items-center mb-4">
                  <DollarSign className="h-5 w-5 text-green-600 mr-2" />
                  <h3 className="text-xl font-semibold text-gray-800">Estimation Result</h3>
                </div>
                <div className="bg-gradient-to-r from-blue-50 to-green-50 p-6 rounded-lg border border-blue-100">
                  <div className="space-y-3">
                    <div>
                      <p className="text-sm text-gray-500 mb-1">SMILES</p>
                      <p className="font-mono text-sm bg-white p-2 rounded border border-gray-200 overflow-x-auto">
                        {result.smiles || "N/A"}
                      </p>
                    </div>

                    <div>
                      <p className="text-sm text-gray-500 mb-1">Estimated Cost</p>
                      <p className="text-3xl font-bold text-green-600">
                        {result.estimatedcost || "N/A"}
                      </p>
                    </div>

                    <div className="space-y-2">
                      <button
                        onClick={toggleResultInfo}
                        className="flex items-center justify-between w-full p-3 bg-amber-50 hover:bg-amber-100 rounded-lg transition-colors duration-200 border border-amber-200 focus:outline-none focus:ring-2 focus:ring-amber-300 focus:ring-offset-2"
                      >
                        <div className="flex items-center space-x-2">
                          <Info className="h-5 w-5 text-amber-600" />
                          <span className="text-base font-semibold text-amber-800">Detailed Analysis <span className="text-sm text-blue-700">(powered by Gemini)</span></span>
                        </div>
                        {isResultInfoOpen ? (
                          <ChevronUp className="h-5 w-5 text-amber-600 transform transition-transform duration-300" />
                        ) : (
                          <ChevronDown className="h-5 w-5 text-amber-600 transform transition-transform duration-300" />
                        )}
                      </button>

                      {isResultInfoOpen && (
                        <div className="space-y-4">
                          <div
                            ref={infoRef}
                            className="p-4 bg-white rounded-lg shadow-lg border border-gray-100 transition-all duration-300"
                          >
                            <div className="prose prose-sm text-gray-600 max-w-none">
                              {renderInformation(result.information)}
                            </div>
                          </div>

                          <div className="flex justify-end">
                            <button
                              onClick={exportToPDF}
                              className="flex items-center px-4 py-2.5 bg-gradient-to-br from-green-500 to-green-600 text-white font-medium rounded-md hover:from-green-600 hover:to-green-700 transition-all duration-200 transform hover:scale-[1.02] shadow-sm hover:shadow-md focus:outline-none focus:ring-2 focus:ring-green-500 focus:ring-offset-2"
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
          <div className="bg-white p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg">
            <div className="flex items-center justify-between mb-6">
              <div className="flex items-center">
                <Clock className="h-6 w-6 text-blue-600 mr-2" />
                <h2 className="text-2xl font-semibold text-gray-800">Estimation History</h2>
              </div>
              <button
                onClick={fetchHistory}
                disabled={historyLoading}
                className="text-blue-600 hover:text-blue-800 transition-colors duration-200"
                title="Refresh history"
              >
                <RefreshCw size={20} className={historyLoading ? "animate-spin" : ""} />
              </button>
            </div>

            {historyLoading ? (
              <div className="flex flex-col items-center justify-center py-12">
                <RefreshCw size={24} className="animate-spin mb-4 text-blue-600" />
                <p className="text-gray-600">Loading history...</p>
              </div>
            ) : history.length > 0 ? (
              <div className="space-y-4 max-h-[500px] overflow-y-auto pr-2 custom-scrollbar">
                {history.map((item) => (
                  <div
                    key={item._id}
                    className="border border-gray-200 p-4 rounded-lg hover:border-blue-300 transition-all duration-200 hover:shadow-sm"
                  >
                    <div className="flex justify-between items-start mb-2">
                      <div className="flex items-center">
                        <DollarSign size={16} className="text-green-600 mr-1" />
                        <span className="font-bold text-lg text-green-600">
                          {item.estimatedcost || "N/A"}
                        </span>
                      </div>
                      <span className="text-xs text-gray-500 bg-gray-100 px-2 py-1 rounded-full">
                        {item.created ? new Date(item.created).toLocaleDateString() : "N/A"}
                      </span>
                    </div>

                    <div className="mt-2">
                      <p className="text-xs text-gray-500 mb-1">SMILES</p>
                      <p className="font-mono text-sm bg-gray-50 p-2 rounded border border-gray-200 overflow-x-auto">
                        {item.smiles || "N/A"}
                      </p>
                    </div>

                    {item.information && (
                      <div className="mt-3 space-y-2">
                        <button
                          onClick={() => toggleHistoryItem(item._id)}
                          className="flex items-center justify-between w-full p-2.5 bg-amber-50 hover:bg-amber-100 rounded-md transition-colors duration-200 border border-amber-200 focus:outline-none focus:ring-2 focus:ring-amber-300 focus:ring-offset-1"
                        >
                          <div className="flex items-center space-x-2">
                            <Info className="h-4 w-4 text-amber-600" />
                            <span className="text-xs font-medium text-amber-800 uppercase tracking-wide">Details <span className="text-xs text-blue-700">(powered by Gemini)</span></span>
                          </div>
                          {openHistoryItems[item._id] ? (
                            <ChevronUp className="h-4 w-4 text-amber-600 transform transition-transform duration-300" />
                          ) : (
                            <ChevronDown className="h-4 w-4 text-amber-600 transform transition-transform duration-300" />
                          )}
                        </button>

                        {openHistoryItems[item._id] && (
                          <div className="space-y-3">
                            <div className="p-3 bg-white rounded-md shadow-sm border border-gray-200 animate-fade-in">
                              <div className="prose prose-sm text-gray-600 max-w-none">
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
                <Database size={32} className="text-gray-400 mb-4" />
                <p className="text-gray-500 mb-2">No previous estimations found.</p>
                <p className="text-sm text-gray-400">
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
          background: #f1f1f1;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb {
          background: #d1d5db;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb:hover {
          background: #9ca3af;
        }
      `}</style>
    </div>
  )
}

export default CostEstimationForm