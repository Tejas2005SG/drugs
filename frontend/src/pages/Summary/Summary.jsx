import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { jsPDF } from 'jspdf';
import { toast } from 'react-hot-toast';
import { useAuthStore } from '../../Store/auth.store.js';
import { Download, FileDown, Loader2, AlertCircle } from 'lucide-react';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode === "development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

const Summary = () => {
  const [summaryData, setSummaryData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const { user } = useAuthStore();

  useEffect(() => {
    fetchSummaryData();
  }, []);

  const fetchSummaryData = async () => {
    try {
      const response = await axiosInstance.get('/summary/summary');
      setSummaryData(response.data.data);
    } catch (err) {
      setError(err.response?.data?.message || 'Failed to fetch summary data');
      toast.error('Failed to fetch summary data');
    } finally {
      setLoading(false);
    }
  };

  const generatePDF = () => {
    if (!summaryData) {
      toast.error('No data available to generate PDF');
      return;
    }

    const doc = new jsPDF();
    let y = 20;
    const margin = 20;
    const pageHeight = doc.internal.pageSize.height;

    const addText = (text, size = 12, isBold = false) => {
      doc.setFontSize(size);
      doc.setFont('helvetica', isBold ? 'bold' : 'normal');
      const lines = doc.splitTextToSize(text, doc.internal.pageSize.width - 2 * margin);
      
      if (y + (lines.length * size * 0.5) > pageHeight - margin) {
        doc.addPage();
        y = margin;
      }
      
      doc.text(lines, margin, y);
      y += lines.length * size * 0.5 + 5;
    };

    // Title
    addText('Drug Discovery Summary Report', 24, true);
    addText(`Generated on: ${new Date().toLocaleString()}`, 10);
    addText(`Researcher: ${user?.firstName} ${user?.lastName}`, 10);
    y += 10;

    // New Drug Section
    if (summaryData.newDrugs?.length > 0) {
      addText('New Drug Discovery', 16, true);
      summaryData.newDrugs.forEach(drug => {
        addText(`Title: ${drug.newmoleculetitle}`, 12, true);
        addText(`SMILES: ${drug.newSmiles}`, 10);
        addText(`IUPAC Name: ${drug.newIupacName || 'Not available'}`, 10);
        if (drug.conversionDetails) addText(`Conversion Details: ${drug.conversionDetails}`, 10);
        y += 5;
      });
    }

    // Cost Estimation Section
    if (summaryData.costEstimations?.length > 0) {
      addText('Cost Estimation', 16, true);
      summaryData.costEstimations.forEach(cost => {
        addText(`SMILES: ${cost.smiles}`, 10);
        addText(`Estimated Cost: ${cost.estimatedcost}`, 10);
        y += 5;
      });
    }

    // AI Naming Section
    if (summaryData.drugNames?.length > 0) {
      addText('AI Generated Names', 16, true);
      summaryData.drugNames.forEach(name => {
        addText(`Suggested Name: ${name.suggestedName}`, 12, true);
        addText(`Status: ${name.status}`, 10);
        if (name.namingDetails) addText(`Details: ${name.namingDetails}`, 10);
        y += 5;
      });
    }

    // Research Papers Section
    if (summaryData.researchPapers?.length > 0) {
      addText('Generated Research Papers', 16, true);
      summaryData.researchPapers.forEach(paper => {
        addText(`Title: ${paper.title}`, 12, true);
        if (paper.abstract) addText(`Abstract: ${paper.abstract}`, 10);
        y += 5;
      });
    }

    // Target Prediction Section
    if (summaryData.targetPredictions?.length > 0) {
      addText('AI-Driven Target Predictions', 16, true);
      summaryData.targetPredictions.forEach(prediction => {
        addText(`SMILES: ${prediction.smiles}`, 10);
        if (prediction.targets) {
          prediction.targets.forEach(target => {
            addText(`Target: ${target.protein} (Confidence: ${target.confidence})`, 10);
          });
        }
        y += 5;
      });
    }

    // Toxicity Section
    if (summaryData.toxicityResults?.length > 0) {
      addText('Toxicity Predictions', 16, true);
      summaryData.toxicityResults.forEach(toxicity => {
        addText(`SMILES: ${toxicity.smiles}`, 10);
        if (toxicity.toxicityResult) {
          addText(`LD50: ${toxicity.toxicityResult.acuteToxicity.LD50}`, 10);
          addText(`Class: ${toxicity.toxicityResult.acuteToxicity.toxicityClass}`, 10);
        }
        y += 5;
      });
    }

    // Save the PDF
    doc.save('drug-discovery-summary.pdf');
    toast.success('PDF generated successfully!');
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <Loader2 className="h-8 w-8 animate-spin text-blue-600" />
        <span className="ml-2">Loading summary data...</span>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <AlertCircle className="h-8 w-8 text-red-600" />
        <span className="ml-2 text-red-600">{error}</span>
      </div>
    );
  }

  return (
    <div className="max-w-7xl mx-auto px-4 py-8">
      <div className="bg-white rounded-xl shadow-lg p-6">
        <div className="flex justify-between items-center mb-6">
          <h1 className="text-2xl font-bold text-gray-800">Drug Discovery Summary</h1>
          <button
            onClick={generatePDF}
            className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
          >
            <Download className="h-5 w-5 mr-2" />
            Download PDF
          </button>
        </div>

        {/* Missing Data Warning */}
        {(!summaryData?.newDrugs?.length || 
          !summaryData?.costEstimations?.length || 
          !summaryData?.drugNames?.length || 
          !summaryData?.researchPapers?.length || 
          !summaryData?.targetPredictions?.length || 
          !summaryData?.toxicityResults?.length) && (
          <div className="bg-yellow-50 border-l-4 border-yellow-400 p-4 mb-6">
            <div className="flex">
              <AlertCircle className="h-5 w-5 text-yellow-400" />
              <div className="ml-3">
                <p className="text-sm text-yellow-700">
                  Some data is missing from your summary. Complete the following sections for a comprehensive report:
                </p>
                <ul className="mt-2 list-disc list-inside text-sm text-yellow-600">
                  {!summaryData?.newDrugs?.length && <li>New Drug Discovery</li>}
                  {!summaryData?.costEstimations?.length && <li>Cost Estimation</li>}
                  {!summaryData?.drugNames?.length && <li>AI Naming Suggestions</li>}
                  {!summaryData?.researchPapers?.length && <li>Research Papers</li>}
                  {!summaryData?.targetPredictions?.length && <li>Target Predictions</li>}
                  {!summaryData?.toxicityResults?.length && <li>Toxicity Predictions</li>}
                </ul>
              </div>
            </div>
          </div>
        )}

        {/* Summary Sections */}
        <div className="space-y-8">
          {/* New Drugs Section */}
          <section>
            <h2 className="text-xl font-semibold text-gray-800 mb-4">New Drug Discovery</h2>
            {summaryData?.newDrugs?.length > 0 ? (
              <div className="space-y-4">
                {summaryData.newDrugs.map((drug, index) => (
                  <div key={index} className="bg-gray-50 p-4 rounded-lg">
                    <h3 className="font-medium text-gray-900">{drug.newmoleculetitle}</h3>
                    <p className="text-sm text-gray-600 mt-1">SMILES: {drug.newSmiles}</p>
                    {drug.newIupacName && (
                      <p className="text-sm text-gray-600">IUPAC: {drug.newIupacName}</p>
                    )}
                  </div>
                ))}
              </div>
            ) : (
              <p className="text-gray-500">No new drugs discovered yet</p>
            )}
          </section>

          {/* Cost Estimation Section */}
          <section>
            <h2 className="text-xl font-semibold text-gray-800 mb-4">Cost Estimation</h2>
            {summaryData?.costEstimations?.length > 0 ? (
              <div className="space-y-4">
                {summaryData.costEstimations.map((cost, index) => (
                  <div key={index} className="bg-gray-50 p-4 rounded-lg">
                    <p className="text-sm text-gray-600">SMILES: {cost.smiles}</p>
                    <p className="font-medium text-gray-900 mt-1">
                      Estimated Cost: {cost.estimatedcost}
                    </p>
                  </div>
                ))}
              </div>
            ) : (
              <p className="text-gray-500">No cost estimations available</p>
            )}
          </section>

          {/* Continue with similar sections for other data types */}
          {/* Add sections for AI Naming, Research Papers, Target Predictions, and Toxicity */}
        </div>
      </div>
    </div>
  );
};

export default Summary;