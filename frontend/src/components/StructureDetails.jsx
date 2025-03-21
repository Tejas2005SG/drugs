import React, { useState, useEffect, useRef } from 'react';
import { jsPDF } from 'jspdf';
import { Canvg } from 'canvg'; // Import Canvg for SVG rendering
import axios from 'axios';
import EnhancedGeminiInfo from './Enhancedgeminiinfo.jsx'; // Import the new component

const StructureDetails = ({ structure, rdkitLoaded }) => {
  const [activeTab, setActiveTab] = useState('parent');
  const [selectedVariant, setSelectedVariant] = useState(null);
  const [variantInfo, setVariantInfo] = useState(null); // Store detailed info for the selected variant
  const [loadingVariantInfo, setLoadingVariantInfo] = useState(false); // Loading state for variant info
  const [variantInfoError, setVariantInfoError] = useState(null); // Error state for variant info
  const parentMolRef = useRef(null);
  const variantMolRefs = useRef({});

  // Reset state when structure changes
  useEffect(() => {
    setSelectedVariant(null);
    setActiveTab('parent');
    setVariantInfo(null); // Reset variant info when structure changes
    setVariantInfoError(null);
  }, [structure]);

  // Render molecule effect for parent structure
  useEffect(() => {
    if (rdkitLoaded && window.RDKit && structure?.smiles && activeTab === 'parent') {
      const mol = window.RDKit.get_mol(structure.smiles);
      if (mol && parentMolRef.current) {
        parentMolRef.current.innerHTML = mol.get_svg(600, 400);
        mol.delete();
      }
    }
  }, [rdkitLoaded, structure, activeTab]);

  // Render molecule effect for selected variant
  useEffect(() => {
    if (
      rdkitLoaded &&
      window.RDKit &&
      selectedVariant?.smiles &&
      activeTab === 'variant' &&
      parentMolRef.current
    ) {
      const mol = window.RDKit.get_mol(selectedVariant.smiles);
      if (mol) {
        parentMolRef.current.innerHTML = mol.get_svg(600, 400);
        mol.delete();
      }
    }
  }, [rdkitLoaded, selectedVariant, activeTab]);

  // Render molecules for variants tab
  useEffect(() => {
    if (rdkitLoaded && window.RDKit && activeTab === 'variants' && structure?.generatedStructures) {
      structure.generatedStructures.forEach((variant, index) => {
        const ref = variantMolRefs.current[index];
        if (ref && variant.smiles) {
          const mol = window.RDKit.get_mol(variant.smiles);
          if (mol) {
            ref.innerHTML = mol.get_svg(160, 120);
            mol.delete();
          }
        }
      });
    }
  }, [rdkitLoaded, structure, activeTab]);

  // Fetch detailed information for the selected variant from Gemini API
  const fetchVariantInfo = async (smiles) => {
    try {
      setLoadingVariantInfo(true);
      setVariantInfoError(null);

      const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;
      if (!GEMINI_API_KEY) {
        throw new Error('Gemini API key is not set in environment variables (VITE_GEMINI_API_KEY)');
      }

      const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";
      const prompt = `Provide detailed information about the molecule with the SMILES string "${smiles}". Include its chemical properties, potential uses, and any relevant biological or medicinal information.`;

      const response = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        {
          contents: [
            {
              parts: [
                {
                  text: prompt
                }
              ]
            }
          ]
        },
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );

      console.log('Gemini API Response for Variant:', JSON.stringify(response.data, null, 2));

      if (response.data.candidates && response.data.candidates.length > 0) {
        const info = response.data.candidates[0].content.parts[0].text;
        setVariantInfo(info);
      } else {
        throw new Error('No content returned from Gemini API');
      }
    } catch (error) {
      console.error('Error fetching variant information from Gemini API:', error);
      setVariantInfoError('Failed to retrieve detailed information about the variant: ' + error.message);
      setVariantInfo(null);
    } finally {
      setLoadingVariantInfo(false);
    }
  };

  const handleVariantSelect = (variant) => {
    setSelectedVariant(variant);
    setActiveTab('variant');
    setVariantInfo(null); // Reset variant info when selecting a new variant
    setVariantInfoError(null);
    fetchVariantInfo(variant.smiles); // Fetch detailed info for the selected variant
  };

  const renderMoleculeFallback = (smiles) => {
    if (!rdkitLoaded || !smiles || !window.RDKit) {
      return (
        <div className="border-2 border-dashed border-gray-300 rounded-lg p-8 flex items-center justify-center">
          <p className="text-gray-500">Molecule visualization unavailable</p>
        </div>
      );
    }
    return <div ref={parentMolRef} className="flex justify-center" />;
  };

  // Export to PDF function
  const exportToPDF = async () => {
    const doc = new jsPDF();
    const pageWidth = doc.internal.pageSize.getWidth();
    const margin = 10;
    let y = margin;

    // Helper function to add text and update y position
    const addText = (text, x, size = 12) => {
      doc.setFontSize(size);
      // Split text into lines if it's too long
      const lines = doc.splitTextToSize(text, pageWidth - 2 * margin);
      doc.text(lines, x, y);
      y += lines.length * size * 0.5; // Adjust line spacing based on number of lines
      return y;
    };

    // Add title
    addText(`${structure.name} Details`, margin, 16);
    y += 5;

    // Add creation date
    addText(`Created: ${new Date(structure.created).toLocaleString()}`, margin, 10);

    // Function to add SVG to PDF
    const addSvgToPDF = async (svgString, x, y, width, height) => {
      if (!svgString) {
        console.error('No SVG string provided for rendering');
        return y;
      }

      try {
        // Log the SVG string for debugging
        console.log('SVG String:', svgString);

        // Create a canvas to render the SVG
        const canvas = document.createElement('canvas');
        canvas.width = width * 2; // Higher resolution for better quality
        canvas.height = height * 2;
        const ctx = canvas.getContext('2d');

        // Use Canvg to render SVG to canvas
        const v = Canvg.fromString(ctx, svgString);
        await v.render();

        // Convert canvas to PNG data URL
        const imgData = canvas.toDataURL('image/png');
        if (!imgData || imgData === 'data:,') {
          console.error('Failed to generate image data from canvas');
          return y;
        }

        // Add image to PDF
        doc.addImage(imgData, 'PNG', x, y, width, height);
        return y + height + 10; // Update y position
      } catch (error) {
        console.error('Error rendering SVG to PDF:', error);
        return y;
      }
    };

    // Add active tab content
    if (activeTab === 'parent') {
      addText('Parent Structure', margin, 14);
      addText(`SMILES: ${structure.smiles}`, margin);
      addText('Generation Parameters:', margin, 12);
      addText(`Algorithm: ${structure.properties.algorithm || 'CMA-ES'}`, margin + 5);
      addText(`Property: ${structure.properties.propertyName || 'QED'}`, margin + 5);
      addText(`Optimization: ${structure.properties.minimize ? 'Minimize' : 'Maximize'}`, margin + 5);
      addText(`Min Similarity: ${structure.properties.minSimilarity || 0.3}`, margin + 5);

      // Add detailed information from Gemini API
      addText('Detailed Information:', margin, 12);
      const infoText = structure.information || 'No detailed information available.';
      y = addText(infoText, margin + 5, 10);

      // Add molecule SVG if RDKit is loaded
      if (rdkitLoaded && window.RDKit && structure.smiles) {
        try {
          const mol = window.RDKit.get_mol(structure.smiles);
          if (mol && mol.is_valid()) {
            const svg = mol.get_svg(300, 200); // Reduced size for better rendering
            y = await addSvgToPDF(svg, margin, y, 90, 60); // Match the size used for variants
            mol.delete();
          } else {
            console.error('Invalid molecule for SMILES:', structure.smiles);
            addText('Molecule visualization unavailable', margin, 10);
          }
        } catch (error) {
          console.error('Error generating SVG for parent structure:', error);
          addText('Molecule visualization unavailable', margin, 10);
        }
      } else {
        addText('Molecule visualization unavailable (RDKit not loaded)', margin, 10);
      }
    } else if (activeTab === 'variant' && selectedVariant) {
      addText('Selected Variant', margin, 14);
      addText(`SMILES: ${selectedVariant.smiles}`, margin);
      addText('Properties:', margin, 12);
      addText(`QED: ${selectedVariant.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5);
      addText(`LogP: ${selectedVariant.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5);
      addText(`Similarity to Parent: ${(selectedVariant.similarity * 100).toFixed(1)}%`, margin + 5);

      // Add detailed information for the variant
      addText('Detailed Information:', margin, 12);
      const variantInfoText = variantInfo || variantInfoError || 'Fetching detailed information...';
      y = addText(variantInfoText, margin + 5, 10);

      // Add molecule SVG if RDKit is loaded
      if (rdkitLoaded && window.RDKit && selectedVariant.smiles) {
        try {
          const mol = window.RDKit.get_mol(selectedVariant.smiles);
          if (mol && mol.is_valid()) {
            const svg = mol.get_svg(300, 200); // Reduced size for consistency
            y = await addSvgToPDF(svg, margin, y, 90, 60);
            mol.delete();
          } else {
            console.error('Invalid molecule for SMILES:', selectedVariant.smiles);
            addText('Molecule visualization unavailable', margin, 10);
          }
        } catch (error) {
          console.error('Error generating SVG for selected variant:', error);
          addText('Molecule visualization unavailable', margin, 10);
        }
      } else {
        addText('Molecule visualization unavailable (RDKit not loaded)', margin, 10);
      }
    } else if (activeTab === 'variants') {
      addText('Generated Variants', margin, 14);
      if (structure.generatedStructures) {
        for (let index = 0; index < structure.generatedStructures.length; index++) {
          const variant = structure.generatedStructures[index];
          addText(`Variant ${index + 1}:`, margin, 12);
          addText(`SMILES: ${variant.smiles}`, margin + 5);
          addText(`QED: ${variant.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5);
          addText(`LogP: ${variant.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5);
          addText(`Similarity: ${(variant.similarity * 100).toFixed(1)}%`, margin + 5);

          // Add SVG for each variant
          if (rdkitLoaded && window.RDKit && variant.smiles) {
            try {
              const mol = window.RDKit.get_mol(variant.smiles);
              if (mol && mol.is_valid()) {
                const svg = mol.get_svg(160, 120);
                y = await addSvgToPDF(svg, margin, y, 90, 60); // Same size as before
                mol.delete();
              } else {
                console.error('Invalid molecule for SMILES:', variant.smiles);
                addText('Molecule visualization unavailable', margin + 5, 10);
              }
            } catch (error) {
              console.error('Error generating SVG for variant:', error);
              addText('Molecule visualization unavailable', margin + 5, 10);
            }
          } else {
            addText('Molecule visualization unavailable (RDKit not loaded)', margin + 5, 10);
          }
          y += 5; // Extra spacing between variants
        }
      }
    }

    // Save the PDF
    doc.save(`${structure.name}-${activeTab}.pdf`);
  };

  return (
    <div>
      <div className="flex justify-between items-center mb-6">
        <h2 className="text-2xl font-bold">{structure.name}</h2>
        <span className="text-sm text-gray-500">
          Created: {new Date(structure.created).toLocaleString()}
        </span>
      </div>

      <div className="border-b border-gray-200 mb-6">
        <nav className="-mb-px flex space-x-8">
          <button
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'parent'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            }`}
            onClick={() => setActiveTab('parent')}
          >
            Parent Structure
          </button>
          <button
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'variants'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            }`}
            onClick={() => setActiveTab('variants')}
          >
            Generated Variants ({structure.generatedStructures?.length || 0})
          </button>
          {selectedVariant && (
            <button
              className={`pb-4 px-1 border-b-2 font-medium text-sm ${
                activeTab === 'variant'
                  ? 'border-blue-500 text-blue-600'
                  : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
              }`}
              onClick={() => setActiveTab('variant')}
            >
              Selected Variant
            </button>
          )}
        </nav>
      </div>

      {activeTab === 'parent' && (
        <div>
          <div className="bg-gray-50 p-4 rounded-lg mb-6">
            <h3 className="text-lg font-semibold mb-2">SMILES</h3>
            <p className="font-mono text-sm break-all">{structure.smiles}</p>
          </div>

          {renderMoleculeFallback(structure.smiles)}

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Detailed Information</h3>
            <EnhancedGeminiInfo 
              information={structure.information}
              isLoading={false}
              error={null}
            />
          </div>

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Generation Parameters</h3>
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">Algorithm</span>
                <p>{structure.properties.algorithm || 'CMA-ES'}</p>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">Property</span>
                <p>{structure.properties.propertyName || 'QED'}</p>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">Optimization</span>
                <p>{structure.properties.minimize ? 'Minimize' : 'Maximize'}</p>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">Min Similarity</span>
                <p>{structure.properties.minSimilarity || 0.3}</p>
              </div>
            </div>
          </div>
          <div className="mt-6 flex justify-center">
            <button
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded focus:outline-none focus:shadow-outline"
              onClick={exportToPDF}
            >
              Export to PDF
            </button>
          </div>
        </div>
      )}

      {activeTab === 'variants' && (
        <div>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            {structure.generatedStructures?.map((variant, index) => (
              <div
                key={index}
                className="border border-gray-200 rounded-lg p-4 hover:shadow-md cursor-pointer"
                onClick={() => handleVariantSelect(variant)}
              >
                <div className="h-40 flex items-center justify-center">
                  {rdkitLoaded && window.RDKit && variant.smiles ? (
                    <div
                      ref={(el) => (variantMolRefs.current[index] = el)}
                      className="w-full h-full"
                    />
                  ) : (
                    <div className="w-full h-full flex items-center justify-center bg-gray-100 rounded">
                      <p className="text-gray-500 text-xs">Visualization unavailable</p>
                    </div>
                  )}
                </div>
                <div className="mt-2">
                  <div className="flex justify-between text-sm">
                    <span className="font-medium">QED:</span>
                    <span>{variant.properties?.qed?.toFixed(3) || 'N/A'}</span>
                  </div>
                  <div className="flex justify-between text-sm">
                    <span className="font-medium">LogP:</span>
                    <span>{variant.properties?.logp?.toFixed(3) || 'N/A'}</span>
                  </div>
                  <div className="flex justify-between text-sm">
                    <span className="font-medium">Similarity:</span>
                    <span>{(variant.similarity * 100).toFixed(1)}%</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <button
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded focus:outline-none focus:shadow-outline"
              onClick={exportToPDF}
            >
              Export to PDF
            </button>
          </div>
        </div>
      )}

      {activeTab === 'variant' && selectedVariant && (
        <div>
          <div className="bg-gray-50 p-4 rounded-lg mb-6">
            <h3 className="text-lg font-semibold mb-2">Variant SMILES</h3>
            <p className="font-mono text-sm break-all">{selectedVariant.smiles}</p>
          </div>

          {renderMoleculeFallback(selectedVariant.smiles)}

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Properties</h3>
            <div className="grid grid-cols-3 gap-4">
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">QED</span>
                <p className="text-xl font-semibold">
                  {selectedVariant.properties?.qed?.toFixed(3) || 'N/A'}
                </p>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">LogP</span>
                <p className="text-xl font-semibold">
                  {selectedVariant.properties?.logp?.toFixed(3) || 'N/A'}
                </p>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <span className="text-sm font-medium text-gray-500">Similarity to Parent</span>
                <p className="text-xl font-semibold">
                  {(selectedVariant.similarity * 100).toFixed(1)}%
                </p>
              </div>
            </div>
          </div>

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Detailed Information</h3>
            <EnhancedGeminiInfo 
              information={variantInfo}
              isLoading={loadingVariantInfo}
              error={variantInfoError}
            />
          </div>

          <div className="mt-6 flex justify-center">
            <button
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded focus:outline-none focus:shadow-outline"
              onClick={exportToPDF}
            >
              Export to PDF
            </button>
          </div>
        </div>
      )}  
    </div>
  );
};

export default StructureDetails;