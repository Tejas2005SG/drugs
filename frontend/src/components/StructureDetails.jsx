import React, { useState, useEffect, useRef } from 'react';
import { jsPDF } from 'jspdf';
import { Canvg } from 'canvg';
import axios from 'axios';
import EnhancedGeminiInfo from './Enhancedgeminiinfo.jsx';

const StructureDetails = ({ structure, rdkitLoaded }) => {
  const [activeTab, setActiveTab] = useState('parent');
  const [selectedVariant, setSelectedVariant] = useState(null);
  const [variantInfo, setVariantInfo] = useState(null);
  const [loadingVariantInfo, setLoadingVariantInfo] = useState(false);
  const [variantInfoError, setVariantInfoError] = useState(null);
  const [showParent3D, setShowParent3D] = useState(false);
  const [showVariant3D, setShowVariant3D] = useState(false);
  const [is3DmolLoaded, setIs3DmolLoaded] = useState(false);
  const parentMolRef = useRef(null);
  const variantMolRefs = useRef({});
  const parent3DRef = useRef(null);
  const variant3DRef = useRef(null);

  useEffect(() => {
    if (showParent3D || showVariant3D) {
      if (!window.$3Dmol && !is3DmolLoaded) {
        const script = document.createElement('script');
        script.src = 'https://3dmol.csb.pitt.edu/build/3Dmol-min.js';
        script.async = true;
        script.onload = () => {
          console.log('3Dmol.js loaded successfully');
          setIs3DmolLoaded(true);
        };
        script.onerror = () => {
          console.error('Failed to load 3Dmol.js');
          setIs3DmolLoaded(false);
        };
        document.head.appendChild(script);
      } else if (window.$3Dmol) {
        setIs3DmolLoaded(true);
      }
    }
  }, [showParent3D, showVariant3D]);

  useEffect(() => {
    setSelectedVariant(null);
    setActiveTab('parent');
    setVariantInfo(null);
    setVariantInfoError(null);
    setShowParent3D(false);
    setShowVariant3D(false);
  }, [structure]);

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && structure?.smiles && activeTab === 'parent') {
      const mol = window.RDKit.get_mol(structure.smiles);
      if (mol && parentMolRef.current) {
        parentMolRef.current.innerHTML = mol.get_svg(600, 400);
        mol.delete();
      }
    }
  }, [rdkitLoaded, structure, activeTab]);

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && selectedVariant?.smiles && activeTab === 'variant') {
      const mol = window.RDKit.get_mol(selectedVariant.smiles);
      if (mol && parentMolRef.current) {
        parentMolRef.current.innerHTML = mol.get_svg(600, 400);
        mol.delete();
      }
    }
  }, [rdkitLoaded, selectedVariant, activeTab]);

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

  useEffect(() => {
    if (showParent3D && parent3DRef.current && structure?.smiles && is3DmolLoaded && window.$3Dmol) {
      parent3DRef.current.innerHTML = ''; // Clear container
      const mol = window.RDKit.get_mol(structure.smiles);
      const molData = mol.get_molblock();
      mol.delete();
      
      const viewer = window.$3Dmol.createViewer(parent3DRef.current, {
        backgroundColor: 'white'
      });
      if (!viewer) {
        console.error('Failed to create parent 3D viewer');
        return;
      }
      viewer.addModel(molData, 'mol');
      viewer.setStyle({}, { stick: { radius: 0.2 }, sphere: { scale: 0.3 } });
      viewer.zoomTo();
      viewer.render();
    }
  }, [showParent3D, structure, is3DmolLoaded]);

  useEffect(() => {
    if (showVariant3D && variant3DRef.current && selectedVariant?.smiles && is3DmolLoaded && window.$3Dmol) {
      variant3DRef.current.innerHTML = ''; // Clear container
      const mol = window.RDKit.get_mol(selectedVariant.smiles);
      const molData = mol.get_molblock();
      mol.delete();
      
      const viewer = window.$3Dmol.createViewer(variant3DRef.current, {
        backgroundColor: 'white'
      });
      if (!viewer) {
        console.error('Failed to create variant 3D viewer');
        return;
      }
      viewer.addModel(molData, 'mol');
      viewer.setStyle({}, { stick: { radius: 0.2 }, sphere: { scale: 0.3 } });
      viewer.zoomTo();
      viewer.render();
    }
  }, [showVariant3D, selectedVariant, is3DmolLoaded]);

  const fetchVariantInfo = async (smiles) => {
    try {
      setLoadingVariantInfo(true);
      setVariantInfoError(null);
      const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;
      if (!GEMINI_API_KEY) throw new Error('Gemini API key not set');
      const response = await axios.post(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: `Provide detailed information about the molecule with SMILES "${smiles}".` }] }] },
        { headers: { 'Content-Type': 'application/json' } }
      );
      if (response.data.candidates?.[0]?.content?.parts?.[0]?.text) {
        setVariantInfo(response.data.candidates[0].content.parts[0].text);
      } else {
        throw new Error('No content from Gemini API');
      }
    } catch (error) {
      setVariantInfoError('Failed to fetch variant info: ' + error.message);
      setVariantInfo(null);
    } finally {
      setLoadingVariantInfo(false);
    }
  };

  const handleVariantSelect = (variant) => {
    setSelectedVariant(variant);
    setActiveTab('variant');
    setVariantInfo(null);
    setVariantInfoError(null);
    fetchVariantInfo(variant.smiles);
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

  const exportToPDF = async () => {
    const doc = new jsPDF();
    const pageWidth = doc.internal.pageSize.getWidth();
    const margin = 10;
    let y = margin;

    const addText = (text, x, size = 12) => {
      doc.setFontSize(size);
      const lines = doc.splitTextToSize(text, pageWidth - 2 * margin);
      doc.text(lines, x, y);
      y += lines.length * size * 0.5;
      return y;
    };

    addText(`${structure.name} Details`, margin, 16);
    y += 5;
    addText(`Created: ${new Date(structure.created).toLocaleString()}`, margin, 10);

    const addSvgToPDF = async (svgString, x, y, width, height) => {
      if (!svgString) return y;
      const canvas = document.createElement('canvas');
      canvas.width = width * 2;
      canvas.height = height * 2;
      const ctx = canvas.getContext('2d');
      const v = Canvg.fromString(ctx, svgString);
      await v.render();
      const imgData = canvas.toDataURL('image/png');
      doc.addImage(imgData, 'PNG', x, y, width, height);
      return y + height + 10;
    };

    if (activeTab === 'parent') {
      addText('Parent Structure', margin, 14);
      addText(`SMILES: ${structure.smiles}`, margin);
      if (rdkitLoaded && window.RDKit && structure.smiles) {
        const mol = window.RDKit.get_mol(structure.smiles);
        if (mol) {
          const svg = mol.get_svg(300, 200);
          y = await addSvgToPDF(svg, margin, y, 90, 60);
          mol.delete();
        }
      }
      addText('Detailed Information:', margin, 12);
      y = addText(structure.information || 'No info available', margin + 5, 10);
    } else if (activeTab === 'variant' && selectedVariant) {
      addText('Selected Variant', margin, 14);
      addText(`SMILES: ${selectedVariant.smiles}`, margin);
      if (rdkitLoaded && window.RDKit && selectedVariant.smiles) {
        const mol = window.RDKit.get_mol(selectedVariant.smiles);
        if (mol) {
          const svg = mol.get_svg(300, 200);
          y = await addSvgToPDF(svg, margin, y, 90, 60);
          mol.delete();
        }
      }
      addText('Detailed Information:', margin, 12);
      y = addText(variantInfo || variantInfoError || 'Fetching info...', margin + 5, 10);
    } else if (activeTab === 'variants') {
      addText('Generated Variants', margin, 14);
      if (structure.generatedStructures) {
        for (let i = 0; i < structure.generatedStructures.length; i++) {
          const variant = structure.generatedStructures[i];
          addText(`Variant ${i + 1}:`, margin, 12);
          addText(`SMILES: ${variant.smiles}`, margin + 5);
          if (rdkitLoaded && window.RDKit && variant.smiles) {
            const mol = window.RDKit.get_mol(variant.smiles);
            if (mol) {
              const svg = mol.get_svg(160, 120);
              y = await addSvgToPDF(svg, margin, y, 90, 60);
              mol.delete();
            }
          }
          y += 5;
        }
      }
    }

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

          <div className="mt-4 flex justify-center">
            <button
              onClick={() => setShowParent3D(!showParent3D)}
              className="bg-blue-500 hover:bg-blue-600 text-white font-bold py-2 px-4 rounded"
              disabled={!rdkitLoaded}
            >
              {showParent3D ? 'Hide 3D View' : 'Show 3D Visualization'}
            </button>
          </div>

          {showParent3D && (
            <div 
              ref={parent3DRef}
              className="w-full h-96 mt-4 border rounded relative"
            >
              {!is3DmolLoaded && <p className="text-center pt-20">Loading 3D viewer...</p>}
            </div>
          )}

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Detailed Information</h3>
            <EnhancedGeminiInfo 
              information={structure.information}
              isLoading={false}
              error={null}
            />
          </div>

          <div className="mt-6 flex justify-center">
            <button
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded"
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
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded"
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

          <div className="mt-4 flex justify-center">
            <button
              onClick={() => setShowVariant3D(!showVariant3D)}
              className="bg-blue-500 hover:bg-blue-600 text-white font-bold py-2 px-4 rounded"
              disabled={!rdkitLoaded}
            >
              {showVariant3D ? 'Hide 3D View' : 'Show 3D Visualization'}
            </button>
          </div>

          {showVariant3D && (
            <div 
              ref={variant3DRef}
              className="w-full h-96 mt-4 border rounded relative"
            >
              {!is3DmolLoaded && <p className="text-center pt-20">Loading 3D viewer...</p>}
            </div>
          )}

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
              className="bg-green-500 hover:bg-green-600 text-white font-bold py-2 px-4 rounded"
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