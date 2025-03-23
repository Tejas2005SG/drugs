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
  const [parentStyle, setParentStyle] = useState('ballstick');
  const [variantStyle, setVariantStyle] = useState('ballstick');
  const [showLabels, setShowLabels] = useState(false);
  const [aiDescription, setAiDescription] = useState('');
  const [spinAnimation, setSpinAnimation] = useState(false);
  const [highlightedAtom, setHighlightedAtom] = useState(null);
  const parentMolRef = useRef(null);
  const variantMolRefs = useRef({});
  const parent3DRef = useRef(null);
  const variant3DRef = useRef(null);
  const parentViewerRef = useRef(null);
  const variantViewerRef = useRef(null);
  const animationFrameId = useRef(null); // For smooth spinning

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
    setVariantInfo(null);
    setVariantInfoError(null);
    setShowParent3D(false);
    setShowVariant3D(false);
    setAiDescription('');
    setSpinAnimation(false);
    setHighlightedAtom(null);
    if (parent3DRef.current) parent3DRef.current.innerHTML = '';
    if (variant3DRef.current) variant3DRef.current.innerHTML = '';
    if (animationFrameId.current) {
      cancelAnimationFrame(animationFrameId.current); // Clean up animation
      animationFrameId.current = null;
    }
  }, [structure, activeTab]);

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

  const fetchAiEnhancements = async (smiles) => {
    const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;
    if (!GEMINI_API_KEY) {
      setAiDescription('Gemini API key not set');
      return;
    }

    try {
      const response = await axios.post(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: `Generate a creative, concise description of potential applications for a molecule with SMILES "${smiles}".` }] }] },
        { headers: { 'Content-Type': 'application/json' }, withCredentials: false }
      );
      if (response.data.candidates?.[0]?.content?.parts?.[0]?.text) {
        setAiDescription(response.data.candidates[0].content.parts[0].text);
      } else {
        setAiDescription('No AI description available');
      }
    } catch (error) {
      console.error('Failed to fetch AI enhancements:', error);
      setAiDescription('Failed to fetch AI description');
    }
  };

  const handleAtomClick = (viewer, atom) => {
    if (highlightedAtom && highlightedAtom.index === atom.index) {
      viewer.setStyle({ index: atom.index }, { sphere: { scale: 0.3 }, stick: { radius: 0.2 } });
      setHighlightedAtom(null);
    } else {
      viewer.setStyle({ index: atom.index }, { sphere: { scale: 0.5, color: 'yellow' }, stick: { radius: 0.3, color: 'yellow' } });
      setHighlightedAtom({ element: atom.elem, x: atom.x, y: atom.y, z: atom.z, index: atom.index });
    }
    viewer.render();
  };

  // Custom spin animation
  const spinModel = (viewer) => {
    if (!spinAnimation || !viewer) return;
    const animate = () => {
      viewer.rotate(2, 'y'); // Rotate 2 degrees around y-axis
      viewer.render();
      animationFrameId.current = requestAnimationFrame(animate);
    };
    animationFrameId.current = requestAnimationFrame(animate);
  };

  // Enhanced 3D rendering for parent
  useEffect(() => {
    if (showParent3D && parent3DRef.current && structure?.smiles && is3DmolLoaded && window.$3Dmol) {
      parent3DRef.current.innerHTML = '';
      const mol = window.RDKit.get_mol(structure.smiles);
      const molData = mol.get_molblock();
      mol.delete();
      
      const viewer = window.$3Dmol.createViewer(parent3DRef.current, { backgroundColor: 'white' });
      if (!viewer) {
        console.error('Failed to create parent 3D viewer');
        return;
      }
      parentViewerRef.current = viewer;
      viewer.addModel(molData, 'mol');

      const styleConfig = parentStyle === 'ballstick' ? { stick: { radius: 0.2 }, sphere: { scale: 0.3 } } :
                         parentStyle === 'stick' ? { stick: { radius: 0.2 } } :
                         parentStyle === 'sphere' ? { sphere: { scale: 0.5 } } :
                         { surface: { opacity: 0.9, color: 'grey' } };
      viewer.setStyle({}, styleConfig);

      if (showLabels) {
        const model = viewer.getModel(0);
        const atoms = model.selectedAtoms({});
        atoms.forEach((atom, index) => {
          viewer.addLabel(
            `${atom.elem}${index + 1}`,
            { fontSize: 10, backgroundColor: 'black', backgroundOpacity: 0.8 },
            { x: atom.x, y: atom.y, z: atom.z }
          );
        });
      }

      viewer.setClickable({}, true, (atom) => handleAtomClick(viewer, atom));
      viewer.zoomTo();
      viewer.render();

      if (spinAnimation) {
        if (animationFrameId.current) cancelAnimationFrame(animationFrameId.current);
        spinModel(viewer);
      } else if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
        animationFrameId.current = null;
      }

      return () => {
        if (animationFrameId.current) {
          cancelAnimationFrame(animationFrameId.current);
          animationFrameId.current = null;
        }
      };
    }
  }, [showParent3D, structure, is3DmolLoaded, parentStyle, showLabels, spinAnimation]);

  // Enhanced 3D rendering for variant
  useEffect(() => {
    if (showVariant3D && variant3DRef.current && selectedVariant?.smiles && is3DmolLoaded && window.$3Dmol) {
      variant3DRef.current.innerHTML = '';
      const mol = window.RDKit.get_mol(selectedVariant.smiles);
      const molData = mol.get_molblock();
      mol.delete();
      
      const viewer = window.$3Dmol.createViewer(variant3DRef.current, { backgroundColor: 'white' });
      if (!viewer) {
        console.error('Failed to create variant 3D viewer');
        return;
      }
      variantViewerRef.current = viewer;
      viewer.addModel(molData, 'mol');

      const styleConfig = variantStyle === 'ballstick' ? { stick: { radius: 0.2 }, sphere: { scale: 0.3 } } :
                         variantStyle === 'stick' ? { stick: { radius: 0.2 } } :
                         variantStyle === 'sphere' ? { sphere: { scale: 0.5 } } :
                         { surface: { opacity: 0.9, color: 'grey' } };
      viewer.setStyle({}, styleConfig);

      if (showLabels) {
        const model = viewer.getModel(0);
        const atoms = model.selectedAtoms({});
        atoms.forEach((atom, index) => {
          viewer.addLabel(
            `${atom.elem}${index + 1}`,
            { fontSize: 10, backgroundColor: 'black', backgroundOpacity: 0.8 },
            { x: atom.x, y: atom.y, z: atom.z }
          );
        });
      }

      viewer.setClickable({}, true, (atom) => handleAtomClick(viewer, atom));
      viewer.zoomTo();
      viewer.render();

      if (spinAnimation) {
        if (animationFrameId.current) cancelAnimationFrame(animationFrameId.current);
        spinModel(viewer);
      } else if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
        animationFrameId.current = null;
      }

      return () => {
        if (animationFrameId.current) {
          cancelAnimationFrame(animationFrameId.current);
          animationFrameId.current = null;
        }
      };
    }
  }, [showVariant3D, selectedVariant, is3DmolLoaded, variantStyle, showLabels, spinAnimation]);

  const download3DView = (viewerRef, filename) => {
    if (viewerRef.current) {
      viewerRef.current.render(() => {
        const canvas = (parent3DRef.current || variant3DRef.current).querySelector('canvas');
        if (canvas) {
          const pngData = canvas.toDataURL('image/png');
          const link = document.createElement('a');
          link.href = pngData;
          link.download = `${filename}.png`;
          link.click();
        } else {
          console.error('Canvas not found for download');
          alert('Failed to download image: Canvas not found');
        }
      });
    }
  };

  const fetchVariantInfo = async (smiles) => {
    try {
      setLoadingVariantInfo(true);
      setVariantInfoError(null);
      const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;
      if (!GEMINI_API_KEY) throw new Error('Gemini API key not set');
      const response = await axios.post(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: `Provide detailed information about the molecule with SMILES "${smiles}".` }] }] },
        { headers: { 'Content-Type': 'application/json' }, withCredentials: false }
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
    setShowVariant3D(false);
    fetchVariantInfo(variant.smiles);
  };

  const handleTabSwitch = (tab) => {
    setActiveTab(tab);
    setShowParent3D(false);
    setShowVariant3D(false);
    setHighlightedAtom(null);
    setSpinAnimation(false);
    if (parent3DRef.current) parent3DRef.current.innerHTML = '';
    if (variant3DRef.current) variant3DRef.current.innerHTML = '';
    if (animationFrameId.current) {
      cancelAnimationFrame(animationFrameId.current);
      animationFrameId.current = null;
    }
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
            onClick={() => handleTabSwitch('parent')}
          >
            Parent Structure
          </button>
          <button
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'variants'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            }`}
            onClick={() => handleTabSwitch('variants')}
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
              onClick={() => handleTabSwitch('variant')}
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
            <>
              <div 
                ref={parent3DRef}
                className="w-full h-96 mt-4 border rounded relative bg-gray-100"
              >
                {!is3DmolLoaded && <p className="text-center pt-20">Loading 3D viewer...</p>}
              </div>
              {is3DmolLoaded && (
                <div className="mt-2">
                  <div className="flex justify-center space-x-4">
                    <select
                      value={parentStyle}
                      onChange={(e) => setParentStyle(e.target.value)}
                      className="border rounded p-1"
                    >
                      <option value="ballstick">Ball & Stick</option>
                      <option value="stick">Stick</option>
                      <option value="sphere">Sphere</option>
                      <option value="surface">Surface</option>
                    </select>
                    <button
                      onClick={() => setShowLabels(!showLabels)}
                      className="bg-gray-500 hover:bg-gray-600 text-white py-1 px-3 rounded"
                    >
                      {showLabels ? 'Hide Labels' : 'Show Labels'}
                    </button>
                    <button
                      onClick={() => setSpinAnimation(!spinAnimation)}
                      className="bg-teal-500 hover:bg-teal-600 text-white py-1 px-3 rounded"
                    >
                      {spinAnimation ? 'Stop Spin' : 'Spin Model'}
                    </button>
                    <button
                      onClick={() => download3DView(parentViewerRef, `${structure.name}-3D`)}
                      className="bg-purple-500 hover:bg-purple-600 text-white py-1 px-3 rounded"
                    >
                      Download PNG
                    </button>
                  </div>
                  {highlightedAtom && (
                    <div className="mt-2 text-center text-sm text-gray-700 bg-yellow-50 p-2 rounded">
                      <strong>Highlighted Atom:</strong> {highlightedAtom.element} (x: {highlightedAtom.x.toFixed(2)}, y: {highlightedAtom.y.toFixed(2)}, z: {highlightedAtom.z.toFixed(2)})
                    </div>
                  )}
                  {aiDescription && (
                    <div className="mt-4 text-center text-sm text-gray-700 bg-blue-50 p-2 rounded">
                      <strong>AI Insight:</strong> {aiDescription}
                    </div>
                  )}
                </div>
              )}
            </>
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
            <>
              <div 
                ref={variant3DRef}
                className="w-full h-96 mt-4 border rounded relative bg-gray-100"
              >
                {!is3DmolLoaded && <p className="text-center pt-20">Loading 3D viewer...</p>}
              </div>
              {is3DmolLoaded && (
                <div className="mt-2">
                  <div className="flex justify-center space-x-4">
                    <select
                      value={variantStyle}
                      onChange={(e) => setVariantStyle(e.target.value)}
                      className="border rounded p-1"
                    >
                      <option value="ballstick">Ball & Stick</option>
                      <option value="stick">Stick</option>
                      <option value="sphere">Sphere</option>
                      <option value="surface">Surface</option>
                    </select>
                    <button
                      onClick={() => setShowLabels(!showLabels)}
                      className="bg-gray-500 hover:bg-gray-600 text-white py-1 px-3 rounded"
                    >
                      {showLabels ? 'Hide Labels' : 'Show Labels'}
                    </button>
                    <button
                      onClick={() => setSpinAnimation(!spinAnimation)}
                      className="bg-teal-500 hover:bg-teal-600 text-white py-1 px-3 rounded"
                    >
                      {spinAnimation ? 'Stop Spin' : 'Spin Model'}
                    </button>
                    <button
                      onClick={() => download3DView(variantViewerRef, `${structure.name}-variant-3D`)}
                      className="bg-purple-500 hover:bg-purple-600 text-white py-1 px-3 rounded"
                    >
                      Download PNG
                    </button>
                  </div>
                  {highlightedAtom && (
                    <div className="mt-2 text-center text-sm text-gray-700 bg-yellow-50 p-2 rounded">
                      <strong>Highlighted Atom:</strong> {highlightedAtom.element} (x: {highlightedAtom.x.toFixed(2)}, y: {highlightedAtom.y.toFixed(2)}, z: {highlightedAtom.z.toFixed(2)})
                    </div>
                  )}
                  {aiDescription && (
                    <div className="mt-4 text-center text-sm text-gray-700 bg-blue-50 p-2 rounded">
                      <strong>AI Insight:</strong> {aiDescription}
                    </div>
                  )}
                </div>
              )}
            </>
          )}

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2">Property Comparison</h3>
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-gray-50 p-3 rounded">
                <h4 className="text-md font-medium">Parent</h4>
                <div className="text-sm">
                  <p><strong>QED:</strong> {structure.properties?.qed?.toFixed(3) || 'N/A'}</p>
                  <p><strong>LogP:</strong> {structure.properties?.logp?.toFixed(3) || 'N/A'}</p>
                </div>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <h4 className="text-md font-medium">Variant</h4>
                <div className="text-sm">
                  <p><strong>QED:</strong> {selectedVariant.properties?.qed?.toFixed(3) || 'N/A'}</p>
                  <p><strong>LogP:</strong> {selectedVariant.properties?.logp?.toFixed(3) || 'N/A'}</p>
                  <p><strong>Similarity:</strong> {(selectedVariant.similarity * 100).toFixed(1)}%</p>
                </div>
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