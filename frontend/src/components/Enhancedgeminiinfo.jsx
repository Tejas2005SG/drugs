import React, { useState } from 'react';
import { Icon } from '@iconify/react';

const EnhancedGeminiInfo = ({ information, isLoading, error }) => {
  const [selectedSection, setSelectedSection] = useState('');
  const [copied, setCopied] = useState(false);

  const handleCopy = async (text) => {
    await navigator.clipboard.writeText(text);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  // Loading State
  if (isLoading) {
    return (
      <div className="animate-pulse bg-white rounded-2xl shadow-md border border-gray-100 p-6">
        <div className="flex items-center gap-3 mb-6">
          <div className="h-10 w-10 bg-gray-200 rounded-full"></div>
          <div className="h-6 bg-gray-200 rounded w-1/3"></div>
        </div>
        <div className="space-y-4">
          <div className="h-4 bg-gray-200 rounded w-full"></div>
          <div className="h-4 bg-gray-200 rounded w-5/6"></div>
          <div className="h-4 bg-gray-200 rounded w-3/4"></div>
        </div>
      </div>
    );
  }

  // Error State
  if (error) {
    return (
      <div className="bg-gradient-to-br from-red-50 to-red-100 border border-red-200 rounded-2xl p-6 flex items-start gap-4 transition-all duration-300 hover:shadow-md">
        <Icon 
          icon="ph:warning-circle-fill" 
          className="w-7 h-7 text-red-600 animate-pulse" 
        />
        <div>
          <h4 className="text-red-900 font-semibold text-lg mb-1.5">Analysis Failed</h4>
          <p className="text-red-800 text-sm leading-relaxed">{error}</p>
        </div>
      </div>
    );
  }

  // Empty State
  if (!information) {
    return (
      <div className="bg-gradient-to-br from-gray-50 to-gray-100 border border-gray-200 rounded-2xl p-8 text-center transition-all duration-300 hover:shadow-md">
        <div className="mx-auto mb-5 w-14 h-14 bg-blue-100 rounded-full flex items-center justify-center">
          <Icon icon="ph:atom" className="w-7 h-7 text-blue-600" />
        </div>
        <h4 className="text-gray-800 font-semibold text-lg mb-2">No Analysis Available</h4>
        <p className="text-gray-600 text-sm max-w-sm mx-auto">
          Provide a SMILES string to generate a comprehensive drug analysis
        </p>
      </div>
    );
  }

  // Parse Information
  const parseInformation = (text) => {
    const sectionRegex = /##\s+(.+?)(?=\n##|$)/gs;
    const sections = [];
    let match;

    while ((match = sectionRegex.exec(text)) !== null) {
      const [fullSection, title] = match;
      const content = fullSection.replace(`## ${title}`, '').trim();
      sections.push({
        id: sections.length,
        title: title.trim(),
        content: content.trim(),
        icon: getSectionIcon(title),
        color: getSectionColor(title),
      });
    }

    const conclusionMatch = text.match(/(Development potential assessment|Priority research questions|Recommended assay cascade)/gi);
    if (conclusionMatch) {
      sections.push({
        id: sections.length,
        title: 'Conclusion',
        content: text.split(conclusionMatch[0])[1].trim(),
        icon: 'ph:check-circle-fill',
        color: 'green',
      });
    }

    return sections;
  };

  const getSectionIcon = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('structural')) return 'ph:cube-fill';
    if (lowerTitle.includes('chemical')) return 'ph:flask-fill';
    if (lowerTitle.includes('target')) return 'ph:target-fill';
    if (lowerTitle.includes('admet')) return 'ph:heart-fill';
    if (lowerTitle.includes('developmental')) return 'ph:book-open-fill';
    if (lowerTitle.includes('computational')) return 'ph:cpu-fill';
    if (lowerTitle.includes('conclusion')) return 'ph:check-circle-fill';
    return 'ph:info-fill';
  };

  const getSectionColor = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('structural')) return 'purple';
    if (lowerTitle.includes('chemical')) return 'blue';
    if (lowerTitle.includes('target')) return 'orange';
    if (lowerTitle.includes('admet')) return 'red';
    if (lowerTitle.includes('developmental')) return 'teal';
    if (lowerTitle.includes('computational')) return 'indigo';
    if (lowerTitle.includes('conclusion')) return 'green';
    return 'gray';
  };

  const renderTechnicalContent = (content) => {
    const lines = content.split('\n');
    let isDiagram = false;
    let isTable = false;
    let tableRows = [];

    return lines.map((line, index) => {
      if (line.match(/[\|\-\+]+/)) {
        isDiagram = true;
        return (
          <pre key={index} className="font-mono text-sm bg-gray-50 p-3 rounded-lg border border-gray-100 mt-3">
            {line}
          </pre>
        );
      }

      if (line.includes('|') && line.includes('-')) {
        isTable = true;
        tableRows.push(line.split('|').map(cell => cell.trim()));
      } else if (isTable && line.trim()) {
        tableRows.push(line.split('|').map(cell => cell.trim()));
      }

      if (isTable && (!line.trim() || index === lines.length - 1)) {
        isTable = false;
        const tableContent = (
          <table key={index} className="w-full border-collapse mt-3 bg-white rounded-lg shadow-sm">
            <tbody>
              {tableRows.map((row, rowIndex) => (
                <tr key={rowIndex} className="border-b border-gray-100 last:border-b-0 hover:bg-gray-50">
                  {row.map((cell, cellIndex) => (
                    <td key={cellIndex} className="p-3 text-sm text-gray-700">
                      {cell}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        );
        tableRows = [];
        return tableContent;
      }

      return (
        <p key={index} className="flex items-start text-sm text-gray-700 leading-relaxed">
          {line.startsWith('-') && <span className="text-blue-500 mr-2">•</span>}
          <span className="flex-1">{line.replace(/^- /, '')}</span>
        </p>
      );
    });
  };

  const sections = parseInformation(information);

  return (
    <div className="bg-white rounded-2xl shadow-lg border border-gray-100 overflow-hidden transition-all duration-300 hover:shadow-xl max-w-3xl mx-auto">
      <div className="bg-gradient-to-r from-blue-700 to-indigo-700 px-6 py-4">
        <h3 className="text-lg font-semibold text-white flex items-center gap-3">
          <Icon icon="ph:atom-fill" className="w-6 h-6 text-white/90" />
          <span>Drug Molecule Analysis Report</span>
        </h3>
      </div>

      <div className="p-6 space-y-6">
        {/* Tabbed Navigation */}
        <div className="flex flex-wrap gap-2 justify-center">
          {sections.map((section) => (
            <button
              key={section.id}
              onClick={() => setSelectedSection(section.id)}
              className={`flex items-center gap-2 px-4 py-2 rounded-lg text-sm font-medium transition-all duration-200 ${
                selectedSection === section.id
                  ? `bg-${section.color}-100 text-${section.color}-800 shadow-md`
                  : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
              }`}
            >
              <Icon 
                icon={section.icon} 
                className={`w-5 h-5 ${selectedSection === section.id ? `text-${section.color}-600` : 'text-gray-500'}`}
              />
              {section.title}
            </button>
          ))}
        </div>

        {/* Section Content */}
        {selectedSection !== '' && (
          <div className="relative">
            <div className="bg-gray-100 rounded-3xl p-6 shadow-sm transition-all duration-300 animate-fade-in">
              <div className="flex items-center justify-between mb-5">
                <div className="flex items-center gap-4">
                  <div className={`p-2.5 bg-${sections[selectedSection].color}-100 rounded-full shadow-sm`}>
                    <Icon 
                      icon={sections[selectedSection].icon} 
                      className={`w-6 h-6 text-${sections[selectedSection].color}-600 animate-pulse`}
                    />
                  </div>
                  <h4 className="text-xl font-semibold text-gray-800">
                    {sections[selectedSection].title}
                  </h4>
                </div>
                <button
                  onClick={() => handleCopy(sections[selectedSection].content)}
                  className="p-2 rounded-full hover:bg-gray-100 transition-all duration-200 relative group"
                  title="Copy to clipboard"
                >
                  <Icon 
                    icon={copied ? 'ph:check-circle-fill' : 'ph:copy-simple'} 
                    className={`w-5 h-5 ${copied ? 'text-green-500' : 'text-gray-500 group-hover:text-gray-700'}`}
                  />
                  {copied && (
                    <span className="absolute -top-8 left-1/2 transform -translate-x-1/2 bg-gray-800 text-white text-xs px-2 py-1 rounded">
                      Copied!
                    </span>
                  )}
                </button>
              </div>

              <div className="space-y-4 w-96">
                {renderTechnicalContent(sections[selectedSection].content)}
                <div className="flex gap-2 flex-wrap">
                  {sections[selectedSection].content.includes('[Experimental]') && (
                    <span className="px-2.5 py-1 text-xs font-medium text-green-700 bg-green-100 rounded-full">
                      Experimental Data
                    </span>
                  )}
                  {sections[selectedSection].content.includes('[Predicted]') && (
                    <span className="px-2.5 py-1 text-xs font-medium text-blue-700 bg-blue-100 rounded-full">
                      Predicted Data
                    </span>
                  )}
                  {sections[selectedSection].content.includes('[Hypothetical]') && (
                    <span className="px-2.5 py-1 text-xs font-medium text-orange-700 bg-orange-100 rounded-full">
                      Hypothetical
                    </span>
                  )}
                </div>
              </div>
            </div>
          </div>
        )}
      </div>

      <div className="px-6 py-3.5 bg-gray-50 border-t border-gray-100">
        <p className="text-xs text-gray-600 flex items-center gap-2">
          <Icon icon="ph:info" className="w-4 h-4 text-gray-500" />
          <span>
            Generated by Gemini AI Engine. 
          </span>
        </p>
      </div>
    </div>
  );
};



export default EnhancedGeminiInfo;