import React, { useState } from 'react';
import { Icon } from '@iconify/react';

const EnhancedGeminiInfo = ({ information, isLoading, error }) => {
  const [expandedSection, setExpandedSection] = useState(null);
  const [copiedSections, setCopiedSections] = useState({});

  const handleCopy = async (text, sectionId) => {
    await navigator.clipboard.writeText(text);
    setCopiedSections(prev => ({ ...prev, [sectionId]: true }));
    setTimeout(() => {
      setCopiedSections(prev => ({ ...prev, [sectionId]: false }));
    }, 2000);
  };

  // Loading state
  if (isLoading) {
    return (
      <div className="animate-pulse bg-white rounded-lg shadow-sm border border-gray-100">
        <div className="flex items-center p-4 border-b border-gray-100">
          <div className="h-4 bg-gray-100 rounded w-32 mr-4"></div>
          <div className="h-4 bg-gray-100 rounded w-24"></div>
        </div>
        <div className="p-4 space-y-4">
          <div className="h-3 bg-gray-100 rounded w-full"></div>
          <div className="h-3 bg-gray-100 rounded w-5/6"></div>
          <div className="h-3 bg-gray-100 rounded w-4/6"></div>
        </div>
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="bg-red-50 border border-red-100 rounded-lg p-4 flex items-start">
        <Icon icon="ph:warning-circle-fill" className="flex-shrink-0 w-5 h-5 text-red-500 mr-3 mt-0.5" />
        <div>
          <h4 className="text-red-800 font-medium mb-1">Error Loading Information</h4>
          <p className="text-red-700 text-sm">{error}</p>
        </div>
      </div>
    );
  }

  // Empty state
  if (!information) {
    return (
      <div className="bg-gray-50 border border-gray-100 rounded-lg p-4 text-center">
        <Icon icon="ph:info-fill" className="w-5 h-5 text-gray-400 mx-auto mb-2" />
        <p className="text-gray-500 text-sm">No detailed information available</p>
      </div>
    );
  }

  // Improved information parsing with better section detection
  const parseInformation = (text) => {
    const sectionRegex = /(\n\s*)?(?=\d+\. |[A-Z][a-zA-Z ]+:\n|## |\*{3,}|[A-Z ]{3,}\n)/g;
    return text.split(sectionRegex)
      .filter(section => section && !section.match(/^\n\s*$/))
      .map((section, index) => {
        const cleanSection = section.replace(/^\n+/, '').trim();
        const [firstLine, ...rest] = cleanSection.split('\n');
        
        const isHeader = firstLine.match(/^(\d+\. |[A-Z][a-zA-Z ]+:|## |\*{3,}|[A-Z ]{3,})/);
        const title = isHeader ? firstLine.replace(/[:#]+$/, '') : `Details ${index + 1}`;
        const content = isHeader ? rest.join('\n') : cleanSection;

        return {
          id: index,
          title: title.replace(/\*/g, '').trim(),
          content: content.trim(),
          icon: getSectionIcon(title),
          color: getSectionColor(title)
        };
      });
  };

  // Get appropriate icon based on section title
  const getSectionIcon = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('property')) return 'ph:flask-fill';
    if (lowerTitle.includes('synthesis')) return 'ph:test-tube-fill';
    if (lowerTitle.includes('medical')) return 'ph:first-aid-kit-fill';
    if (lowerTitle.includes('safety')) return 'ph:warning-fill';
    return 'ph:info-fill';
  };

  // Get color scheme based on section type
  const getSectionColor = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('warning') || lowerTitle.includes('safety')) return 'red';
    if (lowerTitle.includes('medical') || lowerTitle.includes('benefit')) return 'green';
    if (lowerTitle.includes('property') || lowerTitle.includes('chemical')) return 'blue';
    return 'gray';
  };

  const sections = parseInformation(information);

  return (
    <div className="bg-white rounded-lg shadow-sm border border-gray-100 overflow-hidden">
      <div className="bg-gray-50 px-4 py-2.5 border-b border-gray-100">
        <h3 className="text-sm font-semibold text-gray-700 flex items-center">
          <Icon icon="ph:atom-fill" className="w-4 h-4 mr-2 text-blue-500" />
          Compound Analysis Report
        </h3>
      </div>

      <div className="divide-y divide-gray-100">
        {sections.map((section) => (
          <div key={section.id} className="group">
            <button
              onClick={() => setExpandedSection(prev => prev === section.id ? null : section.id)}
              className="w-full flex items-center justify-between px-4 py-3 hover:bg-gray-50 transition-colors"
            >
              <div className="flex items-center space-x-3">
                <div className={`bg-${section.color}-50 p-1.5 rounded-lg`}>
                  <Icon 
                    icon={section.icon} 
                    className={`w-5 h-5 text-${section.color}-600`}
                  />
                </div>
                <h4 className="text-left font-medium text-gray-700">
                  {section.title}
                </h4>
              </div>
              <Icon 
                icon={expandedSection === section.id ? 'ph:caret-up' : 'ph:caret-down'} 
                className="w-5 h-5 text-gray-400 group-hover:text-gray-500"
              />
            </button>

            {(expandedSection === section.id || sections.length === 1) && (
              <div className="px-4 pb-4 pt-2 bg-gray-50 border-t border-gray-100">
                <div className="prose-sm max-w-none relative">
                  <button
                    onClick={() => handleCopy(section.content, section.id)}
                    className="absolute top-0 right-0 p-1.5 hover:bg-gray-200 rounded-lg"
                    title="Copy to clipboard"
                  >
                    <Icon 
                      icon={copiedSections[section.id] ? 'ph:check-fill' : 'ph:copy-simple'} 
                      className={`w-4 h-4 ${copiedSections[section.id] ? 'text-green-500' : 'text-gray-400'}`}
                    />
                  </button>
                  
                  {section.content.split('\n\n').map((paragraph, i) => (
                    <div key={i} className="mb-3 last:mb-0">
                      {paragraph.split('\n').map((line, j) => (
                        <p key={j} className="text-gray-700 mb-2">
                          {line.replace(/(\d+\. )/g, '<br>$1')}
                          {line.match(/^(Molecular Weight|Solubility|LogP):/) && (
                            <span className="ml-2 px-2 py-1 bg-blue-50 text-blue-700 text-xs rounded">
                              Key Property
                            </span>
                          )}
                        </p>
                      ))}
                    </div>
                  ))}

                  {section.content.toLowerCase().includes('warning') && (
                    <div className="mt-3 p-3 bg-red-50 border border-red-100 rounded-lg flex items-start">
                      <Icon icon="ph:warning-fill" className="w-4 h-4 text-red-500 mr-2 mt-0.5" />
                      <span className="text-red-700 text-sm">Important safety information</span>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        ))}
      </div>

      <div className="px-4 py-3 bg-gray-50 border-t border-gray-100">
        <p className="text-xs text-gray-500 flex items-center">
          <Icon icon="ph:info-fill" className="w-3.5 h-3.5 mr-1.5" />
          Information generated by AI - verify critical data with reliable sources
        </p>
      </div>
    </div>
  );
};

export default EnhancedGeminiInfo;