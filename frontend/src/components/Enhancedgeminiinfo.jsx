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
      <div className="animate-pulse bg-white rounded-xl shadow-lg border border-gray-100 overflow-hidden">
        <div className="bg-gradient-to-r from-gray-50 to-gray-100 px-6 py-4 border-b border-gray-100">
          <div className="h-5 bg-gray-200 rounded w-40 mb-1"></div>
          <div className="h-3 bg-gray-100 rounded w-60"></div>
        </div>
        <div className="p-6 space-y-6">
          {[1, 2, 3].map((i) => (
            <div key={i} className="space-y-4">
              <div className="flex items-center space-x-3">
                <div className="h-8 w-8 bg-gray-200 rounded-lg"></div>
                <div className="h-4 bg-gray-200 rounded w-32"></div>
              </div>
              <div className="space-y-2">
                <div className="h-3 bg-gray-100 rounded w-full"></div>
                <div className="h-3 bg-gray-100 rounded w-5/6"></div>
              </div>
            </div>
          ))}
        </div>
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="bg-gradient-to-br from-red-50/80 to-red-100/50 border border-red-200 rounded-xl p-5 flex items-start animate-shake-x animate-duration-300">
        <Icon 
          icon="ph:warning-circle-fill" 
          className="flex-shrink-0 w-6 h-6 text-red-600 mr-3 mt-0.5 animate-pulse" 
        />
        <div>
          <h4 className="text-red-900 font-semibold mb-1 text-lg">Analysis Error</h4>
          <p className="text-red-800/90 text-sm leading-relaxed">{error}</p>
        </div>
      </div>
    );
  }

  // Empty state
  if (!information) {
    return (
      <div className="bg-gradient-to-br from-gray-50/80 to-gray-100/50 border border-gray-200 rounded-xl p-6 text-center">
        <div className="mx-auto mb-4 w-12 h-12 bg-blue-100 rounded-2xl flex items-center justify-center">
          <Icon icon="ph:atom" className="w-6 h-6 text-blue-600/80" />
        </div>
        <h4 className="text-gray-700 font-medium mb-2">No Analysis Available</h4>
        <p className="text-gray-500/90 text-sm max-w-xs mx-auto">
          Detailed compound information will appear here once generated
        </p>
      </div>
    );
  }

  // Information parsing logic
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

  // Icon mapping
  const getSectionIcon = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('property')) return 'ph:flask-fill';
    if (lowerTitle.includes('synthesis')) return 'ph:test-tube-fill';
    if (lowerTitle.includes('medical')) return 'ph:first-aid-kit-fill';
    if (lowerTitle.includes('safety')) return 'ph:warning-fill';
    return 'ph:info-fill';
  };

  // Color mapping
  const getSectionColor = (title) => {
    const lowerTitle = title.toLowerCase();
    if (lowerTitle.includes('warning') || lowerTitle.includes('safety')) return 'red';
    if (lowerTitle.includes('medical') || lowerTitle.includes('benefit')) return 'green';
    if (lowerTitle.includes('property') || lowerTitle.includes('chemical')) return 'blue';
    return 'gray';
  };

  const sections = parseInformation(information);

  return (
    <div className="bg-white rounded-xl shadow-lg border border-gray-200/80 overflow-hidden transition-all duration-200 hover:shadow-xl">
      <div className="bg-gradient-to-r from-blue-600 to-blue-500 px-6 py-4">
        <h3 className="text-sm font-semibold text-white/90 flex items-center space-x-3">
          <Icon 
            icon="ph:atom-fill" 
            className="w-5 h-5 text-white/80" 
          />
          <span>Compound Analysis Report</span>
        </h3>
      </div>

      <div className="divide-y divide-gray-200/60">
        {sections.map((section) => (
          <div 
            key={section.id} 
            className="group transition-all duration-200 hover:bg-gray-50/30"
          >
            <button
              onClick={() => setExpandedSection(prev => prev === section.id ? null : section.id)}
              className="w-full flex items-center justify-between px-6 py-4"
            >
              <div className="flex items-center space-x-4">
                <div className={`bg-${section.color}-100/80 p-2 rounded-xl shadow-sm`}>
                  <Icon 
                    icon={section.icon} 
                    className={`w-6 h-6 text-${section.color}-700/90`}
                  />
                </div>
                <h4 className="text-left font-semibold text-gray-800 text-sm">
                  {section.title}
                </h4>
              </div>
              <Icon 
                icon={expandedSection === section.id ? 'ph:caret-up' : 'ph:caret-down'} 
                className={`w-5 h-5 text-gray-400 transition-transform duration-300 ${
                  expandedSection === section.id ? 'rotate-180' : ''
                }`}
              />
            </button>

            {(expandedSection === section.id || sections.length === 1) && (
              <div className="px-6 pb-5 pt-2 bg-gradient-to-b from-white to-gray-50/60 border-t border-gray-200/50">
                <div className="prose-sm max-w-none relative space-y-4">
                  <button
                    onClick={() => handleCopy(section.content, section.id)}
                    className="absolute top-0 right-0 p-2 hover:bg-gray-200/50 rounded-lg transition-colors duration-200 group/copy"
                  >
                    <Icon 
                      icon={copiedSections[section.id] ? 'ph:check-circle-fill' : 'ph:copy-simple'} 
                      className={`w-5 h-5 transition-all ${
                        copiedSections[section.id] 
                          ? 'text-green-600' 
                          : 'text-gray-500 hover:text-gray-700'
                      }`}
                    />
                  </button>
                  
                  {section.content.split('\n\n').map((paragraph, i) => (
                    <div key={i} className="text-gray-700/90 leading-relaxed space-y-3">
                      {paragraph.split('\n').map((line, j) => (
                        <p key={j} className="flex items-start">
                          {line.startsWith('•') && (
                            <span className="text-blue-500/80 mr-2 mt-1">•</span>
                          )}
                          <span className="flex-1">
                            {line.replace('• ', '')}
                            {line.match(/^(Molecular Weight|Solubility|LogP):/) && (
                              <span className="ml-3 px-2.5 py-1 bg-blue-100 text-blue-700 text-xs font-medium rounded-full">
                                Key Property
                              </span>
                            )}
                          </span>
                        </p>
                      ))}
                    </div>
                  ))}

                  {section.content.toLowerCase().includes('warning') && (
                    <div className="mt-4 p-4 bg-red-50/80 border border-red-200/60 rounded-xl flex items-start space-x-3 animate-pulse-slow">
                      <Icon 
                        icon="ph:warning-fill" 
                        className="w-5 h-5 text-red-600/90 flex-shrink-0 mt-0.5" 
                      />
                      <div>
                        <p className="text-red-800/90 text-sm font-medium mb-1">
                          Safety Notice
                        </p>
                        <p className="text-red-700/90 text-sm leading-relaxed">
                          Critical safety information - handle with proper precautions
                        </p>
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        ))}
      </div>

      <div className="px-6 py-4 bg-gray-100/60 border-t border-gray-200/50">
        <p className="text-xs text-gray-600/90 flex items-center space-x-2">
          <Icon 
            icon="ph:info" 
            className="w-4 h-4 text-gray-500/80" 
          />
          <span>
            AI-generated analysis •{' '}
            <span className="font-medium">Always verify critical data</span>
          </span>
        </p>
      </div>
    </div>
  );
};

export default EnhancedGeminiInfo;