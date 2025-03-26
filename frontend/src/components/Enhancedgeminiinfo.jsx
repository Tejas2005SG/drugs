import React, { useState } from 'react';
import { Icon } from '@iconify/react';

const EnhancedGeminiInfo = ({ information, isLoading, error }) => {
  const [isExpanded, setIsExpanded] = useState(false);

  // Loading state
  if (isLoading) {
    return (
      <div className="animate-pulse bg-white rounded-lg shadow-md border border-gray-100 p-4">
        <div className="h-5 bg-gray-200 rounded w-40 mb-1"></div>
        <div className="h-3 bg-gray-100 rounded w-60"></div>
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="bg-red-50 border border-red-200 rounded-lg p-4 flex items-start">
        <Icon 
          icon="ph:warning-circle-fill" 
          className="w-5 h-5 text-red-600 mr-2 mt-0.5" 
        />
        <div>
          <h4 className="text-red-900 font-semibold text-sm">Analysis Error</h4>
          <p className="text-red-800 text-xs">{error}</p>
        </div>
      </div>
    );
  }

  // Empty state
  if (!information || information.includes('Failed to retrieve detailed information')) {
    return (
      <div className="bg-gray-50 border border-gray-200 rounded-lg p-4 text-center">
        <Icon icon="ph:atom" className="w-6 h-6 text-blue-600 mx-auto mb-2" />
        <h4 className="text-gray-700 font-medium text-sm">No Analysis Available</h4>
        <p className="text-gray-500 text-xs">Detailed compound information will appear here once generated</p>
      </div>
    );
  }

  // Parse the information into sections and bullet points
  const parseInformation = (text) => {
    // Split into sections based on numbered headings (e.g., "1. Structural Analysis:")
    const sectionRegex = /(\d+\.\s[A-Za-z\s]+:)/g;
    const sections = text.split(sectionRegex).filter(Boolean);

    const parsedSections = [];
    let currentSection = null;

    for (let i = 0; i < sections.length; i++) {
      const item = sections[i].trim();
      if (item.match(sectionRegex)) {
        // If this is a section title, store it and prepare to collect its content
        if (currentSection) {
          parsedSections.push(currentSection);
        }
        currentSection = { title: item, bulletPoints: [] };
      } else if (currentSection) {
        // Collect bullet points for the current section
        const bulletPoints = item
          .split('\n')
          .map(line => line.trim())
          .filter(line => line.startsWith('-') || line.length > 0) // Include non-empty lines
          .map(line => line.replace(/^-/, '').trim()); // Remove leading dash
        currentSection.bulletPoints.push(...bulletPoints);
      }
    }

    // Push the last section if it exists
    if (currentSection) {
      parsedSections.push(currentSection);
    }

    // Filter out sections with no bullet points
    return parsedSections.filter(section => section.bulletPoints.length > 0);
  };

  const sections = parseInformation(information);

  return (
    <div className="bg-white rounded-lg shadow-md border border-gray-200">
      {/* Dropdown Header */}
      <button
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full flex items-center justify-between px-4 py-3 bg-gray-50 hover:bg-gray-100 transition-colors"
      >
        <div className="flex items-center space-x-2">
          <Icon icon="ph:info-fill" className="w-5 h-5 text-orange-500" />
          <span className="text-sm font-medium text-gray-800">Details <span className='text-xs text-blue-700'>(Powered by Gemini)</span></span>
        </div>
        <Icon
          icon={isExpanded ? 'ph:caret-up' : 'ph:caret-down'}
          className="w-5 h-5 text-gray-500"
        />
      </button>

      {/* Expanded Content */}
      {isExpanded && (
        <div className="px-4 py-3 border-t border-gray-200">
          {sections.length > 0 ? (
            sections.map((section, index) => (
              <div key={index} className="mb-4">
                <h4 className="text-sm font-semibold text-gray-800 mb-2">{section.title}</h4>
                <ul className="space-y-1">
                  {section.bulletPoints.map((point, i) => (
                    <li key={i} className="flex items-start text-sm text-gray-700">
                      <span className="text-orange-500 mr-2">•</span>
                      <span>{point || 'No detailed information available for this section.'}</span>
                    </li>
                  ))}
                </ul>
              </div>
            ))
          ) : (
            <p className="text-sm text-gray-700">
              No detailed information available. Raw content: <br />
              {information}
            </p>
          )}
        </div>
      )}
    </div>
  );
};

export default EnhancedGeminiInfo;