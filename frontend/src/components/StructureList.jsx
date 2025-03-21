import React from 'react';
import PropTypes from 'prop-types';

const StructureList = ({ structures, onSelect, selected, loading }) => {
  // Defensive check to ensure structures is an array
  const safeStructures = Array.isArray(structures) ? structures : [];

  // Loading state when no structures are present yet
  if (loading && safeStructures.length === 0) {
    return (
      <div className="text-center py-4 text-gray-600">
        Loading structures...
      </div>
    );
  }

  // Empty state when no structures exist
  if (safeStructures.length === 0) {
    return (
      <div className="text-center py-4 text-gray-500">
        No structures yet
      </div>
    );
  }

  return (
    <div className="max-h-96 overflow-y-auto">
      <ul className="divide-y divide-gray-200">
        {safeStructures.map((structure) => (
          <li
            role="listitem" // Added for accessibility
            key={structure._id || Math.random().toString(36).substr(2)} // Fallback key
            className={`py-3 cursor-pointer hover:bg-gray-50 ${
              selected?._id === structure._id ? 'bg-blue-50' : ''
            }`}
            onClick={() => onSelect(structure)}
          >
            <div className="flex items-start">
              <div className="flex-1 min-w-0">
                <p className="text-sm font-medium text-gray-900 truncate">
                  {structure.name || 'Unnamed Structure'} {/* Fallback for name */}
                </p>
                <p className="text-xs text-gray-500 truncate">
                  {structure.created
                    ? new Date(structure.created).toLocaleDateString()
                    : 'Date unknown'} {/* Fallback for date */}
                </p>
              </div>
              <div className="inline-flex items-center text-xs font-semibold text-gray-600">
                {structure.generatedStructures?.length || 0} variants
              </div>
            </div>
          </li>
        ))}
      </ul>
    </div>
  );
};

// PropTypes for type checking
StructureList.propTypes = {
  structures: PropTypes.array.isRequired,
  onSelect: PropTypes.func.isRequired,
  selected: PropTypes.object,
  loading: PropTypes.bool.isRequired,
};

// Memoize the component to prevent unnecessary re-renders
export default React.memo(StructureList);