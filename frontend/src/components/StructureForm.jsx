import React, { useState, useEffect } from 'react';
import debounce from 'lodash/debounce'; // Optional: Install lodash for debouncing

const StructureForm = ({ onSubmit, loading }) => {
  const [formData, setFormData] = useState({
    name: '',
    smiles: '',
    algorithm: 'CMA-ES',
    numMolecules: 30,
    propertyName: 'QED',
    minimize: false,
    minSimilarity: 0.3,
    particles: 30,
    iterations: 10,
  });

  const [suggestions, setSuggestions] = useState([]);
  const [showSuggestions, setShowSuggestions] = useState(false);

  // Fetch suggestions from PubChem based on input
  const fetchSuggestions = async (query) => {
    if (!query) {
      setSuggestions([]);
      setShowSuggestions(false);
      return;
    }

    try {
      const response = await fetch(
        `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(query)}/json`
      );
      const data = await response.json();
      if (data.dictionary_terms && data.dictionary_terms.compound) {
        setSuggestions(data.dictionary_terms.compound.slice(0, 5)); // Limit to 5 suggestions
        setShowSuggestions(true);
      } else {
        setSuggestions([]);
        setShowSuggestions(false);
      }
    } catch (error) {
      console.error('Error fetching suggestions:', error);
      setSuggestions([]);
      setShowSuggestions(false);
    }
  };

  // Debounce the fetchSuggestions function to avoid too many API calls
  const debouncedFetchSuggestions = debounce(fetchSuggestions, 300);

  // Handle input change for the name field
  const handleNameChange = (e) => {
    const value = e.target.value;
    setFormData({ ...formData, name: value, smiles: '' }); // Reset SMILES when typing
    debouncedFetchSuggestions(value);
  };

  // Fetch SMILES from PubChem when a suggestion is selected
  const handleSuggestionSelect = async (selectedName) => {
    try {
      const response = await fetch(
        `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(selectedName)}/property/CanonicalSMILES/JSON`
      );
      const data = await response.json();
      if (data.PropertyTable && data.PropertyTable.Properties[0]) {
        const smiles = data.PropertyTable.Properties[0].CanonicalSMILES;
        setFormData({ ...formData, name: selectedName, smiles });
      } else {
        setFormData({ ...formData, name: selectedName, smiles: '' });
      }
    } catch (error) {
      console.error('Error fetching SMILES:', error);
      setFormData({ ...formData, name: selectedName, smiles: '' });
    }
    setShowSuggestions(false);
  };

  const handleChange = (e) => {
    const { name, value, type, checked } = e.target;
    setFormData({
      ...formData,
      [name]: type === 'checkbox' ? checked : value,
    });
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    onSubmit(formData);
  };

  return (
    <form onSubmit={handleSubmit}>
      <div className="mb-4 relative">
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Structure Name
        </label>
        <input
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
          name="name"
          type="text"
          placeholder="e.g., Ibuprofen"
          value={formData.name}
          onChange={handleNameChange}
          onFocus={() => setShowSuggestions(true)}
          autoComplete="off"
        />
        {showSuggestions && suggestions.length > 0 && (
          <ul className="absolute z-10 bg-white border rounded w-full mt-1 max-h-40 overflow-auto">
            {suggestions.map((suggestion) => (
              <li
                key={suggestion}
                className="px-3 py-2 hover:bg-gray-200 cursor-pointer"
                onClick={() => handleSuggestionSelect(suggestion)}
              >
                {suggestion}
              </li>
            ))}
          </ul>
        )}
      </div>

      <div className="mb-4">
        <label className="block text-gray-700 text-sm font-bold mb-2">
          SMILES String
        </label>
        <textarea
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
          name="smiles"
          placeholder="SMILES will be auto-filled after selecting a compound"
          value={formData.smiles}
          onChange={handleChange}
          rows="3"
          required
        />
        <p className="text-xs text-gray-500 mt-1">
          SMILES will be fetched from PubChem or enter manually
        </p>
      </div>

      <div className="grid grid-cols-2 gap-4 mb-4">
        <div>
          <label className="block text-gray-700 text-sm font-bold mb-2">
            Algorithm
          </label>
          <select
            className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
            name="algorithm"
            value={formData.algorithm}
            onChange={handleChange}
          >
            <option value="CMA-ES">CMA-ES</option>
            <option value="SSD">Sampling Standard Deviation</option>
          </select>
        </div>

        <div>
          <label className="block text-gray-700 text-sm font-bold mb-2">
            Property
          </label>
          <select
            className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
            name="propertyName"
            value={formData.propertyName}
            onChange={handleChange}
          >
            <option value="QED">QED</option>
            <option value="logP">logP</option>
          </select>
        </div>
      </div>

      <div className="grid grid-cols-2 gap-4 mb-4">
        <div>
          <label className="block text-gray-700 text-sm font-bold mb-2">
            Number of Molecules
          </label>
          <input
            className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
            name="numMolecules"
            type="number"
            min="1"
            max="100"
            value={formData.numMolecules}
            onChange={handleChange}
          />
        </div>

        <div>
          <label className="block text-gray-700 text-sm font-bold mb-2">
            Min Similarity
          </label>
          <input
            className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
            name="minSimilarity"
            type="number"
            min="0"
            max="1"
            step="0.1"
            value={formData.minSimilarity}
            onChange={handleChange}
          />
        </div>
      </div>

      <div className="flex items-center mb-6">
        <input
          id="minimize"
          name="minimize"
          type="checkbox"
          className="h-4 w-4 text-blue-600 focus:ring-blue-500 border-gray-300 rounded"
          checked={formData.minimize}
          onChange={handleChange}
        />
        <label htmlFor="minimize" className="ml-2 block text-sm text-gray-700">
          Minimize property (instead of maximize)
        </label>
      </div>

      <button
        className={`bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded focus:outline-none focus:shadow-outline w-full ${
          loading ? 'opacity-50 cursor-not-allowed' : ''
        }`}
        type="submit"
        disabled={loading}
      >
        {loading ? 'Generating...' : 'Generate Structure'}
      </button>
    </form>
  );
};

export default StructureForm;