import React, { useState } from 'react';

const StructureForm = ({ onSubmit, loading }) => {
  // Define a mapping of drug names to SMILES strings
  const drugSmilesMap = {
    Lisuride: '[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34',
    Lenalidomide: 'Nc1cc2c(c(c1)C(=O)N[C@@H]3CCc4c(c(c(c(c4)C3)O)O)O)CC(=O)NC2=O',
    Ensitrelvir: 'CC1=NC(=NN1C)CN2C(=O)N(C(=O)C2(C)C)CC3=NC(=C(C=C3)Cl)NC4=NC=C(C=C4)F',
    Ibuprofen: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
    Floxuridine: 'OC[C@H]1O[C@H](C[C@@H]1O)n2cc(F)c(=O)[nH]c2=O',
    Leflunomide: 'CC1=C(C(=O)NC2=CC=C(C=C2)OC(F)(F)F)C=NO1',
    Tolnaftate: 'CC1=CC=C(C=C1)N(C)C(=S)OC2=C3C=CC=CC3=CC=C2',
    Ifenprodil: 'CC(C)(C1=CC=C(C=C1)O)CCN2CC(C2)C(C3=CC=CC=C3)O',
    Amoxicillin: 'CC1(C(N2C(S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=C(C=C3)O)N)C(=O)O)C',
    Donepezil: 'CC1=CC2=C(C=C1)CC(C2=O)CN3CC4=CC=CC=C4C(C3)CC5=CC=CC=C5',
    Lovastatin: 'CC(C)C1=C[C@H](C[C@@H](O1)CC[C@H]2[C@@H]3CC[C@H]([C@]3(CCC2=O)C)OC(=O)[C@H](C)CC)C',
    Iloperidone: 'CC1=CC2=C(C=C1OC)SC3=C(C(=O)C(CN4CCN(CC4)CCC5=CC=C(C=C5)F)=C(N3)C)O2',
    Dicloxacillin: 'CC1=C(C(=NO1)C2=CC=CC=C2Cl)C(=O)N[C@H]3[C@@H]4N(C3=O)[C@H](C(S4)(C)C)C(=O)O',
    Caffeine: 'CN1C(=O)N(C)C2=C1N=C(N(C2=O)C)N',
    Aspirin: 'CC(=O)OC1=CC=CC=C1C(=O)O',
  };

  const [formData, setFormData] = useState({
    name: '',
    smiles: '',
    algorithm: 'CMA-ES', // Default to CMA-ES
    numMolecules: 30,
    propertyName: 'QED',
    minimize: false,
    minSimilarity: 0.3,
    particles: 30,
    iterations: 10,
  });

  const handleChange = (e) => {
    const { name, value, type, checked } = e.target;
    setFormData({
      ...formData,
      [name]: type === 'checkbox' ? checked : value,
    });
  };

  // Handle drug selection from dropdown
  const handleDrugSelect = (e) => {
    const selectedDrug = e.target.value;
    const smiles = drugSmilesMap[selectedDrug] || ''; // Default to empty if not found
    setFormData({
      ...formData,
      name: selectedDrug, // Optionally set the name to the drug name
      smiles: smiles, // Update the SMILES string
    });
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    onSubmit(formData);
  };

  return (
    <form onSubmit={handleSubmit}>
      <div className="mb-4">
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Structure Name
        </label>
        <select
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
          onChange={handleDrugSelect}
          value={formData.name || ''} // Bind to formData.name
        >
          <option value="">Select a drug...</option>
          {Object.keys(drugSmilesMap).map((drug) => (
            <option key={drug} value={drug}>
              {drug}
            </option>
          ))}
        </select>
      </div>

      <div className="mb-4">
        <label className="block text-gray-700 text-sm font-bold mb-2">
          SMILES String
        </label>
        <textarea
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline"
          name="smiles"
          placeholder="e.g., [H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34"
          value={formData.smiles}
          onChange={handleChange}
          rows="3"
          required
        />
        <p className="text-xs text-gray-500 mt-1">Enter the SMILES representation of your molecule</p>
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