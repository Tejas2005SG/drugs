import React, { useState } from 'react';
import axios from 'axios';

const RDkit = () => {
  const [smiles1, setSmiles1] = useState('');
  const [smiles2, setSmiles2] = useState('');
  const [result, setResult] = useState(null);
  const [error, setError] = useState('');

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setResult(null);

    const apiUrl = 'http://127.0.0.1:5001/api/react'; // Direct Flask call for debugging
    console.log('Sending request to:', apiUrl);
    console.log('Payload:', { smiles1, smiles2 });

    try {
      const response = await axios.post(apiUrl, {
        smiles1,
        smiles2,
      });
      console.log('Response:', response.data);
      setResult(response.data);
    } catch (err) {
      console.error('Error:', err.message, err.response?.data);
      setError(err.response?.data?.error || 'An error occurred while processing the reaction');
    }
  };

  return (
    <div style={{ padding: '20px', fontFamily: 'Arial, sans-serif' }}>
      <h2>Chemical Reaction Simulator for Drug Discovery</h2>
      <form onSubmit={handleSubmit}>
        <div>
          <label>SMILES 1:</label>
          <input
            type="text"
            value={smiles1}
            onChange={(e) => setSmiles1(e.target.value)}
            placeholder="e.g., CCO"
            required
            style={{ margin: '10px', padding: '5px', width: '250px' }}
          />
        </div>
        <div>
          <label>SMILES 2:</label>
          <input
            type="text"
            value={smiles2}
            onChange={(e) => setSmiles2(e.target.value)}
            placeholder="e.g., CC(=O)O"
            required
            style={{ margin: '10px', padding: '5px', width: '250px' }}
          />
        </div>
        <button type="submit" style={{ margin: '10px', padding: '10px', backgroundColor: '#4CAF50', color: 'white', border: 'none' }}>
          Discover All Possible Reactions
        </button>
      </form>

      {error && <p style={{ color: 'red', margin: '10px' }}>{error}</p>}
      {result && (
        <div style={{ margin: '20px' }}>
          <h3>Reaction Analysis</h3>
          <h4>Reactants:</h4>
          <ul>
            {result.reactants.map((reactant, index) => (
              <li key={index}>
                SMILES: {reactant.smiles}, MW: {reactant.properties.molecular_weight.toFixed(3)}, 
                LogP: {reactant.properties.logP.toFixed(2)}, TPSA: {reactant.properties.tpsa.toFixed(2)}
              </li>
            ))}
          </ul>
          <h4>Possible Reactions:</h4>
          {result.reactionResults.length > 0 ? (
            result.reactionResults.map((reaction, index) => (
              <div key={index} style={{ marginBottom: '20px', border: '1px solid #ccc', padding: '10px' }}>
                <p><strong>Reaction {index + 1}: {reaction.reactionType}</strong></p>
                <p><strong>Description:</strong> {reaction.description}</p>
                <p><strong>Drug Relevance:</strong> {reaction.drugRelevance}</p>
                <p><strong>Product Sets:</strong></p>
                {reaction.productSets.map((productSet, setIndex) => (
                  <div key={setIndex} style={{ marginLeft: '20px' }}>
                    <p>Product Set {setIndex + 1}:</p>
                    <ul>
                      {productSet.map((product, pIndex) => (
                        <li key={pIndex}>
                          SMILES: {product.smiles}, MW: {product.molecular_weight.toFixed(3)}, 
                          LogP: {product.logP.toFixed(2)}, TPSA: {product.tpsa.toFixed(2)}, 
                          H-Donors: {product.num_h_donors}, H-Acceptors: {product.num_h_acceptors}
                        </li>
                      ))}
                    </ul>
                  </div>
                ))}
              </div>
            ))
          ) : (
            <p>No valid reactions found.</p>
          )}
          <h4>Statistics:</h4>
          <p>Mean Molecular Weight: {result.statistics.mean_mw.toFixed(3)}</p>
          <p>Std Dev Molecular Weight: {result.statistics.std_mw.toFixed(3)}</p>
          <p>Min Molecular Weight: {result.statistics.min_mw.toFixed(3)}</p>
          <p>Max Molecular Weight: {result.statistics.max_mw.toFixed(3)}</p>
        </div>
      )}
    </div>
  );
};

export default RDkit;