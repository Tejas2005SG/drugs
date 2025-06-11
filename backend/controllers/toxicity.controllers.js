import axios from 'axios';
import Toxicity from '../models/toxicity.model.js';

const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent';
const PUBCHEM_API_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

// Expanded symptom-to-endpoint mapping
const symptomToEndpointMap = {
  liver: ['hepatotoxicity'],
  heart: ['cardiotoxicity'],
  respiratory: ['respiratory_toxicity'],
  neurological: ['neurotoxicity'],
  kidney: ['nephrotoxicity'],
  cancer: ['carcinogenicity'],
  skin: ['skin_sensitization'],
  thyroid: ['thyroid_toxicity'],
  joint: ['musculoskeletal_toxicity'],
  ear: ['ototoxicity'],
  confusion: ['neurotoxicity'],
  heartbeat: ['cardiotoxicity'],
  nail: ['dermatotoxicity'],
  knuckle: ['musculoskeletal_toxicity'],
  blockage: ['ototoxicity']
};

// Fetch compound data from PubChem
const fetchPubChemData = async (compound) => {
  try {
    const response = await axios.get(`${PUBCHEM_API_URL}/compound/name/${encodeURIComponent(compound)}/JSON`, {
      timeout: 10000
    });
    const data = response.data.PC_Compounds?.[0];
    if (!data) return null;
    return {
      smiles: data.props?.find(prop => prop.urn.label === 'SMILES')?.value.sval,
      molecularWeight: data.props?.find(prop => prop.urn.label === 'Molecular Weight')?.value.fval,
      logP: data.props?.find(prop => prop.urn.label === 'LogP')?.value.fval,
      cid: data.cid
    };
  } catch (error) {
    console.error('PubChem API error:', error.message);
    return null;
  }
};

export const predictToxicityController = async (req, res) => {
  try {
    const { compound, symptoms } = req.body;
    const userId = req.user?._id;

    // Validate inputs
    if (!compound || !symptoms || typeof compound !== 'string' || typeof symptoms !== 'string') {
      return res.status(400).json({ message: 'Valid compound (SMILES or molecule name) and symptoms are required' });
    }
    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    // Normalize symptoms
    const symptomsArray = symptoms
      .split(',')
      .map((s) => s.trim().toLowerCase())
      .filter((s) => s.length > 0);

    if (symptomsArray.length === 0) {
      return res.status(400).json({ message: 'At least one valid symptom is required' });
    }

    // Check for duplicate entry
    const existingEntry = await Toxicity.findOne({
      smiles: compound,
      symptoms,
      userId
    }).sort({ created: -1 });
    if (existingEntry) {
      return res.status(200).json({
        message: 'Prediction already exists',
        result: existingEntry.toxicityResult,
        historyId: existingEntry._id
      });
    }

    // Enhanced Gemini prompt
    const geminiPrompt = `
      Perform an advanced toxicity analysis for the compound "${compound}" (SMILES or molecule name, e.g., "aspirin") proposed for treating symptoms "${symptoms}". Use web searches to fetch data from authentic sources like PubChem (https://pubchem.ncbi.nlm.nih.gov), EPA's TEST (https://www.epa.gov), ToxValDB (https://www.epa.gov), or ATSDR (https://www.atsdr.cdc.gov). Follow these steps:

      1. **Input Validation**
        - Validate "${compound}" as a SMILES string or molecule name using PubChem or ChemSpider.
        - If a molecule name, convert to canonical SMILES (e.g., "aspirin" → "CC(=O)OC1=CC=CC=C1C(=O)O").
        - If invalid, infer a related compound by name similarity, substructure, or symptom context (e.g., "Iodoxil" → iodixanol).
        - Provide error details and suggestions for invalid inputs.

      2. **Molecular Properties**
        - Fetch physicochemical properties: molecular weight, logP, solubility, H-bond donors/acceptors.
        - Cite sources (e.g., PubChem CID).

      3. **Overall Toxicity Profile**
        - Estimate acute toxicity: LD50 (mg/kg, route/species), GHS toxicity class (1–5).
        - Use QSAR models (e.g., EPA’s TEST) if experimental data is unavailable.

      4. **Toxicological Endpoints**
        - Assess hepatotoxicity, cardiotoxicity, neurotoxicity, nephrotoxicity, carcinogenicity, skin sensitization, thyroid_toxicity, ototoxicity, musculoskeletal_toxicity, dermatotoxicity.
        - Prioritize endpoints linked to symptoms "${symptoms}".
        - Provide likelihood (high, medium, low, inactive) and rationale.

      5. **Toxicophores**
        - Identify substructures linked to toxicity (e.g., iodine atoms, epoxides).
        - Note prevalence in related compounds.

      6. **Mechanisms of Toxicity**
        - Describe functional groups and biochemical pathways (e.g., CYP450 inhibition, thyroid disruption).
        - Link to symptoms where possible.

      7. **Side Effects**
        - List adverse effects based on structure and symptoms.
        - Specify organ systems and severity (mild, moderate, severe).

      8. **Structure-Based Concerns**
        - Compare to toxic molecules with similar structures.
        - Highlight substructures causing concern.

      9. **Safety Recommendations**
        - Suggest preclinical tests (e.g., Ames test, hERG assay).
        - Propose structural modifications to reduce toxicity.

      10. **Sources**
        - List all sources used (e.g., PubChem, EPA reports) with URLs or identifiers (format: .web:<number> <source>).

      **Output Format**
      Return a plain JSON string (no Markdown, no code blocks) with:
      {
        "isInputValid": boolean,
        "inputError": string | null,
        "smiles": string | null,
        "moleculeName": string | null,
        "inferredCompound": string | null,
        "inferredSmiles": string | null,
        "errorDetails": { "reason": string, "suggestions": [string] } | null,
        "acuteToxicity": { "LD50": string, "toxicityClass": string },
        "endpoints": {
          "hepatotoxicity": string,
          "carcinogenicity": string,
          "cardiotoxicity": string,
          "neurotoxicity": string,
          "nephrotoxicity": string,
          "skin_sensitization": string,
          "thyroid_toxicity": string,
          "ototoxicity": string,
          "musculoskeletal_toxicity": string,
          "dermatotoxicity": string
        },
        "sideEffects": [{ "name": string, "description": string, "severity": string }],
        "mechanisms": [{ "feature": string, "pathway": string }],
        "structureConcerns": [{ "substructure": string, "similarCompound": string, "concern": string }],
        "safetyRecommendations": [{ "test": string, "modification": string }],
        "qsarAnalysis": {
          "properties": [{ "name": string, "value": string, "implication": string }],
          "prediction": string,
          "symptomContext": [{ "symptom": string, "insight": string }]
        },
        "toxicophoreAnalysis": [{ "substructure": string, "concern": string, "prevalence": string }],
        "sources": [string]
      }

      **Fallback Protocol**
      If "${compound}" is invalid or data is missing:
      - Set "isInputValid" to false, provide "inputError" and "errorDetails".
      - Search for similar compounds by name or substructure.
      - Provide symptom-driven insights based on toxins linked to "${symptoms}" (e.g., iodine compounds for thyroid symptoms).
      - Cite sources for all data.
    `;

    let geminiAnalysis;
    try {
      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: geminiPrompt }] }] },
        { headers: { 'Content-Type': 'application/json' }, timeout: 15000 }
      );

      let analysisText = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
      if (!analysisText) {
        throw new Error('No content returned from Gemini API');
      }

      analysisText = analysisText.replace(/```json\n?/, '').replace(/\n?```/, '').trim();

      try {
        geminiAnalysis = JSON.parse(analysisText);
        console.log('Parsed Gemini analysis:', JSON.stringify(geminiAnalysis, null, 2));
      } catch (parseError) {
        console.error('JSON parsing error:', parseError.message, 'Raw response:', analysisText);
        throw parseError;
      }
    } catch (error) {
      console.error('Gemini API error:', error.message);
      // Fallback to PubChem API
      let pubChemData = await fetchPubChemData(compound);
      let inferredCompound = null;
      if (!pubChemData && compound.toLowerCase().includes('iodoxil')) {
        inferredCompound = 'iodixanol';
        pubChemData = await fetchPubChemData(inferredCompound);
      }

      geminiAnalysis = {
        isInputValid: false,
        inputError: `Failed to fetch data: ${error.message}`,
        smiles: pubChemData?.smiles || null,
        moleculeName: pubChemData ? compound : null,
        inferredCompound: inferredCompound || null,
        inferredSmiles: pubChemData?.smiles || null,
        errorDetails: {
          reason: `Gemini/PubChem error: ${error.message}`,
          suggestions: ['Verify compound name/SMILES', 'Try a known compound (e.g., aspirin)', 'Check API connectivity']
        },
        acuteToxicity: { LD50: 'Unknown', toxicityClass: 'Unknown' },
        endpoints: Object.keys(symptomToEndpointMap).reduce((acc, key) => {
          acc[symptomToEndpointMap[key][0]] = symptomsArray.some(s => s.includes(key)) ? 'Potential risk based on symptoms' : 'Unknown';
          return acc;
        }, {
          hepatotoxicity: 'Unknown',
          carcinogenicity: 'Unknown',
          cardiotoxicity: 'Unknown',
          neurotoxicity: 'Unknown',
          nephrotoxicity: 'Unknown',
          skin_sensitization: 'Unknown',
          thyroid_toxicity: 'Unknown',
          ototoxicity: 'Unknown',
          musculoskeletal_toxicity: 'Unknown',
          dermatotoxicity: 'Unknown'
        }),
        sideEffects: symptomsArray.map(symptom => ({
          name: symptom.replace(/\b\w/g, c => c.toUpperCase()),
          description: `Potential adverse effect related to ${symptom}`,
          severity: 'Unknown'
        })),
        mechanisms: symptomsArray.map(symptom => ({
          feature: `Symptom-based (${symptom})`,
          pathway: `Potential toxicity related to ${symptom} organ system`
        })),
        structureConcerns: [],
        safetyRecommendations: [
          { test: 'Validate compound identity', modification: 'Provide a valid SMILES or known molecule name' }
        ],
        qsarAnalysis: {
          properties: pubChemData ? [
            { name: 'Molecular Weight', value: `${pubChemData.molecularWeight} g/mol`, implication: 'Affects bioavailability' },
            { name: 'LogP', value: pubChemData.logP?.toString() || 'Unknown', implication: 'Indicates hydrophobicity' }
          ] : [],
          prediction: 'Insufficient data for QSAR prediction',
          symptomContext: symptomsArray.map(symptom => ({
            symptom,
            insight: `Symptom may indicate toxicity in ${symptomToEndpointMap[Object.keys(symptomToEndpointMap).find(key => symptom.includes(key))] || 'unknown'} system`
          }))
        },
        toxicophoreAnalysis: [],
        sources: pubChemData ? [`.web:1 PubChem CID: ${pubChemData.cid}`] : []
      };
    }

    // Format toxicity result
    const toxicityResult = {
      smiles: geminiAnalysis.smiles || compound,
      symptoms,
      acuteToxicity: geminiAnalysis.acuteToxicity,
      endpoints: geminiAnalysis.endpoints,
      sideEffects: geminiAnalysis.sideEffects
    };

    // Save to database
    try {
      const toxicityEntry = new Toxicity({
        smiles: geminiAnalysis.smiles || compound,
        symptoms,
        toxicityResult,
        userId,
        geminiAnalysis
      });

      await toxicityEntry.validate();
      await toxicityEntry.save();

      res.status(200).json({
        message: geminiAnalysis.isInputValid ? 'Toxicity predicted successfully' : 'Toxicity analysis performed with limitations',
        result: {
          ...toxicityResult,
          _id: toxicityEntry._id,
          moleculeName: geminiAnalysis.moleculeName,
          inferredCompound: geminiAnalysis.inferredCompound,
          inferredSmiles: geminiAnalysis.inferredSmiles,
          errorDetails: geminiAnalysis.errorDetails,
          mechanisms: geminiAnalysis.mechanisms,
          structureConcerns: geminiAnalysis.structureConcerns,
          safetyRecommendations: geminiAnalysis.safetyRecommendations,
          qsarAnalysis: geminiAnalysis.qsarAnalysis,
          toxicophoreAnalysis: geminiAnalysis.toxicophoreAnalysis,
          sources: geminiAnalysis.sources
        }
      });
    } catch (validationError) {
      console.error('Database validation error:', validationError.message);
      return res.status(500).json({
        message: 'Failed to save toxicity data',
        error: `Database validation failed: ${validationError.message}`
      });
    }
  } catch (error) {
    console.error('Error predicting toxicity:', error.message);
    res.status(500).json({
      message: 'Failed to predict toxicity',
      error: error.message
    });
  }
};

export const getGeminiAnalysis = async (req, res) => {
  try {
    const { compound, symptoms } = req.body;
    const userId = req.user?._id;

    if (!compound || !symptoms) {
      return res.status(400).json({ message: 'Compound and symptoms are required' });
    }
    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    // Normalize symptoms
    const symptomsArray = symptoms
      .split(',')
      .map((s) => s.trim().toLowerCase())
      .filter((s) => s.length > 0);

    // Align prompt with predictToxicityController
    const prompt = `
      Perform an advanced toxicity analysis for the compound "${compound}" (SMILES or molecule name, e.g., "aspirin") proposed for treating symptoms "${symptoms}". Use web searches to fetch data from authentic sources like PubChem (https://pubchem.ncbi.nlm.nih.gov), EPA's TEST (https://www.epa.gov), ToxValDB (https://www.ema.europa.eu), or ATSDR (https://www.atsdr.gov). Follow these steps:

      1. **Input Validation**
        - Validate "${compound}" as a SMILES string or molecule name using PubChem or ChemSpider.
        - If a molecule name, convert to canonical SMILES (e.g., "aspirin" → "CC(=O)OC1=CC=CC=C1C1C(=O)O)O").
        - If invalid, infer a related compound by name similarity, substructure, or symptom context (e.g., "Iodoxil" → iodixanol).
        - Provide error details and suggestions for invalid inputs.

      2. **Molecular Properties**
        - Fetch physicochemical properties: molecular weight, logP, solubility, H-bond donors/acceptors.
        - Cite sources (e.g., PubChem CID).

      3. **Overall Toxicity Profile**
        - Estimate acute toxicity: LD50 (mg/kg, route/species), GHS toxicity class (1–5).
        - Use QSAR models (e.g., EPA’s TEST) if experimental data is unavailable.

      4. **Toxicological Endpoints**
        - Assess hepatotoxicity, cardiotoxicity, neurotoxicity, nephrotoxicity, carcinogenicity, skin sensitization, thyroid_toxicity, ototoxicity, musculoskeletal_toxicity, dermatotoxicity.
        - Prioritize endpoints linked to symptoms "${symptoms}".
        - Provide likelihood (high, medium, low, inactive) and rationale.

      5. **Toxicophores**
        - Identify substructures linked to toxicity (e.g., iodine atoms, epoxides).
        - Note prevalence in related compounds.

      6. **Mechanisms of Toxicity**
        - Describe functional groups and biochemical pathways (e.g., CYP450 inhibition, thyroid disruption).
        - Link to symptoms where possible.

      7. **Side Effects**
        - List adverse effects based on structure and symptoms.
        - Specify organ systems and severity (mild, moderate, severe).

      8. **Structure-Based Concerns**
        - Compare to toxic molecules with similar structures.
        - Highlight substructures causing concern.

      9. **Safety Recommendations**
        - Suggest preclinical tests (e.g., Ames test, hERG assay).
        - Propose structural modifications to reduce toxicity.

      10. **Sources**
        - List all sources used (e.g., PubChem, EPA reports) with URLs or identifiers (format: .web:<number> <source>).

      **Output Format**
      Return a plain JSON string (no Markdown, no code blocks) with:
      {
        "isInputValid": boolean,
        "inputError": string | null,
        "smiles": string | null,
        "moleculeName": string | null,
        "inferredCompound": string | null,
        "inferredSmiles": string | null,
        "errorDetails": { "reason": string, "suggestions": [string] } | null,
        "acuteToxicity": { "LD50": string, "toxicityClass": string },
        "endpoints": {
          "hepatotoxicity": string,
          "carcinogenicity": string,
          "cardiotoxicity": string,
          "neurotoxicity": string,
          "nephrotoxicity": string,
          "skin_sensitization": string,
          "thyroid_toxicity": string,
          "ototoxicity": string,
          "musculoskeletal_toxicity": string,
          "dermatotoxicity": string
        },
        "sideEffects": [{ "name": string, "description": string, "severity": string }],
        "mechanisms": [{ "feature": string, "pathway": string }],
        "structureConcerns": [{ "substructure": string, "similarCompound": string, "concern": string }],
        "safetyRecommendations": [{ "test": string, "modification": string }],
        "qsarAnalysis": {
          "properties": [{ "name": string, "value": string, "implication": string }],
          "prediction": string,
          "symptomContext": [{ "symptom": string, "insight": string }]
        },
        "toxicophoreAnalysis": [{ "substructure": string, "concern": string, "prevalence": string }],
        "sources": [string]
      }

      **Fallback Protocol**
      If "${compound}" is invalid or data is missing:
      - Set "isInputValid" to false, provide "inputError" and "errorDetails".
      - Search for similar compounds by name or substructure.
      - Provide symptom-driven insights based on toxins linked to "${symptoms}" (e.g., iodine compounds for thyroid symptoms).
      - Cite sources for all data.
    `;

    let geminiAnalysis;
    try {
      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: prompt }] }] },
        { headers: { 'Content-Type': 'application/json' }, timeout: 15000 }
      );

      let analysisText = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
      if (!analysisText) {
        throw new Error('No content returned from Gemini API');
      }

      analysisText = analysisText.replace(/```json\n?/, '').replace(/\n?```/, '').trim();

      try {
        geminiAnalysis = JSON.parse(analysisText);
        delete geminiAnalysis.rawText;
        console.log('Parsed Gemini analysis:', JSON.stringify(geminiAnalysis, null, 2));
      } catch (parseError) {
        console.error('JSON parsing error:', parseError.message, 'Raw response:', analysisText);
        throw parseError;
      }
    } catch (error) {
      console.error('Gemini API error:', error.message);
      // Fallback to PubChem API
      let pubChemData = await fetchPubChemData(compound);
      let inferredCompound = null;
      if (!pubChemData && compound.toLowerCase().includes('iodoxil')) {
        inferredCompound = 'iodixanol';
        pubChemData = await fetchPubChemData(inferredCompound);
      }

      geminiAnalysis = {
        isInputValid: false,
        inputError: `Failed to fetch data: ${error.message}`,
        smiles: pubChemData?.smiles || null,
        moleculeName: pubChemData ? compound : null,
        inferredCompound: inferredCompound || null,
        inferredSmiles: pubChemData?.smiles || null,
        errorDetails: {
          reason: `Gemini/PubChem error: ${error.message}`,
          suggestions: ['Verify compound name/SMILES', 'Try a known compound like aspirin', 'Check API connectivity']
        },
        acuteToxicity: { LD50: 'Unknown', toxicityClass: 'Unknown' },
        endpoints: {
          hepatotoxicity: 'Unknown',
          carcinogenicity: 'Unknown',
          cardiotoxicity: 'Unknown',
          neurotoxicity: 'Unknown',
          nephrotoxicity: 'Unknown',
          skin_sensitization: 'Unknown',
          thyroid_toxicity: 'Unknown',
          ototoxicity: 'Unknown',
          musculoskeletal_toxicity: 'Unknown',
          dermatotoxicity: 'Unknown'
        },
        sideEffects: symptomsArray.map(symptom => ({
          name: symptom.replace(/\b\w/g, c => c.toUpperCase()),
          description: `Potential adverse effect related to ${symptom}`,
          severity: 'Unknown'
        })),
        mechanisms: [],
        structureConcerns: [],
        safetyRecommendations: [
          { test: 'Validate compound identity', modification: 'Provide a valid SMILES or molecule name' }
        ],
        qsarAnalysis: {
          properties: pubChemData ? [
            { name: 'Molecular Weight', value: `${pubChemData.molecularWeight} g/mol`, implication: 'Affects bioavailability' },
            { name: 'LogP', value: pubChemData.logP?.toString() || 'Unknown', implication: 'Indicates hydrophobicity' }
          ] : [],
          prediction: 'Insufficient data for QSAR prediction',
          symptomContext: symptomsArray.map(symptom => ({
            symptom: symptom,
            insight: `Symptom may indicate toxicity in ${symptomToEndpointMap[Object.keys(symptomToEndpointMap).find(key => symptom.includes(key))] || 'unknown'} system`
          }))
        },
        toxicophoreAnalysis: [],
        sources: pubChemData ? [`.web:1 PubChem CID: ${pubChemData.cid}`] : []
      };
    }

    res.status(200).json({
      message: 'Gemini analysis generated successfully',
      analysis: geminiAnalysis
    });
  } catch (error) {
    console.error('Error generating Gemini analysis:', error.message);
    res.status(500).json({
      message: 'Failed to generate Gemini analysis',
      error: error.message
    });
  }
};

export const saveGeminiAnalysis = async (req, res) => {
  try {
    const { smiles, symptoms, geminiAnalysis } = req.body;
    const userId = req.user?._id;

    if (!smiles || !symptoms || !geminiAnalysis) {
      return res.status(400).json({ message: 'SMILES, symptoms, and analysis are required' });
    }
    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const toxicityEntry = await Toxicity.findOne({ smiles, symptoms, userId }).sort({ created: -1 });

    if (!toxicityEntry) {
      return res.status(404).json({ message: 'No matching toxicity prediction found' });
    }

    // Validate geminiAnalysis
    const validatedAnalysis = {
      isInputValid: geminiAnalysis.isInputValid || false,
      inputError: geminiAnalysis.inputError || null,
      smiles: geminiAnalysis.smiles || null,
      moleculeName: geminiAnalysis.moleculeName || smiles,
      inferredCompound: geminiAnalysis.inferredCompound || null,
      inferredSmiles: geminiAnalysis.inferredSmiles || null,
      errorDetails: geminiAnalysis.errorDetails || null,
      acuteToxicity: geminiAnalysis.acuteToxicity || { LD50: 'Unknown', toxicityClass: 'Unknown' },
      endpoints: geminiAnalysis.endpoints || {
        hepatotoxicity: 'Unknown',
        carcinogenicity: 'Unknown',
        cardiotoxicity: 'Unknown',
        neurotoxicity: 'Unknown',
        nephrotoxicity: 'Unknown',
        skin_sensitization: 'Unknown',
        thyroid_toxicity: 'Unknown',
        ototoxicity: 'Unknown',
        musculoskeletal_toxicity: 'Unknown',
        dermatotoxicity: 'Unknown'
      },
      sideEffects: geminiAnalysis.sideEffects || [],
      mechanisms: geminiAnalysis.mechanisms || [],
      structureConcerns: geminiAnalysis.structureConcerns || [],
      safetyRecommendations: geminiAnalysis.safetyRecommendations || [],
      qsarAnalysis: geminiAnalysis.qsarAnalysis || {
        properties: [],
        prediction: '',
        symptomContext: []
      },
      toxicophoreAnalysis: geminiAnalysis.toxicophoreAnalysis || [],
      sources: geminiAnalysis.sources || []
    };

    toxicityEntry.geminiAnalysis = validatedAnalysis;
    await toxicityEntry.validate();
    await toxicityEntry.save();

    res.status(200).json({ message: 'Analysis saved successfully', historyId: toxicityEntry._id });
  } catch (error) {
    console.error('Error saving Gemini analysis:', error.message, error);
    res.status(500).json({
      message: 'Failed to save analysis',
      error: `Toxicity validation failed: ${error.message}`
    });
  }
};

export const getToxicityHistory = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const history = await Toxicity.find({ userId }).sort({ created: -1 });

    res.status(200).json({ history });
  } catch (error) {
    console.error('Error fetching toxicity history:', error.message);
    res.status(500).json({
      message: 'Failed to fetch history',
      error: error.message
    });
  }
};