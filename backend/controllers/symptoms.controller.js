import axios from 'axios';

const API_MEDIC_URL = 'https://sandbox-healthservice.priaid.ch/';
const API_KEY = process.env.API_MEDIC_KEY;
const SECRET_KEY = process.env.API_MEDIC_SECRET;

// Get auth token
const getToken = async () => {
  try {
    const response = await axios.post('https://sandbox-authservice.priaid.ch/login', {}, {
      headers: {
        'Authorization': `Bearer ${API_KEY}:${SECRET_KEY}`
      }
    });
    return response.data.Token;
  } catch (error) {
    console.error('Auth error:', error);
    throw new Error('Failed to authenticate with ApiMedic');
  }
};

// Get list of symptoms
export const getSymptoms = async (req, res) => {
  try {
    const token = await getToken();
    const response = await axios.get(`${API_MEDIC_URL}symptoms`, {
      params: {
        token,
        language: 'en-gb'
      }
    });
    res.json(response.data);
  } catch (error) {
    console.error('Symptoms error:', error);
    res.status(500).json({ error: 'Failed to fetch symptoms' });
  }
};

// Get diagnosis based on symptoms
export const getDiagnosis = async (req, res) => {
  try {
    const { symptoms } = req.body;
    if (!symptoms || !Array.isArray(symptoms)) {
      return res.status(400).json({ error: 'Invalid symptoms provided' });
    }

    const token = await getToken();
    const response = await axios.get(`${API_MEDIC_URL}diagnosis`, {
      params: {
        token,
        symptoms: JSON.stringify(symptoms),
        language: 'en-gb',
        gender: 'male',
        year_of_birth: 1990
      }
    });
    res.json(response.data);
  } catch (error) {
    console.error('Diagnosis error:', error);
    res.status(500).json({ error: 'Failed to get diagnosis' });
  }
};