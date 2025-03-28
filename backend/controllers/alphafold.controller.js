import axios from 'axios';
import { Job } from '../models/alphafold.model.js';
import fs from 'fs/promises';
import path from 'path';
import { v4 as uuidv4 } from 'uuid';

const UNIPROT_API = 'https://www.ebi.ac.uk/proteins/api';
const ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/files';
const TIMEOUT = 30000;

const validateUniprotId = (uniprotId) => {
  if (!uniprotId || typeof uniprotId !== 'string') {
    return { valid: false, message: 'UniProt ID must be a non-empty string' };
  }
  const validUniprot = /^[A-Z0-9]{6,10}$/i.test(uniprotId);
  if (!validUniprot) {
    return { valid: false, message: 'Invalid UniProt ID format' };
  }
  return { valid: true };
};

const fetchUniprotSequence = async (uniprotId) => {
  try {
    const response = await axios.get(`${UNIPROT_API}/proteins/${uniprotId}`, { timeout: TIMEOUT });
    const sequence = response.data.sequence;
    if (!sequence) throw new Error('No sequence found in UniProt response');
    return sequence;
  } catch (error) {
    console.error(`Failed to fetch UniProt sequence for ${uniprotId}: ${error.message}`);
    throw new Error(`UniProt fetch failed: ${error.message}`);
  }
};

const fetchAlphaFoldPDB = async (uniprotId) => {
  try {
    const pdbUrl = `${ALPHAFOLD_API}/AF-${uniprotId}-F1-model_v4.pdb`;
    const response = await axios.get(pdbUrl, { responseType: 'text', timeout: TIMEOUT });
    return response.data;
  } catch (error) {
    console.error(`Failed to fetch AlphaFold PDB for ${uniprotId}: ${error.message}`);
    throw new Error(`AlphaFold PDB fetch failed: ${error.message}`);
  }
};

const runAlphaFold = async (job) => {
  try {
    const sequence = await fetchUniprotSequence(job.uniprot_id);
    const workDir = path.join('jobs', job.pythonJobId);
    await fs.mkdir(workDir, { recursive: true });

    const fastaPath = path.join(workDir, 'input.fasta');
    await fs.writeFile(fastaPath, `>query\n${sequence}`);

    const outputDir = path.join(workDir, 'output');
    await fs.mkdir(outputDir, { recursive: true });

    const pdbPath = path.join(outputDir, 'ranked_0.pdb');
    const pdbData = await fetchAlphaFoldPDB(job.uniprot_id);
    await fs.writeFile(pdbPath, pdbData);
    console.log(`PDB file created at: ${pdbPath}`);

    job.status = 'completed';
    job.pdbUrl = `/results/${job.pythonJobId}/output/ranked_0.pdb`;
    await job.save();
    console.log(`Job ${job.pythonJobId} completed successfully`);
  } catch (error) {
    console.error(`Error processing job ${job.pythonJobId}: ${error.message}`);
    job.status = 'failed';
    job.error = error.message;
    await job.save();
  }
};

export const submitPrediction = async (req, res) => {
  let job;
  try {
    const { uniprot_id } = req.body;
    console.log('Received prediction request:', uniprot_id);

    const validation = validateUniprotId(uniprot_id);
    if (!validation.valid) {
      return res.status(400).json({ error: validation.message });
    }

    const jobId = uuidv4();
    job = new Job({
      uniprot_id,
      status: 'pending',
      createdAt: new Date(),
      pythonJobId: jobId,
    });
    await job.save();

    job.status = 'processing';
    await job.save();

    runAlphaFold(job).catch(async (err) => {
      console.error(`Background AlphaFold error: ${err.message}`);
      job.status = 'failed';
      job.error = err.message;
      await job.save();
    });

    res.status(202).json({
      jobId: job._id,
      status: job.status,
      uniprotId: uniprot_id,
    });
  } catch (error) {
    console.error('Submission error:', error);
    if (job) {
      job.status = 'failed';
      job.error = error.message;
      await job.save();
    }
    res.status(500).json({ error: error.message || 'Failed to submit prediction' });
  }
};

export const getUniprotSummary = async (req, res) => {
  try {
    const { uniprotId } = req.params;
    const response = await axios.get(`${UNIPROT_API}/proteins/${uniprotId}`, { timeout: TIMEOUT });
    res.json(response.data);
  } catch (error) {
    console.error(`Failed to fetch UniProt summary: ${error.message}`);
    res.status(error.response?.status || 500).json({
      error: error.message || 'Failed to fetch UniProt summary',
    });
  }
};

export const getUniprotAnnotations = async (req, res) => {
  try {
    const { uniprotId } = req.params;
    const response = await axios.get(`${UNIPROT_API}/features/${uniprotId}`, { timeout: TIMEOUT });
    res.json(response.data);
  } catch (error) {
    console.error(`Failed to fetch UniProt annotations: ${error.message}`);
    res.status(error.response?.status || 500).json({
      error: error.message || 'Failed to fetch annotations',
    });
  }
};

export const getJobStatus = async (req, res) => {
  try {
    const { jobId } = req.params;
    const job = await Job.findById(jobId);
    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }
    res.json({
      jobId: job._id,
      status: job.status,
      uniprotId: job.uniprot_id,
      pdb_url: job.pdbUrl,
      error: job.error,
    });
  } catch (error) {
    console.error('Job status error:', error);
    res.status(500).json({ error: 'Failed to fetch job status' });
  }
};