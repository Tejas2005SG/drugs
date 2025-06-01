import React, { useState, useEffect } from 'react';
import {
  Container,
  TextField,
  Button,
  Typography,
  Box,
  CircularProgress,
  Alert,
  List,
  ListItem,
  ListItemText,
  Paper,
  Grid,
  Collapse,
  IconButton,
  Tooltip,
} from '@mui/material';
import { ExpandMore, ExpandLess, Science } from '@mui/icons-material';
import axios from 'axios';
import { Check } from 'lucide-react';

const timelineSteps = {
  predict: [
    'Fetching Symptoms',
    'Predicting Diseases',
    'Identifying Target Proteins',
    'Retrieving Ligands',
    'Finalizing Analysis',
  ],
  react: [
    'Fetching Ligands from Database',
    'Detecting Functional Groups',
    'Validating Reactants',
    'Running Chemical Reactions',
    'Finalizing Reaction Outputs',
  ],
};

const loadingMessages = {
  predict: [
    'Processing your symptoms to identify potential diseases...',
    'Analyzing symptom patterns to predict diseases...',
    'Identifying target proteins for therapeutic intervention...',
    'Retrieving ligand structures from the database...',
    'Finalizing disease and target analysis for reaction processing...',
  ],
  react: [
    'Fetching ligand SMILES from MongoDB...',
    'Analyzing ligand structures for functional groups...',
    'Validating reactants for chemical compatibility...',
    'Simulating reactions with RDChiral and RDKit...',
    'Generating and scoring reaction products...',
  ],
};

const RDkit = () => {
  const [symptoms, setSymptoms] = useState('');
  const [loading, setLoading] = useState(false);
  const [reactionLoading, setReactionLoading] = useState(false);
  const [error, setError] = useState('');
  const [result, setResult] = useState(null);
  const [reactionResult, setReactionResult] = useState(null);
  const [expanded, setExpanded] = useState({});
  const [currentStep, setCurrentStep] = useState(0);

  // Handle timeline progress
  useEffect(() => {
    if (!loading && !reactionLoading) return;

    const steps = loading ? timelineSteps.predict : timelineSteps.react;
    const interval = setInterval(() => {
      setCurrentStep((prev) => {
        if (prev < steps.length - 1) {
          return prev + 1;
        }
        return prev;
      });
    }, 2000);

    return () => clearInterval(interval);
  }, [loading, reactionLoading]);

  // Jump to final step when loading finishes
  useEffect(() => {
    if (!loading && !reactionLoading && currentStep > 0) {
      const steps = loading ? timelineSteps.predict : timelineSteps.react;
      setCurrentStep(steps.length - 1);
    }
  }, [loading, reactionLoading]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setResult(null);
    setReactionResult(null);
    setLoading(true);
    setCurrentStep(0);

    try {
      const symptomList = symptoms
        .split(',')
        .map((s) => s.trim())
        .filter((s) => s);

      if (symptomList.length === 0) {
        setError('Please enter at least one symptom.');
        setLoading(false);
        return;
      }

      const diseaseResponse = await axios.post(
        'http://localhost:5000/api/newdrug/predictDisease',
        { symptoms: symptomList },
        {
          headers: { Authorization: `Bearer ${localStorage.getItem('token')}` },
        }
      );

      const proteinResponse = await axios.post(
        'http://localhost:5000/api/newdrug/predictTargetProtein',
        {},
        {
          headers: { Authorization: `Bearer ${localStorage.getItem('token')}` },
        }
      );

      setResult({
        disease: diseaseResponse.data,
        proteins: proteinResponse.data,
      });
    } catch (err) {
      setError(
        err.response?.data?.error || 'Error processing disease/protein prediction.'
      );
    } finally {
      setLoading(false);
    }
  };

  const handleRunReactions = async () => {
    setReactionLoading(true);
    setError('');
    setReactionResult(null);
    setCurrentStep(0);

    try {
      const reactionResponse = await axios.get('http://127.0.0.1:5001/api/react');
      setReactionResult(reactionResponse.data);
    } catch (err) {
      setError(err.response?.data?.error || 'Error running reactions.');
    } finally {
      setReactionLoading(false);
    }
  };

  const handleToggle = (key) => {
    setExpanded((prev) => ({ ...prev, [key]: !prev[key] }));
  };

  const getFunctionalGroups = (ligandSmile) => {
    if (ligandSmile === 'Not applicable') {
      return 'N/A (Protein)';
    }
    if (!reactionResult?.reactants) {
      return 'Pending reaction';
    }
    const reactant = reactionResult.reactants.find((r) => r.smiles === ligandSmile);
    if (!reactant || !Array.isArray(reactant.properties?.functional_groups)) {
      console.warn(`No functional groups found for ligand SMILES: ${ligandSmile}`);
      return 'N/A';
    }
    return reactant.properties.functional_groups.join(', ');
  };

  return (
    <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
      <Paper elevation={3} sx={{ p: 4 }}>
        <Typography variant="h4" gutterBottom>
          Drug Discovery Pipeline
        </Typography>
        <Typography variant="subtitle1" color="textSecondary" gutterBottom>
          Enter symptoms to predict diseases, identify targets, and run optimized ligand reactions.
        </Typography>

        <Box component="form" onSubmit={handleSubmit} sx={{ mt: 2 }}>
          <TextField
            fullWidth
            label="Symptoms (comma-separated)"
            value={symptoms}
            onChange={(e) => setSymptoms(e.target.value)}
            placeholder="e.g., fever, cough, fatigue"
            variant="outlined"
            margin="normal"
          />
          <Button
            type="submit"
            variant="contained"
            color="primary"
            disabled={loading}
            sx={{ mt: 2, mr: 2 }}
          >
            {loading ? <CircularProgress size={24} /> : 'Predict and Analyze'}
          </Button>
          {result && (
            <Button
              variant="contained"
              color="secondary"
              onClick={handleRunReactions}
              disabled={reactionLoading}
              sx={{ mt: 2 }}
            >
              {reactionLoading ? <CircularProgress size={24} /> : 'Run Reactions'}
            </Button>
          )}
        </Box>

        {(loading || reactionLoading) && (
          <Box sx={{ mt: 4, p: 3, bgcolor: 'grey.100', borderRadius: 2, boxShadow: 1 }}>
            <Typography variant="h6" gutterBottom>
              Processing Pipeline
            </Typography>
            <Grid container spacing={2}>
              {/* Left: Vertical Timeline */}
              <Grid item xs={12} md={6}>
                {(loading ? timelineSteps.predict : timelineSteps.react).map((step, index) => (
                  <Box key={index} sx={{ display: 'flex', alignItems: 'center', mb: 2, position: 'relative' }}>
                    {/* Step Dot */}
                    <Box
                      sx={{
                        width: 24,
                        height: 24,
                        borderRadius: '50%',
                        bgcolor: index <= currentStep ? 'primary.main' : 'grey.300',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        transition: 'all 0.5s',
                      }}
                    >
                      {index < currentStep ? (
                        <Check sx={{ color: 'white', fontSize: 16 }} />
                      ) : (
                        <Typography sx={{ color: index <= currentStep ? 'white' : 'grey.600', fontSize: 12 }}>
                          {index + 1}
                        </Typography>
                      )}
                    </Box>
                    {/* Step Label */}
                    <Typography
                      sx={{
                        ml: 2,
                        fontWeight: index === currentStep ? 'bold' : 'normal',
                        color: index === currentStep ? 'primary.main' : 'text.secondary',
                        transition: 'color 0.5s',
                      }}
                      aria-live={index === currentStep ? 'polite' : 'off'}
                    >
                      {step}
                    </Typography>
                    {/* Connecting Line */}
                    {index < (loading ? timelineSteps.predict : timelineSteps.react).length - 1 && (
                      <Box
                        sx={{
                          position: 'absolute',
                          left: 12,
                          top: 24,
                          width: 2,
                          height: 40,
                          bgcolor: index < currentStep ? 'primary.main' : 'grey.300',
                          transition: 'background-color 0.5s',
                        }}
                      />
                    )}
                  </Box>
                ))}
              </Grid>
              {/* Right: Dynamic Message */}
              <Grid item xs={12} md={6} sx={{ display: 'flex', alignItems: 'center' }}>
                <Box sx={{ display: 'flex', alignItems: 'center', animate: { opacity: [0, 1], transition: { duration: 0.5 } } }}>
                  <Science sx={{ mr: 2, color: 'primary.main', animation: 'spin 2s linear infinite' }} />
                  <Typography
                    variant="body1"
                    sx={{ color: 'text.primary', fontWeight: 'medium' }}
                    aria-live="polite"
                  >
                    {loading ? loadingMessages.predict[currentStep] : loadingMessages.react[currentStep]}
                  </Typography>
                </Box>
              </Grid>
            </Grid>
          </Box>
        )}

        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}

        {result && (
          <Box sx={{ mt: 4 }}>
            <Typography variant="h5" gutterBottom>
              Results
            </Typography>

            {/* Disease Prediction */}
            <Paper elevation={2} sx={{ p: 2, mb: 2 }}>
              <Typography variant="h6">Predicted Disease</Typography>
              <List>
                {result.disease.predictedDiseases.map((disease, index) => (
                  <ListItem key={index}>
                    <ListItemText
                      primary={`Disease: ${disease.diseaseName}`}
                      secondary={
                        <>
                          Match: ${disease.DiseaseMatchness}
                          <br />
                          Cautions: ${disease.diseaseCautions.join(', ')}
                        </>
                      }
                    />
                  </ListItem>
                ))}
              </List>
              <Typography variant="subtitle1" sx={{ mt: 2 }}>
                Disease Analysis
              </Typography>
              {Object.entries(result.disease.DiseaseAnalysis).map(([key, items]) => (
                <Box key={key} sx={{ mt: 1 }}>
                  <Box sx={{ display: 'flex', alignItems: 'center' }}>
                    <Typography variant="body1" fontWeight="bold">
                      ${key.replace(/([A-Z])/g, ' $1').trim()}
                    </Typography>
                    <IconButton onClick={() => handleToggle(key)} aria-label={expanded[key] ? 'Collapse' : 'Expand'}>
                      {expanded[key] ? <ExpandLess /> : <ExpandMore />}
                    </IconButton>
                  </Box>
                  <Collapse in={expanded[key]}>
                    <List dense>
                      {items.map((item, idx) => (
                        <ListItem key={idx}>
                          <ListItemText
                            primary={item.summary}
                            secondary={
                              <>
                                Source: ${item.source}
                                <br />
                                {item.url && (
                                  <a href={item.url} target="_blank" rel="noopener noreferrer">
                                    Link
                                  </a>
                                )}
                              </>
                            }
                          />
                        </ListItem>
                      ))}
                    </List>
                  </Collapse>
                </Box>
              ))}
            </Paper>

            {/* Target Proteins and Ligands */}
            <Paper elevation={2} sx={{ p: 2, mb: 2 }}>
              <Typography variant="h6">Target Proteins and Ligands</Typography>
              <Grid container spacing={2}>
                <Grid item xs={12} md={6}>
                  <Typography variant="subtitle1">Target Proteins</Typography>
                  <List dense>
                    {result.proteins.TargetProteins.map((protein, index) => (
                      <ListItem key={index}>
                        <ListItemText
                          primary={protein.proteinName}
                          secondary={
                            <>
                              Function: ${protein.proteinFunction}
                              <br />
                              Description: ${protein.ProtienDiscription}
                              <br />
                              Detailed: ${protein.proteinDetailedDiscription}
                              <br />
                              UniProt ID: ${protein.proteinUniport}
                            </>
                          }
                        />
                      </ListItem>
                    ))}
                  </List>
                </Grid>
                <Grid item xs={12} md={6}>
                  <Typography variant="subtitle1">Target Ligands</Typography>
                  <List dense>
                    {result.proteins.TargetLigands.map((ligand, index) => (
                      <ListItem key={index}>
                        <ListItemText
                          primary={ligand.ligandName}
                          secondary={
                            <>
                              Function: ${ligand.ligandFunction}
                              <br />
                              Description: ${ligand.LigandDiscription}
                              <br />
                              Detailed: ${ligand.LigandDetailedDiscription}
                              <br />
                              SMILES: ${ligand.LigandSmile}
                              <br />
                              Functional Groups: ${getFunctionalGroups(ligand.LigandSmile)}
                              <br />
                              DrugBank ID: ${ligand.ligandDrugBankID}
                              {ligand.LigandSmile === 'Not applicable' && (
                                <Alert severity="warning" sx={{ mt: 1 }}>
                                  This ligand is a protein and cannot be processed for reactions.
                                </Alert>
                              )}
                            </>
                          }
                        />
                      </ListItem>
                    ))}
                  </List>
                </Grid>
              </Grid>
            </Paper>

            {/* Reaction Results */}
            {reactionResult && (
              <Paper elevation={2} sx={{ p: 2 }}>
                <Typography variant="h6">Reaction Results</Typography>
                <Grid container spacing={2}>
                  <Grid item xs={12} md={6}>
                    <Typography variant="subtitle1">Successful Reactions</Typography>
                    {reactionResult.reactionResults.length === 0 ? (
                      <Typography variant="body2" color="textSecondary">
                        No successful reactions found.
                      </Typography>
                    ) : (
                      reactionResult.reactionResults.map((reaction, index) => (
                        <Box key={index} sx={{ mb: 2 }}>
                          <Box sx={{ display: 'flex', alignItems: 'center' }}>
                            <Tooltip title={`Confidence: ${(reaction.confidence * 100).toFixed(1)}%`}>
                              <Typography variant="body1" fontWeight="bold">
                                ${reaction.reactionType}
                              </Typography>
                            </Tooltip>
                            <IconButton
                              onClick={() => handleToggle(`reaction-${index}`)}
                              aria-label={expanded[`reaction-${index}`] ? 'Collapse' : 'Expand'}
                            >
                              {expanded[`reaction-${index}`] ? <ExpandLess /> : <ExpandMore />}
                            </IconButton>
                          </Box>
                          <Collapse in={expanded[`reaction-${index}`]}>
                            <List dense>
                              <ListItem>
                                <ListItemText
                                  primary={`Reactants: ${reaction.reactants.join(', ')}`}
                                  secondary={
                                    <>
                                      Functional Groups: ${reaction.reactantGroups.map((g) => g.join(', ')).join(' | ')}
                                      <br />
                                      Description: ${reaction.description}
                                      <br />
                                      Drug Relevance: ${reaction.drugRelevance}
                                      <br />
                                      Confidence: ${(reaction.confidence * 100).toFixed(1)}%
                                    </>
                                  }
                                />
                              </ListItem>
                              {reaction.productSets.map((set, setIndex) => (
                                <ListItem key={setIndex}>
                                  <ListItemText
                                    primary={`Product Set ${setIndex + 1}`}
                                    secondary={set.map((product, prodIndex) => (
                                      <Box key={prodIndex} sx={{ mb: 1 }}>
                                        SMILES: ${product.smiles}
                                        <br />
                                        MW: ${product.molecular_weight.toFixed(2)}
                                        <br />
                                        LogP: ${product.logP.toFixed(2)}
                                        <br />
                                        TPSA: ${product.tpsa.toFixed(2)}
                                        <br />
                                        H-Donors: ${product.num_h_donors}
                                        <br />
                                        H-Acceptors: ${product.num_h_acceptors}
                                        <br />
                                        Rotatable Bonds: ${product.num_rotatable_bonds}
                                        <br />
                                        Stereochemistry: ${product.has_stereochemistry ? 'Yes' : 'No'}
                                      </Box>
                                    ))}
                                  />
                                </ListItem>
                              ))}
                            </List>
                          </Collapse>
                        </Box>
                      ))
                    )}
                  </Grid>
                  <Grid item xs={12} md={6}>
                    <Typography variant="subtitle1">Failed Reactions</Typography>
                    {reactionResult.failedReactions.length === 0 ? (
                      <Typography variant="body2" color="textSecondary">
                        No failed reactions.
                      </Typography>
                    ) : (
                      <List dense>
                        {reactionResult.failedReactions.map((fail, index) => (
                          <ListItem key={index}>
                            <ListItemText
                              primary={`Reactants: ${fail.reactants.join(', ')}`}
                              secondary={
                                <>
                                  Reaction: ${fail.reactionType || 'N/A'}
                                  <br />
                                  Reason: ${fail.reason}
                                </>
                              }
                            />
                          </ListItem>
                        ))}
                      </List>
                    )}
                  </Grid>
                </Grid>
                <Typography variant="subtitle1" sx={{ mt: 2 }}>
                  Statistics
                </Typography>
                <Box sx={{ pl: 2 }}>
                  <Typography variant="body2">
                    Mean MW: ${reactionResult.statistics.mean_mw.toFixed(2)}
                    <br />
                    Std MW: ${reactionResult.statistics.std_mw.toFixed(2)}
                    <br />
                    Min MW: ${reactionResult.statistics.min_mw.toFixed(2)}
                    <br />
                    Max MW: ${reactionResult.statistics.max_mw.toFixed(2)}
                  </Typography>
                </Box>
              </Paper>
            )}
          </Box>
        )}
      </Paper>
    </Container>
  );
};

export default RDkit;