"use client"

import { useState, useEffect } from "react"
import {
  Container,
  TextField,
  Button,
  Typography,
  Box,
  CircularProgress,
  Alert,
  Paper,
  Grid,
  Collapse,
  IconButton,
  Tooltip,
  Menu,
  MenuItem,
  Chip,
  Divider,
  LinearProgress,
  Card,
  CardContent,
  CardHeader,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
} from "@mui/material"
import {
  ExpandMore,
  ExpandLess,
  Science,
  CheckCircle,
  LocalHospital,
  Biotech,
  Psychology,
  Warning,
  Info,
  Visibility,
  RestartAlt,
} from "@mui/icons-material"
import axios from "axios"
import { Radar } from "react-chartjs-2"
import Chart from "chart.js/auto"

const timelineSteps = {
  predict: [
    "Fetching Symptoms",
    "Predicting Diseases",
    "Identifying Target Proteins",
    "Retrieving Ligands",
    "Finalizing Analysis",
  ],
  react: [
    "Fetching Ligands from Database",
    "Detecting Functional Groups",
    "Validating Reactants",
    "Running Chemical Reactions",
    "Computing ADMET Properties",
    "Finalizing Reaction Outputs",
  ],
}

const loadingMessages = {
  predict: [
    "Processing your symptoms to identify potential diseases...",
    "Analyzing symptom patterns to predict diseases...",
    "Identifying target proteins for therapeutic intervention...",
    "Retrieving ligand structures from the database...",
    "Finalizing disease and target analysis for reaction processing...",
  ],
  react: [
    "Fetching ligand SMILES from MongoDB...",
    "Analyzing ligand structures for functional groups...",
    "Validating reactants for chemical compatibility...",
    "Simulating reactions with RDChiral and RDKit...",
    "Computing ADMET properties for reactants and products...",
    "Generating and scoring reaction products...",
  ],
}

const RDkit = () => {
  const [symptoms, setSymptoms] = useState("")
  const [loading, setLoading] = useState(false)
  const [reactionLoading, setReactionLoading] = useState(false)
  const [error, setError] = useState("")
  const [result, setResult] = useState(null)
  const [reactionResult, setReactionResult] = useState(null)
  const [expanded, setExpanded] = useState({})
  const [currentStep, setCurrentStep] = useState(0)
  const [availableReactions, setAvailableReactions] = useState([])
  const [reactionsLoading, setReactionsLoading] = useState(false)
  const [anchorEl, setAnchorEl] = useState(null)
  const [smilesImages, setSmilesImages] = useState({})

  useEffect(() => {
    const fetchAvailableReactions = async () => {
      setReactionsLoading(true)
      try {
        const response = await axios.get("http://127.0.0.1:5001/api/reactions")
        setAvailableReactions(response.data.reactions)
      } catch (err) {
        setError(err.response?.data?.error || "Error fetching available reactions.")
      } finally {
        setReactionsLoading(false)
      }
    }

    fetchAvailableReactions()
  }, [])

  useEffect(() => {
    if (!loading && !reactionLoading) return

    const steps = loading ? timelineSteps.predict : timelineSteps.react
    const interval = setInterval(() => {
      setCurrentStep((prev) => {
        if (prev < steps.length - 1) {
          return prev + 1
        }
        return prev
      })
    }, 2000)

    return () => clearInterval(interval)
  }, [loading, reactionLoading])

  useEffect(() => {
    if (!loading && !reactionLoading) {
      if (error) {
        setCurrentStep(0)
      } else if (currentStep > 0) {
        const steps = loading ? timelineSteps.predict : timelineSteps.react
        setCurrentStep(steps.length - 1)
      }
    }
  }, [loading, reactionLoading, error, currentStep])

  useEffect(() => {
    const fetchSmilesImages = async () => {
      if (!reactionResult) return

      const allSmiles = new Set()

      reactionResult.reactants.forEach((reactant) => {
        if (reactant.smiles && reactant.smiles !== "Not applicable" && reactant.smiles.length > 1) {
          allSmiles.add(reactant.smiles)
        }
      })

      reactionResult.reactionResults.forEach((reaction) => {
        reaction.products.forEach((product) => {
          if (product.smiles && product.smiles.length > 1) {
            allSmiles.add(product.smiles)
          }
        })
      })

      const newImages = { ...smilesImages }
      for (const smiles of allSmiles) {
        if (!newImages[smiles]) {
          try {
            const response = await axios.get(`http://127.0.0.1:5001/api/smiles_to_image?smiles=${encodeURIComponent(smiles)}`)
            newImages[smiles] = response.data.image
          } catch (err) {
            console.warn(`Failed to fetch image for SMILES: ${smiles} - ${err.message}`)
            newImages[smiles] = null
          }
        }
      }
      setSmilesImages(newImages)
    }

    fetchSmilesImages()
  }, [reactionResult])

  const handleSubmit = async (e) => {
    e.preventDefault()
    setError("")
    setResult(null)
    setLoading(true)
    setCurrentStep(0)

    try {
      const symptomsList = symptoms
        .split(",")
        .map((s) => s.trim())
        .filter((s) => s)

      if (symptomsList.length === 0) {
        setError("Please enter at least one symptom.")
        setLoading(false)
        return
      }

      const diseaseResponse = await axios.post(
        "http://localhost:5000/api/newdrug/predictDisease",
        { symptoms: symptomsList },
        {
          headers: { Authorization: `Bearer ${localStorage.getItem("token")}` },
        },
      )

      const proteinResponse = await axios.post(
        "http://localhost:5000/api/newdrug/predictTargetProtein",
        {},
        {
          headers: { Authorization: `Bearer ${localStorage.getItem("token")}` },
        },
      )

      setResult({
        disease: diseaseResponse.data,
        proteins: proteinResponse.data,
      })
    } catch (err) {
      setError(err.response?.data?.error || "Error processing disease/protein prediction.")
      setLoading(false)
      setCurrentStep(0)
    } finally {
      setLoading(false)
    }
  }

  const handleRunReactions = async () => {
    setReactionLoading(true)
    setError("")
    setReactionResult(null)
    setCurrentStep(0)
    setSmilesImages({})

    try {
      const reactionResponse = await axios.get("http://127.0.0.1:5001/api/react?include_admet=true")
      setReactionResult(reactionResponse.data)
    } catch (err) {
      setError(err.response?.data?.error || "Error running reactions.")
      setReactionLoading(false)
      setCurrentStep(0)
    } finally {
      setReactionLoading(false)
    }
  }

  const handleReset = () => {
    setSymptoms("")
    setError("")
    setResult(null)
    setReactionResult(null)
    setLoading(false)
    setReactionLoading(false)
    setCurrentStep(0)
    setExpanded({})
    setAnchorEl(null)
    setSmilesImages({})
  }

  const handleToggle = (key) => {
    setExpanded((prev) => ({ ...prev, [key]: !prev[key] }))
  }

  const getFunctionalGroups = (ligandSmiles) => {
    if (ligandSmiles === "Not applicable") {
      return "N/A (Protein)"
    }
    if (!reactionResult || !reactionResult.reactants) {
      return "Run reactions to analyze functional groups"
    }
    const reactant = reactionResult.reactants.find((r) => r.smiles === ligandSmiles)
    if (!reactant || !Array.isArray(reactant.properties?.functional_groups)) {
      console.warn(`No functional groups found for ligand SMILES: ${ligandSmiles}`)
      return "Not identified"
    }
    return reactant.properties.functional_groups.join(", ") || "None"
  }

  const handleMenuOpen = (event) => {
    setAnchorEl(event.currentTarget)
  }

  const handleMenuClose = () => {
    setAnchorEl(null)
  }

  const deduplicateProducts = (products) => {
    const seenSmiles = new Set()
    return products.filter((product) => {
      if (seenSmiles.has(product.smiles)) {
        return false
      }
      seenSmiles.add(product.smiles)
      return true
    })
  }

  const prepareAdmetGraphData = (reactionIndex) => {
    if (!reactionResult || !reactionResult.reactionResults[reactionIndex]) return null

    const reaction = reactionResult.reactionResults[reactionIndex]
    const products = reaction.products

    const propertiesToShow = ["Non-Toxic", "Soluble", "Bioavailable", "hERG Safe"]
    const datasets = []

    products.forEach((product, idx) => {
      if (!product.admet_properties) return
      const admetSections = product.admet_properties.find((section) => section.section === "ADMET")
      if (!admetSections) return

      const data = propertiesToShow.map((prop) => {
        const property = admetSections.properties.find((p) => p.name === prop)
        return property ? Number(property.prediction) : 0
      })

      datasets.push({
        label: `Product ${idx + 1} (${product.smiles.slice(0, 10)}...)`,
        data,
        backgroundColor: `rgba(255, 99, 132, 0.2)`,
        borderColor: `rgba(255, 99, 132, 1)`,
        borderWidth: 2,
        pointBackgroundColor: `rgba(255, 99, 132, 1)`,
        pointBorderColor: '#fff',
        pointHoverBackgroundColor: '#fff',
        pointHoverBorderColor: `rgba(255, 99, 132, 1)`,
      })
    })

    return {
      labels: propertiesToShow,
      datasets,
    }
  }

  return (
    <Box sx={{ minHeight: "100vh", bgcolor: "#f5f5f5", py: 4 }}>
      <Container maxWidth="xl">
        <Paper elevation={2} sx={{ p: 4, mb: 4, borderRadius: 1 }}>
          <Box sx={{ textAlign: "center", mb: 4 }}>
            <Typography variant="h4" sx={{ fontWeight: 700, color: "#2c3e50", mb: 1 }}>
              Drug Discovery Pipeline
            </Typography>
            <Typography variant="subtitle1" sx={{ color: "#5d6d7e", maxWidth: 800, mx: "auto" }}>
              Advanced molecular analysis system for symptom-based disease prediction, target identification, and
              optimized ligand reaction simulation
            </Typography>
          </Box>

          {/* Input Form */}
          <Card variant="outlined" sx={{ mb: 4, borderRadius: 1 }}>
            <CardContent sx={{ p: 3 }}>
              <Box component="form" onSubmit={handleSubmit}>
                <Grid container spacing={3} alignItems="center">
                  <Grid item xs={12} md={6}>
                    <TextField
                      fullWidth
                      label="Enter Symptoms"
                      value={symptoms}
                      onChange={(e) => setSymptoms(e.target.value)}
                      placeholder="e.g., fever, cough, fatigue, headache"
                      variant="outlined"
                      InputProps={{
                        startAdornment: <LocalHospital sx={{ mr: 1, color: "text.secondary" }} />,
                      }}
                    />
                  </Grid>
                  <Grid item xs={12} md={6}>
                    <Box sx={{ display: "flex", gap: 2, flexWrap: "wrap" }}>
                      <Button
                        type="submit"
                        variant="contained"
                        size="large"
                        disabled={loading || reactionLoading}
                        startIcon={loading ? <CircularProgress size={20} /> : <Psychology />}
                        sx={{
                          bgcolor: "#2c3e50",
                          "&:hover": {
                            bgcolor: "#1a252f",
                          },
                          textTransform: "none",
                        }}
                      >
                        {loading ? "Analyzing..." : "Predict & Analyze"}
                      </Button>
                      {result && (
                        <Button
                          variant="contained"
                          size="large"
                          onClick={handleRunReactions}
                          disabled={loading || reactionLoading}
                          startIcon={reactionLoading ? <CircularProgress size={20} /> : <Biotech />}
                          sx={{
                            bgcolor: "#34495e",
                            "&:hover": {
                              bgcolor: "#2c3e50",
                            },
                            textTransform: "none",
                          }}
                        >
                          {reactionLoading ? "Processing..." : "Run Reactions"}
                        </Button>
                      )}
                      <Button
                        variant="outlined"
                        size="large"
                        onClick={handleReset}
                        disabled={loading || reactionLoading}
                        startIcon={<RestartAlt />}
                        sx={{
                          color: "#2c3e50",
                          borderColor: "#2c3e50",
                          "&:hover": {
                            borderColor: "#1a252f",
                            bgcolor: "rgba(44, 62, 80, 0.04)",
                          },
                          textTransform: "none",
                        }}
                      >
                        Reset
                      </Button>
                    </Box>
                  </Grid>
                </Grid>
              </Box>
            </CardContent>
          </Card>

          {/* Processing Timeline */}
          {(loading || reactionLoading) && (
            <Card variant="outlined" sx={{ mb: 4, borderRadius: 1 }}>
              <CardHeader
                title={
                  <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                    Processing Pipeline
                  </Typography>
                }
                subheader={
                  <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                    Real-time analysis progress tracking
                  </Typography>
                }
              />
              <CardContent sx={{ pt: 0 }}>
                <LinearProgress
                  variant="determinate"
                  value={(currentStep / (loading ? timelineSteps.predict : timelineSteps.react).length) * 100}
                  sx={{
                    mb: 4,
                    height: 6,
                    borderRadius: 3,
                    bgcolor: "#e0e0e0",
                  }}
                />
                <Grid container spacing={4}>
                  <Grid item xs={12} lg={7}>
                    <Box sx={{ position: "relative" }}>
                      {(loading ? timelineSteps.predict : timelineSteps.react).map((step, index) => (
                        <Box key={index} sx={{ display: "flex", alignItems: "center", mb: 3, position: "relative" }}>
                          <Box
                            sx={{
                              width: 36,
                              height: 36,
                              borderRadius: "50%",
                              bgcolor: index <= currentStep ? "#2c3e50" : "#e0e0e0",
                              display: "flex",
                              alignItems: "center",
                              justifyContent: "center",
                              transition: "all 0.3s ease",
                              zIndex: 2,
                              position: "relative",
                            }}
                          >
                            {index < currentStep ? (
                              <CheckCircle sx={{ color: "white", fontSize: 20 }} />
                            ) : index === currentStep ? (
                              <CircularProgress size={18} sx={{ color: "white" }} />
                            ) : (
                              <Typography sx={{ color: "grey.600", fontSize: 14, fontWeight: 600 }}>
                                {index + 1}
                              </Typography>
                            )}
                          </Box>
                          <Box sx={{ ml: 2, flex: 1 }}>
                            <Typography
                              variant="body1"
                              sx={{
                                fontWeight: index === currentStep ? 600 : 400,
                                color: index === currentStep ? "#2c3e50" : index < currentStep ? "#2c3e50" : "#5d6d7e",
                                transition: "all 0.3s ease",
                              }}
                            >
                              {step}
                            </Typography>
                          </Box>
                          {index < (loading ? timelineSteps.predict : timelineSteps.react).length - 1 && (
                            <Box
                              sx={{
                                position: "absolute",
                                left: 18,
                                top: 36,
                                width: 2,
                                height: 40,
                                bgcolor: index < currentStep ? "#2c3e50" : "#e0e0e0",
                                transition: "all 0.3s ease",
                                zIndex: 1,
                              }}
                            />
                          )}
                        </Box>
                      ))}
                    </Box>
                  </Grid>
                  <Grid item xs={12} lg={5}>
                    <Card variant="outlined" sx={{ height: "100%", display: "flex", alignItems: "center" }}>
                      <CardContent sx={{ textAlign: "center", width: "100%" }}>
                        <Science sx={{ fontSize: 40, color: "#2c3e50", mb: 2 }} />
                        <Typography variant="subtitle1" sx={{ color: "#2c3e50", fontWeight: 600, mb: 1 }}>
                          Current Process
                        </Typography>
                        <Typography variant="body2" sx={{ color: "#5d6d7e" }}>
                          {loading ? loadingMessages.predict[currentStep] : loadingMessages.react[currentStep]}
                        </Typography>
                      </CardContent>
                    </Card>
                  </Grid>
                </Grid>
              </CardContent>
            </Card>
          )}

          {/* Error Alert */}
          {error && (
            <Alert
              severity="error"
              sx={{
                mb: 4,
                "& .MuiAlert-icon": {
                  fontSize: 20,
                },
              }}
              onClose={() => setError("")}
              icon={<Warning />}
            >
              <Typography variant="body2" sx={{ fontWeight: 500 }}>
                {error}
              </Typography>
            </Alert>
          )}

          {/* Available Reactions Section */}
          <Card variant="outlined" sx={{ mb: 4, borderRadius: 1 }}>
            <CardHeader
              title={
                <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                  Available Chemical Reactions
                </Typography>
              }
              subheader={
                <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                  Browse supported reaction types and mechanisms
                </Typography>
              }
            />
            <CardContent>
              {reactionsLoading ? (
                <Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
                  <CircularProgress size={24} />
                  <Typography>Loading available reactions...</Typography>
                </Box>
              ) : availableReactions.length > 0 ? (
                <Box>
                  <Button
                    variant="contained"
                    onClick={handleMenuOpen}
                    startIcon={<Visibility />}
                    sx={{
                      bgcolor: "#2c3e50",
                      "&:hover": {
                        bgcolor: "#1a252f",
                      },
                      textTransform: "none",
                    }}
                  >
                    View Available Reactions ({availableReactions.length})
                  </Button>
                  <Menu
                    anchorEl={anchorEl}
                    open={Boolean(anchorEl)}
                    onClose={handleMenuClose}
                    PaperProps={{
                      style: {
                        maxHeight: 500,
                        width: "700px",
                      },
                    }}
                  >
                    {availableReactions.map((reaction, index) => (
                      <MenuItem key={index} onClick={handleMenuClose} sx={{ py: 2 }}>
                        <Box sx={{ width: "100%" }}>
                          <Box sx={{ display: "flex", alignItems: "center", mb: 1 }}>
                            <Chip
                              label={reaction.type}
                              size="small"
                              sx={{
                                mr: 2,
                                bgcolor: "#2c3e50",
                                color: "white",
                                fontWeight: 600,
                              }}
                            />
                            <Chip
                              label={`Priority: ${reaction.priority}`}
                              size="small"
                              variant="outlined"
                              sx={{ borderColor: "#5d6d7e", color: "#5d6d7e" }}
                            />
                          </Box>
                          <Typography variant="body2" sx={{ mb: 1 }}>
                            {reaction.description}
                          </Typography>
                          <Typography
                            variant="body2"
                            sx={{
                              fontFamily: "monospace",
                              bgcolor: "#f5f5f5",
                              p: 1,
                              borderRadius: 1,
                              fontSize: "0.85rem",
                            }}
                          >
                            <strong>SMARTS:</strong> {reaction.smarts}
                          </Typography>
                        </Box>
                      </MenuItem>
                    ))}
                  </Menu>
                </Box>
              ) : (
                <Alert severity="info" icon={<Info />}>
                  <Typography variant="body2">
                    No reactions currently available. Please check your connection or try again later.
                  </Typography>
                </Alert>
              )}
            </CardContent>
          </Card>
        </Paper>

        {/* Results Section - Vertical Layout */}
        {result && (
          <Box sx={{ display: "flex", flexDirection: "column", gap: 3 }}>
            {/* Disease Prediction Results */}
            <Paper elevation={2} sx={{ borderRadius: 1 }}>
              <CardHeader
                title={
                  <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                    Disease Prediction Results
                  </Typography>
                }
                subheader={
                  <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                    AI-powered disease identification and analysis
                  </Typography>
                }
              />
              <CardContent>
                <Box sx={{ mb: 4 }}>
                  <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2, color: "#2c3e50" }}>
                    Predicted Diseases
                  </Typography>
                  {result.disease.predictedDiseases.map((disease, index) => (
                    <Card key={index} variant="outlined" sx={{ mb: 2, borderRadius: 1 }}>
                      <CardContent sx={{ p: 2 }}>
                        <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                          <Typography variant="subtitle1" sx={{ fontWeight: 600, flex: 1 }}>
                            {disease.diseaseName || "Unknown Disease"}
                          </Typography>
                          <Chip
                            label={`Match: ${disease.DiseaseMatchness || "N/A"}`}
                            size="small"
                            sx={{
                              bgcolor: "#e8f5e9",
                              color: "#2e7d32",
                              fontWeight: 500,
                            }}
                          />
                        </Box>
                        {disease.diseaseCautions && disease.diseaseCautions.length > 0 && (
                          <Box>
                            <Typography variant="body2" sx={{ fontWeight: 600, mb: 1, color: "#5d6d7e" }}>
                              Cautions:
                            </Typography>
                            <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                              {disease.diseaseCautions.map((caution, idx) => (
                                <Chip
                                  key={idx}
                                  label={caution}
                                  size="small"
                                  variant="outlined"
                                  sx={{ borderColor: "#d32f2f", color: "#d32f2f" }}
                                />
                              ))}
                            </Box>
                          </Box>
                        )}
                      </CardContent>
                    </Card>
                  ))}
                </Box>

                <Divider sx={{ my: 3 }} />

                <Box>
                  <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2, color: "#2c3e50" }}>
                    Detailed Disease Analysis
                  </Typography>
                  {Object.entries(result.disease.DiseaseAnalysis || {}).map(([key, items]) => (
                    <Card key={key} variant="outlined" sx={{ mb: 2, borderRadius: 1 }}>
                      <CardContent sx={{ p: 2 }}>
                        <Box
                          sx={{ display: "flex", alignItems: "center", cursor: "pointer" }}
                          onClick={() => handleToggle(key)}
                        >
                          <Typography
                            variant="subtitle1"
                            sx={{ fontWeight: 600, flex: 1, textTransform: "capitalize" }}
                          >
                            {key.replace(/([A-Z])/g, " $1").trim()}
                          </Typography>
                          <Chip label={items.length} size="small" sx={{ mr: 1, bgcolor: "#e0e0e0" }} />
                          <IconButton size="small">{expanded[key] ? <ExpandLess /> : <ExpandMore />}</IconButton>
                        </Box>
                        <Collapse in={expanded[key]}>
                          <Box sx={{ mt: 2 }}>
                            {items.map((item, idx) => (
                              <Paper key={idx} variant="outlined" sx={{ p: 2, mb: 2, bgcolor: "#fafafa" }}>
                                <Typography variant="body2" sx={{ fontWeight: 500, mb: 1 }}>
                                  {item.summary || "No summary available"}
                                </Typography>
                                <Box sx={{ display: "flex", alignItems: "center", gap: 2, mt: 1 }}>
                                  <Chip label={`Source: ${item.source || "N/A"}`} size="small" variant="outlined" />
                                  {item.url && (
                                    <Button
                                      size="small"
                                      href={item.url}
                                      target="_blank"
                                      rel="noopener noreferrer"
                                      sx={{ textTransform: "none", color: "#2c3e50" }}
                                    >
                                      View Source
                                    </Button>
                                  )}
                                </Box>
                              </Paper>
                            ))}
                          </Box>
                        </Collapse>
                      </CardContent>
                    </Card>
                  ))}
                </Box>
              </CardContent>
            </Paper>

            {/* Target Proteins */}
            <Paper elevation={2} sx={{ borderRadius: 1 }}>
              <CardHeader
                title={
                  <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                    Target Proteins
                  </Typography>
                }
                subheader={
                  <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                    Identified protein targets for therapeutic intervention
                  </Typography>
                }
              />
              <CardContent>
                {result.proteins.TargetProteins.map((protein, index) => (
                  <Card key={index} variant="outlined" sx={{ mb: 2, borderRadius: 1 }}>
                    <CardContent sx={{ p: 2 }}>
                      <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2 }}>
                        {protein.proteinName || "Unknown Protein"}
                      </Typography>
                      <Grid container spacing={2}>
                        <Grid item xs={12} sm={6}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Function:
                          </Typography>
                          <Typography variant="body2" sx={{ mb: 1 }}>
                            {protein.proteinFunction || "N/A"}
                          </Typography>
                        </Grid>
                        <Grid item xs={12} sm={6}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            UniProt ID:
                          </Typography>
                          <Chip label={protein.proteinUniport || "N/A"} size="small" sx={{ fontFamily: "monospace" }} />
                        </Grid>
                        <Grid item xs={12}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Description:
                          </Typography>
                          <Typography variant="body2" sx={{ mb: 1 }}>
                            {protein.ProtienDiscription || "N/A"}
                          </Typography>
                        </Grid>
                        <Grid item xs={12}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Detailed Information:
                          </Typography>
                          <Typography variant="body2">{protein.proteinDetailedDiscription || "N/A"}</Typography>
                        </Grid>
                      </Grid>
                    </CardContent>
                  </Card>
                ))}
              </CardContent>
            </Paper>

            {/* Target Ligands */}
            <Paper elevation={2} sx={{ borderRadius: 1 }}>
              <CardHeader
                title={
                  <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                    Target Ligands
                  </Typography>
                }
                subheader={
                  <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                    Therapeutic compounds and molecular ligands
                  </Typography>
                }
              />
              <CardContent>
                {result.proteins.TargetLigands.map((ligand, index) => (
                  <Card key={index} variant="outlined" sx={{ mb: 2, borderRadius: 1 }}>
                    <CardContent sx={{ p: 2 }}>
                      <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                        <Typography variant="subtitle1" sx={{ fontWeight: 600, flex: 1 }}>
                          {ligand.ligandName || "Unknown Ligand"}
                        </Typography>
                        {ligand.LigandSmile === "Not applicable" && (
                          <Chip label="Protein" size="small" sx={{ bgcolor: "#fff3e0", color: "#e65100" }} />
                        )}
                      </Box>

                      <Grid container spacing={2}>
                        <Grid item xs={12} sm={6}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Function:
                          </Typography>
                          <Typography variant="body2" sx={{ mb: 1 }}>
                            {ligand.ligandFunction || "N/A"}
                          </Typography>
                        </Grid>
                        <Grid item xs={12} sm={6}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            DrugBank ID:
                          </Typography>
                          <Chip
                            label={ligand.ligandDrugBankID || "N/A"}
                            size="small"
                            sx={{ fontFamily: "monospace" }}
                          />
                        </Grid>
                        <Grid item xs={12}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            SMILES:
                          </Typography>
                          <Paper
                            variant="outlined"
                            sx={{
                              p: 1,
                              bgcolor: "#fafafa",
                              fontFamily: "monospace",
                              fontSize: "0.85rem",
                              wordBreak: "break-all",
                            }}
                          >
                            {ligand.LigandSmile || "N/A"}
                          </Paper>
                        </Grid>
                        {ligand.LigandSmile && ligand.LigandSmile !== "Not applicable" && smilesImages[ligand.LigandSmile] ? (
                          <Grid item xs={12}>
                            <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                              Structure:
                            </Typography>
                            <img
                              src={`data:image/png;base64,${smilesImages[ligand.LigandSmile]}`}
                              alt={`Structure of ${ligand.ligandName}`}
                              style={{ maxWidth: "200px", marginTop: "8px" }}
                            />
                          </Grid>
                        ) : (
                          <Grid item xs={12}>
                            <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                              Structure:
                            </Typography>
                            <Typography variant="caption" sx={{ color: "text.secondary" }}>
                              {ligand.LigandSmile === "Not applicable"
                                ? "Structure not applicable (Protein)"
                                : "Unable to load structure image"}
                            </Typography>
                          </Grid>
                        )}
                        <Grid item xs={12}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Functional Groups:
                          </Typography>
                          <Typography variant="body2" sx={{ mb: 1 }}>
                            {getFunctionalGroups(ligand.LigandSmile)}
                          </Typography>
                        </Grid>
                        <Grid item xs={12}>
                          <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                            Description:
                          </Typography>
                          <Typography variant="body2">{ligand.LigandDiscription || "N/A"}</Typography>
                        </Grid>
                        {ligand.LigandSmile !== "Not applicable" && reactionResult && reactionResult.reactants && (
                          <Grid item xs={12}>
                            <Typography variant="body2" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                              ADMET Properties:
                            </Typography>
                            {(() => {
                              const reactant = reactionResult.reactants.find((r) => r.smiles === ligand.LigandSmile)
                              if (!reactant || !reactant.admet_properties || reactant.admet_properties.length === 0) {
                                return <Typography variant="body2">Not available</Typography>
                              }
                              return reactant.admet_properties.map((section, secIdx) => (
                                <Box key={secIdx} sx={{ mt: 2 }}>
                                  <Typography variant="subtitle2" sx={{ fontWeight: 600, mb: 1 }}>
                                    {section.section}
                                  </Typography>
                                  <Table size="small">
                                    <TableHead>
                                      <TableRow>
                                        <TableCell sx={{ fontWeight: 600 }}>Property</TableCell>
                                        <TableCell sx={{ fontWeight: 600 }}>Prediction</TableCell>
                                        <TableCell sx={{ fontWeight: 600 }}>Units</TableCell>
                                      </TableRow>
                                    </TableHead>
                                    <TableBody>
                                      {section.properties.map((prop, propIdx) => (
                                        <TableRow key={propIdx}>
                                          <TableCell>{prop.name}</TableCell>
                                          <TableCell>{prop.prediction}</TableCell>
                                          <TableCell>{prop.units || "-"}</TableCell>
                                        </TableRow>
                                      ))}
                                    </TableBody>
                                  </Table>
                                </Box>
                              ))
                            })()}
                          </Grid>
                        )}
                      </Grid>

                      {ligand.LigandSmile === "Not applicable" && (
                        <Alert severity="warning" sx={{ mt: 2 }} icon={<Warning />}>
                          <Typography variant="body2">
                            This ligand is a protein and cannot be processed for chemical reactions.
                          </Typography>
                        </Alert>
                      )}
                    </CardContent>
                  </Card>
                ))}
              </CardContent>
            </Paper>

            {/* Chemical Reaction Results */}
            {reactionResult && (
              <Paper elevation={2} sx={{ borderRadius: 1 }}>
                <CardHeader
                  title={
                    <Typography variant="h6" sx={{ fontWeight: 600, color: "#2c3e50" }}>
                      Chemical Reaction Results
                    </Typography>
                  }
                  subheader={
                    <Typography variant="body2" sx={{ color: "text.secondary", mt: 0.5 }}>
                      Computational chemistry analysis and product generation
                    </Typography>
                  }
                />
                <CardContent>
                  <Grid container spacing={3}>
                    {/* Successful Reactions */}
                    <Grid item xs={12} lg={8}>
                      <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2, color: "#2c3e50" }}>
                        Successful Reactions
                      </Typography>
                      {reactionResult.reactionResults.length === 0 ? (
                        <Alert severity="info" icon={<Info />}>
                          <Typography variant="body2">
                            No successful reactions were found with the current parameters.
                          </Typography>
                        </Alert>
                      ) : (
                        reactionResult.reactionResults.map((reaction, index) => {
                          const mainProducts = deduplicateProducts(reaction.products.slice(0, 1))
                          return (
                            <Card key={index} variant="outlined" sx={{ mb: 3, borderRadius: 1 }}>
                              <CardContent sx={{ p: 2 }}>
                                <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                                  <Typography variant="subtitle1" sx={{ fontWeight: 600, flex: 1 }}>
                                    {reaction.reactionType}
                                  </Typography>
                                  <Tooltip title={`Confidence: ${(reaction.confidence * 100).toFixed(1)}%`}>
                                    <Chip
                                      label={`${(reaction.confidence * 100).toFixed(1)}% Confidence`}
                                      size="small"
                                      sx={{
                                        bgcolor:
                                          reaction.confidence > 0.7
                                            ? "#e8f5e9"
                                            : reaction.confidence > 0.4
                                              ? "#fff3e0"
                                              : "#ffebee",
                                        color:
                                          reaction.confidence > 0.7
                                            ? "#2e7d32"
                                            : reaction.confidence > 0.4
                                              ? "#e65100"
                                              : "#c62828",
                                        fontWeight: 500,
                                      }}
                                    />
                                  </Tooltip>
                                  <IconButton size="small" onClick={() => handleToggle(`reaction-${index}`)}>
                                    {expanded[`reaction-${index}`] ? <ExpandLess /> : <ExpandMore />}
                                  </IconButton>
                                </Box>

                                <Typography variant="body2" sx={{ mb: 2, color: "#5d6d7e" }}>
                                  {reaction.description}
                                </Typography>

                                <Box sx={{ mb: 2 }}>
                                  <Typography variant="body2" sx={{ fontWeight: 600, mb: 1, color: "#5d6d7e" }}>
                                    Reactants:
                                  </Typography>
                                  <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                                    {reaction.reactants.map((reactant, idx) => (
                                      <Box key={idx} sx={{ mb: 1 }}>
                                        <Chip
                                          label={reactant}
                                          size="small"
                                          variant="outlined"
                                          sx={{ fontFamily: "monospace" }}
                                        />
                                        {smilesImages[reactant] ? (
                                          <img
                                            src={`data:image/png;base64,${smilesImages[reactant]}`}
                                            alt={`Reactant ${idx + 1}`}
                                            style={{ maxWidth: "150px", marginTop: "8px", marginLeft: "8px" }}
                                          />
                                        ) : (
                                          <Typography variant="caption" sx={{ color: "text.secondary", ml: 1 }}>
                                            Unable to load structure image
                                          </Typography>
                                        )}
                                      </Box>
                                    ))}
                                  </Box>
                                </Box>

                                <Collapse in={expanded[`reaction-${index}`]}>
                                  <Divider sx={{ my: 2 }} />
                                  <Typography variant="subtitle2" sx={{ fontWeight: 600, mb: 2 }}>
                                    Main Product
                                  </Typography>
                                  <Grid container spacing={2}>
                                    {mainProducts.map((product, prodIndex) => (
                                      <Grid item xs={12} md={6} key={prodIndex}>
                                        <Paper variant="outlined" sx={{ p: 2, borderRadius: 1, bgcolor: "#fafafa" }}>
                                          <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                                            <Typography variant="subtitle2" sx={{ fontWeight: 600, flex: 1 }}>
                                              Product {prodIndex + 1}
                                            </Typography>
                                          </Box>
                                          <Grid container spacing={1}>
                                            <Grid item xs={12}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                SMILES:
                                              </Typography>
                                              <Paper
                                                variant="outlined"
                                                sx={{
                                                  p: 1,
                                                  mt: 0.5,
                                                  bgcolor: "white",
                                                  fontFamily: "monospace",
                                                  fontSize: "0.8rem",
                                                  wordBreak: "break-all",
                                                }}
                                              >
                                                {product.smiles}
                                              </Paper>
                                            </Grid>
                                            <Grid item xs={12}>
                                              {smilesImages[product.smiles] ? (
                                                <img
                                                  src={`data:image/png;base64,${smilesImages[product.smiles]}`}
                                                  alt={`Product ${prodIndex + 1}`}
                                                  style={{ maxWidth: "200px", marginTop: "8px" }}
                                                />
                                              ) : (
                                                <Typography variant="caption" sx={{ color: "text.secondary" }}>
                                                  Unable to load structure image
                                                </Typography>
                                              )}
                                            </Grid>
                                            <Grid item xs={6}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                Molecular Weight:
                                              </Typography>
                                              <Typography variant="body2">
                                                {product.molecular_weight.toFixed(2)} g/mol
                                              </Typography>
                                            </Grid>
                                            <Grid item xs={6}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                LogP:
                                              </Typography>
                                              <Typography variant="body2">{product.logP.toFixed(2)}</Typography>
                                            </Grid>
                                            <Grid item xs={6}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                TPSA:
                                              </Typography>
                                              <Typography variant="body2">{product.tpsa.toFixed(2)} Å²</Typography>
                                            </Grid>
                                            <Grid item xs={6}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                H-Donors/Acceptors:
                                              </Typography>
                                              <Typography variant="body2">
                                                {product.num_h_donors}/{product.num_h_acceptors}
                                              </Typography>
                                            </Grid>
                                            <Grid item xs={12}>
                                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                                Functional Groups:
                                              </Typography>
                                              <Box sx={{ display: "flex", flexWrap: "wrap", gap: 0.5, mt: 0.5 }}>
                                                {product.functional_groups.length > 0 ? (
                                                  product.functional_groups.map((group, groupIdx) => (
                                                    <Chip
                                                      key={groupIdx}
                                                      label={group}
                                                      size="small"
                                                      variant="outlined"
                                                      sx={{ fontSize: "0.7rem" }}
                                                    />
                                                  ))
                                                ) : (
                                                  <Typography variant="body2" sx={{ color: "text.secondary" }}>
                                                    None identified
                                                  </Typography>
                                                )}
                                              </Box>
                                            </Grid>
                                            {product.admet_properties && product.admet_properties.length > 0 && (
                                              <Grid item xs={12}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{ fontWeight: 600, color: "#5d6d7e", mt: 2 }}
                                                >
                                                  ADMET Properties:
                                                </Typography>
                                                {product.admet_properties.map((section, secIdx) => (
                                                  <Box key={secIdx} sx={{ mt: 2 }}>
                                                    <Typography variant="subtitle2" sx={{ fontWeight: 600, mb: 1 }}>
                                                      {section.section}
                                                    </Typography>
                                                    <Table size="small">
                                                      <TableHead>
                                                        <TableRow>
                                                          <TableCell sx={{ fontWeight: 600 }}>Property</TableCell>
                                                          <TableCell sx={{ fontWeight: 600 }}>Prediction</TableCell>
                                                          <TableCell sx={{ fontWeight: 600 }}>Units</TableCell>
                                                        </TableRow>
                                                      </TableHead>
                                                      <TableBody>
                                                        {section.properties.map((prop, propIdx) => (
                                                          <TableRow key={propIdx}>
                                                            <TableCell>{prop.name}</TableCell>
                                                            <TableCell>{prop.prediction}</TableCell>
                                                            <TableCell>{prop.units || "-"}</TableCell>
                                                          </TableRow>
                                                        ))}
                                                      </TableBody>
                                                    </Table>
                                                  </Box>
                                                ))}
                                              </Grid>
                                            )}
                                          </Grid>
                                        </Paper>
                                      </Grid>
                                    ))}
                                  </Grid>
                                  {/* ADMET Graph (Updated to Radar Chart) */}
                                  <Box sx={{ mt: 4 }}>
                                    <Typography variant="subtitle2" sx={{ fontWeight: 600, mb: 2 }}>
                                      ADMET Properties Comparison
                                    </Typography>
                                    {prepareAdmetGraphData(index) ? (
                                      <Radar
                                        data={prepareAdmetGraphData(index)}
                                        options={{
                                          responsive: true,
                                          plugins: {
                                            legend: {
                                              position: "top",
                                            },
                                            title: {
                                              display: true,
                                              text: "Blood-Brain Barrier Safe",
                                            },
                                            tooltip: {
                                              callbacks: {
                                                label: (tooltipItem) => {
                                                  return `${tooltipItem.dataset.label}: ${tooltipItem.raw}`;
                                                },
                                              },
                                            },
                                          },
                                          scales: {
                                            r: {
                                              beginAtZero: true,
                                              min: 0,
                                              max: 100,
                                              ticks: {
                                                stepSize: 25,
                                              },
                                              pointLabels: {
                                                font: {
                                                  size: 12,
                                                },
                                              },
                                              grid: {
                                                color: "rgba(0, 0, 0, 0.1)",
                                              },
                                            },
                                          },
                                        }}
                                        height={100}
                                      />
                                    ) : (
                                      <Typography variant="body2" sx={{ color: "text.secondary" }}>
                                        No ADMET data available for comparison.
                                      </Typography>
                                    )}
                                  </Box>
                                </Collapse>
                              </CardContent>
                            </Card>
                          )
                        })
                      )}
                    </Grid>

                    {/* Failed Reactions & Statistics */}
                    <Grid item xs={12} lg={4}>
                      <Box sx={{ mb: 4 }}>
                        <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2, color: "#2c3e50" }}>
                          Failed Reactions
                        </Typography>
                        {reactionResult.failedReactions.length === 0 ? (
                          <Alert severity="success" icon={<CheckCircle />}>
                            <Typography variant="body2">All reactions completed successfully!</Typography>
                          </Alert>
                        ) : (
                          reactionResult.failedReactions.map((fail, index) => (
                            <Card key={index} variant="outlined" sx={{ mb: 2, borderRadius: 1 }}>
                              <CardContent sx={{ p: 2 }}>
                                <Typography variant="subtitle2" sx={{ fontWeight: 600, mb: 1 }}>
                                  {fail.reactionType || "Unknown Reaction"}
                                </Typography>
                                <Typography variant="body2" sx={{ mb: 1, color: "#5d6d7e" }}>
                                  Reactants: {fail.reactants.join(", ")}
                                </Typography>
                                <Typography variant="body2" sx={{ color: "#d32f2f" }}>
                                  Reason: {fail.reason}
                                </Typography>
                              </CardContent>
                            </Card>
                          ))
                        )}
                      </Box>

                      <Card variant="outlined" sx={{ borderRadius: 1, bgcolor: "#fafafa" }}>
                        <CardHeader
                          title={
                            <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                              Statistics
                            </Typography>
                          }
                        />
                        <CardContent sx={{ pt: 0 }}>
                          <Grid container spacing={2}>
                            <Grid item xs={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                Mean MW:
                              </Typography>
                              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                                {reactionResult.statistics.mean_mw.toFixed(2)}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                Std MW:
                              </Typography>
                              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                                {reactionResult.statistics.std_mw.toFixed(2)}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                Min MW:
                              </Typography>
                              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                                {reactionResult.statistics.min_mw.toFixed(2)}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, color: "#5d6d7e" }}>
                                Max MW:
                              </Typography>
                              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                                {reactionResult.statistics.max_mw.toFixed(2)}
                              </Typography>
                            </Grid>
                          </Grid>
                        </CardContent>
                      </Card>
                    </Grid>
                  </Grid>
                </CardContent>
              </Paper>
            )}
          </Box>
        )}
      </Container>
    </Box>
  )
}

export default RDkit