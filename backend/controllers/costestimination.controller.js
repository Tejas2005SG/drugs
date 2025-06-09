import CostEstimation from "../models/costestimination.model.js";
import axios from "axios";

import { MongoClient } from "mongodb";

const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL =
  "https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-05-20:generateContent";

// Enhanced cost estimation algorithm with multiple factors
const calculateMolecularComplexity = (smiles) => {
  if (!smiles || typeof smiles !== 'string') return 1.0;

  let complexity = 1.0;
  
  // Base complexity from molecular size
  const length = smiles.length;
  if (length > 100) complexity *= 1.5;
  else if (length > 50) complexity *= 1.2;
  
  // Functional group complexity factors
  const functionalGroups = {
    // Halides (expensive but common)
    'F': 1.1, 'Cl': 1.1, 'Br': 1.3, 'I': 1.8,
    // Heteroatoms
    'N': 1.1, 'O': 1.05, 'S': 1.2, 'P': 1.4,
    // Ring systems
    'c': 1.1, // aromatic carbon
    // Stereochemistry
    '@': 1.3, '@@': 1.3,
    // Multiple bonds
    '=': 1.1, '#': 1.2,
    // Branching - need to escape special regex characters
    '\\(': 1.05, '\\)': 1.05
  };
  
  for (const [group, factor] of Object.entries(functionalGroups)) {
    const count = (smiles.match(new RegExp(group, 'g')) || []).length;
    complexity *= Math.pow(factor, count * 0.1); // Damped effect
  }
  
  // Cap complexity to prevent extreme values
  return Math.min(complexity, 3.0);
};

const estimateStepCount = (smiles) => {
  if (!smiles || typeof smiles !== 'string') return 3;
  
  const length = smiles.length;
  let steps = 3; // Base minimum steps
  
  // Estimate steps based on molecular complexity
  if (length > 80) steps += 4;
  else if (length > 60) steps += 3;
  else if (length > 40) steps += 2;
  else if (length > 20) steps += 1;
  
  // Additional steps for complex functional groups
  const complexGroups = ['I', 'Br', '@', '@@', '#'];
  for (const group of complexGroups) {
    steps += (smiles.match(new RegExp(group, 'g')) || []).length * 0.5;
  }
  
  return Math.min(Math.ceil(steps), 12); // Cap at 12 steps
};

const calculateBaseCost = (complexity, stepCount) => {
  // Base cost factors (realistic pharmaceutical synthesis costs)
  const baseCostPerStep = 15; // $15 per gram per synthesis step
  const complexityMultiplier = 1 + (complexity - 1) * 0.5; // Reduce impact
  
  const rawCost = baseCostPerStep * stepCount * complexityMultiplier;
  
  // Apply realistic cost caps and floors
  const minCost = 5; // $5/gram minimum
  const maxCost = 300; // $300/gram maximum for most compounds
  
  return Math.max(minCost, Math.min(rawCost, maxCost));
};

const createCostRange = (baseCost) => {
  // Create a reasonable range around the base cost
  const variability = 0.3; // ±30% variability
  const lowerBound = Math.max(1, Math.floor(baseCost * (1 - variability)));
  const upperBound = Math.ceil(baseCost * (1 + variability));
  
  return { lowerBound, upperBound };
};

// Enhanced Gemini cost estimation with algorithm-based validation
const geminiCostEstimation = async (smiles) => {
  try {
    // Calculate algorithmic estimate first
    const complexity = calculateMolecularComplexity(smiles);
    const stepCount = estimateStepCount(smiles);
    const algorithmicCost = calculateBaseCost(complexity, stepCount);
    const costRange = createCostRange(algorithmicCost);
    
    const prompt = `
You are a pharmaceutical cost estimation AI. Analyze the SMILES string "${smiles}" and provide a realistic cost estimate.

IMPORTANT CONSTRAINTS:
- Most pharmaceutical compounds cost between $5-$300 per gram for synthesis
- Avoid extreme estimates above $300/gram unless truly exceptional
- Consider that the algorithmic estimate suggests $${algorithmicCost.toFixed(0)} per gram
- Use this as a baseline and adjust based on specific structural challenges

**Required Output Format:**
**Estimated Cost:** \`$${costRange.lowerBound}–$${costRange.upperBound} per gram\`
*(Algorithmic baseline: $${algorithmicCost.toFixed(0)}/gram, Complexity: ${complexity.toFixed(2)}, Steps: ${stepCount})*

**Molecular Analysis:**
- **IUPAC Name:** [From PubChem or estimate]
- **Molecular Weight:** [Calculate from SMILES]
- **Key Functional Groups:** [List major groups]
- **Synthesis Complexity:** [Rate 1-5 scale]

**Cost Breakdown:**
- **Raw Materials:** $X-Y per gram
- **Labor & Equipment:** $X-Y per gram  
- **Purification:** $X-Y per gram
- **Quality Control:** $X-Y per gram

**Synthesis Route:**
- **Primary Route:** [Describe likely synthetic approach]
- **Critical Steps:** [2-3 key challenging steps]
- **Expected Yield:** X% overall

**Market Comparison:**
- **Similar Compounds:** [Price references if available]
- **Commercial Availability:** [Check if precursors are available]
- **Scale Considerations:** [Bulk pricing potential]

**Cost Drivers & Risks:**
- **Major Cost Factors:** [Top 2-3 cost drivers]
- **Technical Risks:** [Synthesis challenges]
- **Regulatory Considerations:** [If any]

VALIDATION RULES:
1. Final cost MUST be reasonable ($5-$300/gram range)
2. Justify any costs above $100/gram with specific technical challenges
3. Consider that most APIs are produced at $10-$80/gram
4. Account for economies of scale in larger batches
`;

    const response = await axios.post(
      GEMINI_API_URL,
      {
        contents: [{ role: "user", parts: [{ text: prompt }] }],
      },
      {
        headers: {
          "Content-Type": "application/json",
          "x-goog-api-key": GEMINI_API_KEY,
        },
      }
    );

    let generatedText = response.data?.candidates?.[0]?.content?.parts?.[0]?.text || 
                       "No response from AI";

    // Enhanced cost extraction with validation
    let estimatedCost = `$${costRange.lowerBound}–$${costRange.upperBound} per gram`;
    let information = generatedText;

    // Try to extract cost from AI response
    const costRegex = /Estimated Cost: \$(\d+)(?:–\$(\d+))? per gram/i;
    const costMatch = generatedText.match(costRegex);

    if (costMatch) {
      const aiLowerBound = parseInt(costMatch[1]);
      const aiUpperBound = parseInt(costMatch[2] || costMatch[1]);
      
      // Validate AI estimate against reasonable bounds
      if (aiLowerBound <= 500 && aiUpperBound <= 500) {
        // Use AI estimate if reasonable
        estimatedCost = `$${aiLowerBound}–$${aiUpperBound} per gram`;
      } else {
        // Fall back to algorithmic estimate if AI estimate is too high
        console.warn(`AI estimate too high (${aiLowerBound}-${aiUpperBound}), using algorithmic estimate`);
        estimatedCost = `$${costRange.lowerBound}–$${costRange.upperBound} per gram`;
      }
    }

    // Clean up the information
    information = information
      .replace(/\*+/g, "")
      .replace(/(\d+\.\s[A-Za-z\s]+:)/g, "\n$1")
      .replace(/-+/g, "-")
      .trim();

    // Add algorithmic insights
    const algorithmicInsights = `
    
**Algorithmic Analysis:**
- Molecular Complexity Score: ${complexity.toFixed(2)}/3.0
- Estimated Synthesis Steps: ${stepCount}
- Base Cost Calculation: $${algorithmicCost.toFixed(0)} per gram
- Cost Range: $${costRange.lowerBound}–$${costRange.upperBound} per gram
    `;

    return {
      estimatedcost: estimatedCost,
      information: information + algorithmicInsights,
      algorithmicData: {
        complexity,
        stepCount,
        baseCost: algorithmicCost,
        costRange
      }
    };

  } catch (error) {
    console.error("Error calling Gemini AI:", error.message);
    
    // Fallback to pure algorithmic estimation
    const complexity = calculateMolecularComplexity(smiles);
    const stepCount = estimateStepCount(smiles);
    const algorithmicCost = calculateBaseCost(complexity, stepCount);
    const costRange = createCostRange(algorithmicCost);
    
    return {
      estimatedcost: `$${costRange.lowerBound}–$${costRange.upperBound} per gram`,
      information: `Algorithmic estimation (AI service unavailable):
      - Molecular Complexity: ${complexity.toFixed(2)}/3.0
      - Estimated Steps: ${stepCount}
      - Base Cost: $${algorithmicCost.toFixed(0)} per gram
      - This estimate is based on molecular complexity analysis and typical pharmaceutical synthesis costs.`,
      algorithmicData: {
        complexity,
        stepCount,
        baseCost: algorithmicCost,
        costRange
      }
    };
  }
};

// Updated post cost estimation function
export const postCostEstimation = async (req, res) => {
  try {
    const { smiles } = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: "User not authenticated" });
    }
    if (!smiles) {
      return res.status(400).json({ message: "SMILES string is required" });
    }

    // Validate SMILES string basic format
    if (typeof smiles !== 'string' || smiles.length < 2) {
      return res.status(400).json({ message: "Invalid SMILES string format" });
    }

    const { estimatedcost, information, algorithmicData } = await geminiCostEstimation(smiles);

    const costEstimation = new CostEstimation({
      smiles,
      estimatedcost,
      information,
      userId,
      algorithmicData, // Store algorithmic data for future analysis
    });

    await costEstimation.save();

    res.status(201).json({
      success: true,
      data: costEstimation,
      message: "Cost estimation completed using hybrid AI-algorithmic approach"
    });
  } catch (error) {
    res.status(500).json({
      success: false,
      message: "Server error",
      error: error.message,
    });
  }
};

// Get cost estimations for a user (unchanged)
export const getCostEstimation = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: "User not authenticated" });
    }

    const estimations = await CostEstimation.find({ userId }).sort({
      created: -1,
    });

    if (!estimations.length) {
      return res
        .status(404)
        .json({ message: "No cost estimations found for this user" });
    }

    res.status(200).json({
      success: true,
      data: estimations,
    });
  } catch (error) {
    res.status(500).json({
      success: false,
      message: "Server error",
      error: error.message,
    });
  }
};

// Rest of the code remains the same...
