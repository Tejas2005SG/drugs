import CostEstimation from '../models/costestimination.model.js';
import axios from 'axios';

const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

// Function to call Gemini AI for cost estimation
// Function to call Gemini AI for cost estimation
// Function to call Gemini AI for cost estimation
const geminiCostEstimation = async (smiles) => {
    try {
      const prompt = `
        Analyze the drug molecule represented by the SMILES string "${smiles}" and estimate its synthesis and production cost in USD per gram. Follow the steps below to ensure a precise and realistic cost estimation:

### 1. Structural Analysis
- Identify the molecular structure, including:
  - Number of chiral centers and their impact on stereoselective synthesis.
  - Types and number of functional groups (e.g., amines, esters, halogens).
  - Presence of fused ring systems, heterocycles, or macrocycles.
- Assess how these structural features affect synthetic difficulty:
  - Chiral centers may require expensive chiral catalysts or resolution techniques.
  - Complex ring systems may need multi-step cyclization reactions.
  - Functional groups may require specific reagents or conditions (e.g., high-pressure reactions for certain reductions).

### 2. Synthetic Route Analysis
- Perform a detailed retrosynthetic analysis to propose a feasible synthetic pathway:
  - Break down the molecule into simpler precursors.
  - Estimate the number of synthetic steps (e.g., 3–5 steps for simple molecules, 10+ steps for complex ones).
  - Identify key reagents, catalysts, and reaction conditions (e.g., temperature, pressure, inert atmosphere).
- Estimate the cost impact of the synthetic route:
  - Consider the cost of reagents (e.g., palladium catalysts can be expensive).
  - Account for specialized conditions (e.g., cryogenic temperatures or high-pressure reactors).
  - Evaluate the yield of each step and its effect on overall cost (e.g., low-yield steps increase material costs).

### 3. Commercial Availability and Market Data
- Check the availability of the molecule or similar compounds in chemical databases:
  - Query databases like ChemSpace, ZINC, or Sigma-Aldrich for pricing data.
  - If the exact molecule is unavailable, identify structurally similar molecules (e.g., same core scaffold or functional groups) and use their pricing as a benchmark.
- Consider market factors:
  - Assess demand for the molecule (e.g., high demand for pharmaceuticals may increase cost).
  - Evaluate supply chain factors (e.g., scarcity of starting materials or intermediates).
  - Account for regional pricing differences (e.g., costs may vary between the US, Europe, and Asia).

### 4. Cost Factors
- Analyze additional factors that influence the cost:
  - **Starting Materials**: Are the precursors commercially available and affordable, or do they need to be synthesized?
  - **Scale of Production**: Lab-scale synthesis (mg quantities) is more expensive per gram than industrial-scale (kg quantities).
  - **Purity Requirements**: Higher purity (e.g., >99% for pharmaceuticals) requires additional purification steps, increasing costs.
  - **Regulatory and Safety Considerations**: Handling hazardous reagents (e.g., toxic or flammable compounds) may require specialized equipment or safety protocols.
  - **Labor and Expertise**: Complex syntheses may require skilled chemists or specialized facilities, increasing labor costs.
- Provide a breakdown of how each factor contributes to the final cost.

### 5. Benchmarking with Known Molecules
- Compare the molecule to well-known compounds with similar complexity to guide the cost estimate:
  - For example, simple molecules like aspirin (SMILES: CC(=O)OC1=CC=CC=C1C(=O)O) typically cost $5–$20 per gram.
  - Moderately complex molecules like caffeine (SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C) may cost $20–$50 per gram.
  - Complex molecules like paclitaxel (a cancer drug) can cost $500–$2000 per gram due to their synthetic complexity.
- Use these benchmarks to ensure the estimated cost is realistic and aligns with industry standards.

### Output Requirements
- Provide an estimated price range in the format: "Estimated Cost: $X–$Y per gram" (e.g., "Estimated Cost: $30–$80 per gram").
  - Ensure the range is as narrow as possible (e.g., within a factor of 2–3, not 10) based on the analysis.
- Follow this with detailed reasoning for the estimation in the following format:
  - Start with a concise introductory paragraph (1–2 sentences) about the SMILES string (e.g., "The SMILES string '...' represents a molecule with...") without bullet points.
  - Use numbered points for main sections (e.g., "1. Structural Analysis").
  - Use bullet points for sub-details under each section (e.g., "- The molecule has two chiral centers, increasing synthesis costs.").
  - Ensure each bullet point is concise (1 sentence) and avoids redundancy.
  - Do not use bullet points for the introductory paragraph or between main sections.
  - Do not use stars (*) or other special characters (e.g., **, ---) in the description.
  - Use proper grammar and ensure clarity in each point.
- Highlight any assumptions made during the estimation (e.g., "Assuming lab-scale production due to lack of market data.").

### Example Output
For the SMILES string "CC(=O)OC1=CC=CC=C1C(=O)O" (aspirin):
- Estimated Cost: $5–$15 per gram
- The SMILES string "CC(=O)OC1=CC=CC=C1C(=O)O" represents aspirin, a simple aromatic molecule with an ester and carboxylic acid group.
- 1. Structural Analysis:
  - The molecule has no chiral centers, simplifying synthesis.
  - It contains an ester and carboxylic acid, which are common functional groups.
  - The aromatic ring is straightforward to synthesize.
- 2. Synthetic Route:
  - Retrosynthesis suggests starting from salicylic acid and acetylating it with acetic anhydride.
  - The synthesis requires 1–2 steps with high yields.
  - Reagents like acetic anhydride are inexpensive and widely available.
- 3. Commercial Availability:
  - Aspirin is widely available, with prices around $5–$15 per gram in chemical catalogs.
  - Similar aromatic esters have comparable pricing.
- 4. Cost Factors:
  - Starting materials (salicylic acid) are cheap and abundant.
  - No specialized equipment is required.
  - Industrial-scale production further reduces costs.
- Assumptions: The estimate assumes small-scale production for research purposes.
      `;
  
      const response = await axios.post(GEMINI_API_URL, {
        contents: [{ role: "user", parts: [{ text: prompt }] }],
      }, {
        headers: { "Content-Type": "application/json", "x-goog-api-key": GEMINI_API_KEY },
      });
  
      let generatedText = response.data?.candidates?.[0]?.content?.parts?.[0]?.text || "No response from AI";
  
      // Parse the response to extract the cost
      let estimatedCost = "N/A";
      let information = generatedText;
  
      // Look for a line like "Estimated Cost: $X–$Y per gram"
      const costRegex = /Estimated Cost: \$(\d+)(?:–\$(\d+))? per gram/i;
      const costMatch = generatedText.match(costRegex);
  
      if (costMatch) {
        const lowerBound = costMatch[1];
        const upperBound = costMatch[2] || lowerBound; // If no range, use the single value
        estimatedCost = `$${lowerBound}–$${upperBound} per gram`;
        // Remove the "Estimated Cost" line from the information
        information = generatedText.replace(costRegex, "").trim();
      } else {
        // Fallback: Look for any price range in the text (e.g., "$50–$100 per gram")
        const fallbackRegex = /\$(\d+)(?:–\$(\d+))? per gram/i;
        const fallbackMatch = generatedText.match(fallbackRegex);
        if (fallbackMatch) {
          const lowerBound = fallbackMatch[1];
          const upperBound = fallbackMatch[2] || lowerBound;
          estimatedCost = `$${lowerBound}–$${upperBound} per gram`;
        }
      }
  
      // Clean up the information: Remove stars and ensure proper formatting
      information = information
        .replace(/\*+/g, "") // Remove all stars
        .replace(/(\d+\.\s[A-Za-z\s]+:)/g, "\n$1") // Ensure numbered sections are on new lines
        .replace(/-+/g, "-") // Normalize dashes for bullet points
        .trim();
  
      return {
        estimatedcost: estimatedCost,
        information: information || "No additional information provided",
      };
    } catch (error) {
      console.error("Error calling Gemini AI:", error.message);
      return { estimatedcost: "Error estimating cost", information: "AI service unavailable" };
    }
  };


// Post cost estimation using Gemini AI
export const postCostEstimation = async (req, res) => {
  try {
    const { smiles } = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }
    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }

    const { estimatedcost, information } = await geminiCostEstimation(smiles);

    const costEstimation = new CostEstimation({
      smiles,
      estimatedcost,
      information,
      userId,
    });

    await costEstimation.save();

    res.status(201).json({
      success: true,
      data: costEstimation,
    });
  } catch (error) {
    res.status(500).json({
      success: false,
      message: 'Server error',
      error: error.message,
    });
  }
};

// Get cost estimations for a user
export const getCostEstimation = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: 'User not authenticated' });
    }

    const estimations = await CostEstimation.find({ userId }).sort({ created: -1 });

    if (!estimations.length) {
      return res.status(404).json({ message: 'No cost estimations found for this user' });
    }

    res.status(200).json({
      success: true,
      data: estimations,
    });
  } catch (error) {
    res.status(500).json({
      success: false,
      message: 'Server error',
      error: error.message,
    });
  }
};
