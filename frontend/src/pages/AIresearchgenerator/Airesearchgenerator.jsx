import React, { useState, useEffect ,useRef} from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import jsPDF from "jspdf";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const CROSSREF_API_URL = "https://api.crossref.org/works";

// Axios instance for backend requests (with credentials)
const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

// Axios instance for CrossRef requests (no credentials)
const crossrefAxios = axios.create({
  baseURL: CROSSREF_API_URL,
  withCredentials: false,
});

// Utility to clean HTML tags from abstracts
const cleanAbstract = (abstract) => {
  if (!abstract || abstract === "Abstract not available.") return abstract;
  return abstract
    .replace(/<\/?[^>]+(>|$)/g, "") // Remove HTML tags
    .replace(/\s+/g, " ") // Normalize whitespace
    .trim();
};

function Airesearchgenerator() {
  const [activeTab, setActiveTab] = useState("related");
  const [molecules, setMolecules] = useState([]);
  const [selectedTitle, setSelectedTitle] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [researchPapers, setResearchPapers] = useState([]); // Newly fetched papers
  const [savedPapers, setSavedPapers] = useState([]); // Previously saved related papers
  const [savedGeneratedPapers, setSavedGeneratedPapers] = useState([]); // Previously saved generated papers
  const [researchSummary, setResearchSummary] = useState("");
  const [generatedPaper, setGeneratedPaper] = useState(null); // Store the generated research paper
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
 const toastShown = useRef(false);
  const { user, checkAuth, checkingAuth } = useAuthStore();

      useEffect(() => {
      if (!toastShown.current) {
        toast(
          "First , Generate a new molecule ",
          {
            position: "top-right",
            duration: 8000, // 6 seconds
            style: {
              background: "#fefcbf", // Yellow background for warning
              color: "#92400e", // Dark yellow text
              border: "1px solid #f59e0b",
            },
            icon: "⚠️", // Warning icon
          }
        );
        toastShown.current = true; // Mark the toast as shown
      }
    }, []); 

  useEffect(() => {
    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError("Authentication failed. Please log in.");
        return;
      }
      await fetchAllMolecules();
      await fetchSavedPapers();
      await fetchSavedGeneratedPapers();
    };
    initializeApp();
  }, [checkAuth]);

  const fetchAllMolecules = async () => {
    if (!user?._id) return;

    setLoading(true);
    try {
      const response = await axiosInstance.get("/protein/generatednewmolecule");
      const fetchedMolecules = response.data.molecules || [];
      setMolecules(fetchedMolecules);
      if (fetchedMolecules.length > 0 && !selectedTitle && !selectedSmiles) {
        setSelectedTitle(fetchedMolecules[0].newmoleculetitle);
        setSelectedSmiles(fetchedMolecules[0].newSmiles);
      }
    } catch (err) {
      console.error("Error fetching molecules:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch molecules");
    } finally {
      setLoading(false);
    }
  };

  const fetchSavedPapers = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/protein/saved-research-papers");
      setSavedPapers(response.data.papers || []);
    } catch (err) {
      console.error("Error fetching saved papers:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch saved papers");
    }
  };

  const fetchSavedGeneratedPapers = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/protein/saved-generated-research-papers");
      setSavedGeneratedPapers(response.data.papers || []);
    } catch (err) {
      console.error("Error fetching saved generated papers:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch saved generated papers");
    }
  };

  const checkIfPapersExist = async (title, smiles) => {
    if (!user?._id || !title || !smiles) return false;

    try {
      const response = await axiosInstance.get("/protein/check-saved-papers", {
        params: { title, smiles },
      });
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved papers:", err.response?.data || err.message);
      return false;
    }
  };

  const checkIfGeneratedPaperExists = async (title, smiles) => {
    if (!user?._id || !title || !smiles) return false;

    try {
      const response = await axiosInstance.get("/protein/check-saved-generated-papers", {
        params: { title, smiles },
      });
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved generated papers:", err.response?.data || err.message);
      return false;
    }
  };

  const savePapers = async (molecule, papers) => {
    if (!user?._id || !papers.length) return;

    const payload = {
      userId: user._id,
      molecule: { title: molecule.newmoleculetitle, smiles: molecule.newSmiles },
      papers,
    };
    console.log("Saving papers with payload:", payload);

    try {
      await axiosInstance.post("/protein/save-research-papers", payload);
    } catch (err) {
      console.error("Error saving papers:", err.response?.data || err.message);
      toast.error("Failed to save research papers");
    }
  };

  const saveGeneratedPaper = async (molecule, paper) => {
    if (!user?._id || !paper) return;

    const payload = {
      userId: user._id,
      molecule: { title: molecule.newmoleculetitle, smiles: molecule.newSmiles },
      paper,
    };
    console.log("Saving generated paper with payload:", payload);

    try {
      await axiosInstance.post("/protein/save-generated-research-paper", payload);
      await fetchSavedGeneratedPapers(); // Refresh the saved generated papers list
    } catch (err) {
      console.error("Error saving generated paper:", err.response?.data || err.message);
      toast.error("Failed to save generated research paper");
    }
  };

  
  // const handleResearchClick = async () => {
  //   if (!selectedTitle || !selectedSmiles) {
  //     toast.error("Please select both a title and SMILES string");
  //     return;
  //   }
  
  //   const papersExist = await checkIfPapersExist(selectedTitle, selectedSmiles);
  //   if (papersExist) {
  //     toast("Research papers for this molecule are already saved. Redirecting to Saved Research Papers.", {
  //       type: "info",
  //     });
  //     setActiveTab("saved");
  //     return;
  //   }
  
  //   setLoading(true);
  //   setError(null);
  //   setResearchPapers([]);
  //   setResearchSummary("");
  
  //   const selectedMol = molecules.find(
  //     (mol) => mol.newmoleculetitle === selectedTitle && mol.newSmiles === selectedSmiles
  //   );
  //   if (!selectedMol) {
  //     setError("Selected molecule not found");
  //     setLoading(false);
  //     return;
  //   }
  
  //   const uniqueId = `${Date.now()}-${Math.random().toString(36).substring(2, 9)}`; // Prevent caching
  
  //   const geminiPrompt = `
  //     You are a world-class expert in cheminformatics, molecular modeling, and scientific literature synthesis. The molecule provided is newly generated, with no exact research papers existing. Your task is to infer and generate **unique research paper metadata** in IEEE format, strictly based on the SMILES string's chemical structure. You will assist in querying the IEEE Xplore API (using the provided API key) to fetch related research papers based on structural analogs or chemical properties derived from the SMILES.
  
  //     ### Molecule Details:
  //     - **SMILES String**: "${selectedMol.newSmiles}"
  //     - **IUPAC Name**: "${selectedMol.newIupacName}" (use only for validation, not primary input)
  //     - **Conversion Details (Synthesis Pathway)**: "${selectedMol.conversionDetails}"
  //     - **Potential Diseases (Therapeutic Relevance)**: "${selectedMol.potentialDiseases}"
  //     - **Additional Information**: "${selectedMol.information}"
  //     - **Unique Request ID**: "${uniqueId}" (ensure uniqueness)
  
  //     ### Deep Training Requirements:
  //     1. **SMILES Structural Analysis**:
  //        - Perform a detailed breakdown of the SMILES string "${selectedMol.newSmiles}" to identify:
  //          - Functional groups (e.g., hydroxyl, carbonyl, amine)
  //          - Ring systems (e.g., aromatic, aliphatic)
  //          - Substituents and stereochemistry
  //          - Potential reactivity or bioactivity based on these features.
  //        - Use this analysis as the sole basis for generating paper titles, abstracts, and relevance.
  
  //     2. **Generate Unique Research Paper Metadata**:
  //        - Provide details for 2-3 hypothetical research papers derived from the SMILES structure.
  //        - Ensure unique authors (3 per paper), publication years (2015-2025), and IEEE-compliant DOIs.
  //        - Abstracts (50-100 words) must reflect SMILES-derived structure and inferred applications (e.g., drug design, synthesis).
  //        - Mark these as simulated if no API data is available.
  
  //     3. **IEEE Xplore API Integration**:
  //        - Use the API key [INSERT_YOUR_API_KEY_HERE] to query the IEEE Xplore API.
  //        - Construct a query based on SMILES-derived keywords (e.g., functional groups, chemical classes).
  //        - Example API endpoint: https://ieeexploreapi.ieee.org/api/v1/search/articles?apikey=[INSERT_YOUR_API_KEY_HERE]&query=SMILES_KEYWORDS.
  
  //     4. **Output Format**:
  //        Return the response in a structured JSON format **and nothing else**:
  //        {
  //          "summary": "A unique 2-3 sentence explanation connecting the SMILES structure to research areas.",
  //          "generated_papers": [
  //            {
  //              "title": "Unique Paper Title 1 Based on SMILES",
  //              "authors": "A. Kim, B. Patel, C. Nguyen",
  //              "year": "2021",
  //              "abstract": "Abstract summarizing SMILES-derived features and applications.",
  //              "doi": "10.1109/ACCESS.2021.7654321",
  //              "url": "https://ieeexplore.ieee.org/document/7654321",
  //              "is_simulated": true
  //            }
  //          ],
  //          "ieee_api_papers": [
  //            {
  //              "title": "Real Paper Title from API",
  //              "authors": "Author List from API",
  //              "year": "Year from API",
  //              "abstract": "Abstract from API",
  //              "doi": "DOI from API",
  //              "url": "URL from API",
  //              "is_simulated": false
  //            }
  //          ]
  //        }
  //   `;
  
  //   try {
  //     const geminiResponse = await axiosInstance.post("/protein/proxy/gemini", {
  //       prompt: geminiPrompt,
  //     });
  //     const geminiContent = geminiResponse.data.content;
  //     console.log("Gemini raw response:", geminiContent); // Debug log
  
  //     const jsonMatch = geminiContent.match(/{[\s\S]*}/);
  //     if (!jsonMatch) {
  //       throw new Error("No valid JSON found in Gemini response");
  //     }
  
  //     const parsedGemini = JSON.parse(jsonMatch[0]);
  //     setResearchSummary(parsedGemini.summary || "");
  
  //     // Process generated papers
  //     let papers = parsedGemini.generated_papers || [];
  //     papers = papers.map((paper) => ({ ...paper, is_simulated: true }));
  
  //     // Fetch from IEEE Xplore API if API key is provided
  //     const apiKey = "[INSERT_YOUR_API_KEY_HERE]"; // Replace with your actual API key
  //     if (apiKey && apiKey !== "[INSERT_YOUR_API_KEY_HERE]") {
  //       // Derive keywords from SMILES for API query
  //       const smilesKeywords = selectedMol.newSmiles
  //         .replace(/[^A-Za-z0-9]/g, " ")
  //         .split(" ")
  //         .filter((kw) => kw.length > 2)
  //         .join(" OR ");
  //       const ieeeApiUrl = `https://ieeexploreapi.ieee.org/api/v1/search/articles?apikey=${apiKey}&query=${encodeURIComponent(smilesKeywords)}&max_records=3&sort_order=relevance`;
  
  //       try {
  //         const ieeeResponse = await axios.get(ieeeApiUrl);
  //         const ieeePapers = (ieeeResponse.data.articles || []).map((article) => ({
  //           title: article.title || "Untitled",
  //           authors: article.authors?.map((a) => a.full_name).join(", ") || "Unknown Authors",
  //           year: article.publication_year || "Unknown Year",
  //           abstract: article.abstract || "Abstract not available",
  //           doi: article.doi || "No DOI available",
  //           url: article.url || (article.doi ? `https://ieeexplore.ieee.org/document/${article.doi.split("/").pop()}` : "No URL available"),
  //           is_simulated: false,
  //         }));
  
  //         papers = [...papers, ...ieeePapers].slice(0, 3); // Combine and limit to 3 papers
  //         console.log("IEEE API papers fetched:", ieeePapers); // Debug log
  //       } catch (apiErr) {
  //         console.error("IEEE API error:", apiErr.message);
  //         toast.error("Failed to fetch IEEE papers; using simulated data.");
  //       }
  //     } else {
  //       console.warn("No valid API key provided. Using only Gemini-generated papers.");
  //     }
  
  //     // Validate and adjust URLs
  //     papers = papers.map((paper, index) => {
  //       let validUrl = paper.url;
  //       if (!validUrl || !validUrl.startsWith("http")) {
  //         const ieeeDoi = paper.doi || `10.1109/ACCESS.${Date.now() + index}${Math.floor(Math.random() * 1000)}`;
  //         validUrl = `https://ieeexplore.ieee.org/document/${ieeeDoi.split("/").pop()}`;
  //       }
  //       // Test URL (optional, requires additional logic or external service)
  //       return {
  //         ...paper,
  //         doi: paper.doi || validUrl.split("/document/")[1],
  //         url: validUrl,
  //         is_valid: paper.is_simulated ? false : true, // Assume API URLs are valid unless proven otherwise
  //       };
  //     });
  
  //     setResearchPapers(papers);
  //     await savePapers(selectedMol, papers);
  //     toast.success("Related research papers fetched and saved successfully!");
  //   } catch (err) {
  //     console.error("Error fetching research:", err);
  //     setError(err.response?.data?.message || "Failed to fetch research papers");
  //     toast.error("Failed to fetch research papers");
  //   } finally {
  //     setLoading(false);
  //   }
  // };

  const handleResearchClick = async () => {
    if (!selectedTitle || !selectedSmiles) {
      toast.error("Please select both a title and SMILES string");
      return;
    }
  
    const papersExist = await checkIfPapersExist(selectedTitle, selectedSmiles);
    if (papersExist) {
      toast("Research papers for this molecule are already saved. Redirecting to Saved Research Papers.", {
        type: "info",
      });
      setActiveTab("saved");
      return;
    }
  
    setLoading(true);
    setError(null);
    setResearchPapers([]);
    setResearchSummary("");
  
    const selectedMol = molecules.find(
      (mol) => mol.newmoleculetitle === selectedTitle && mol.newSmiles === selectedSmiles
    );
    if (!selectedMol) {
      setError("Selected molecule not found");
      setLoading(false);
      return;
    }
  
    const uniqueId = `${Date.now()}-${Math.random().toString(36).substring(2, 9)}`; // Prevent caching
  
    const geminiPrompt = `
      You are a world-class expert in cheminformatics, molecular modeling, and scientific literature synthesis with deep knowledge of IEEE Xplore conventions. The molecule provided is newly generated, with no exact research papers existing. Your task is to:
      1. Perform a detailed SMILES analysis to infer chemical properties and generate 2-3 unique, plausible research paper metadata entries in IEEE format, based solely on the SMILES structure.
      2. Assist in constructing a query for the IEEE Xplore API to fetch up to 3 real research papers related to structural analogs or chemical properties derived from the SMILES.
      3. Ensure all generated and fetched metadata includes valid, IEEE-compliant DOIs and URLs that align with https://ieeexplore.ieee.org/document/ format.
  
      ### Molecule Details:
      - **SMILES String**: "${selectedMol.newSmiles}"
      - **IUPAC Name**: "${selectedMol.newIupacName || 'Not Available'}" (for validation only)
      - **Conversion Details**: "${selectedMol.conversionDetails || 'Not Available'}"
      - **Potential Diseases**: "${selectedMol.potentialDiseases || 'Not Available'}"
      - **Additional Information**: "${selectedMol.information || 'Not Available'}"
      - **Unique Request ID**: "${uniqueId}"
  
      ### Detailed Instructions:
      1. **SMILES Structural Analysis**:
         - Break down "${selectedMol.newSmiles}" into:
           - Functional groups (e.g., -OH, -COOH, -NH2)
           - Ring systems (e.g., benzene, pyridine)
           - Substituents and stereochemistry
           - Inferred bioactivity (e.g., anti-inflammatory, antimicrobial).
         - Use this to craft unique paper titles, abstracts (50-100 words), and relevance.
  
      2. **Generate Unique Research Paper Metadata**:
         - Create 2-3 hypothetical papers with:
           - Unique authors (e.g., "J. Smith, L. Chen, M. Garcia")
           - Publication years (2015-2025)
           - IEEE-compliant DOIs (e.g., 10.1109/TBME.2022.9876543)
           - Abstracts reflecting SMILES-derived features and applications.
         - Assign URLs as https://ieeexplore.ieee.org/document/{DOI_suffix} using the DOI’s numeric part.
         - Mark these as "is_simulated": true.
  
      3. **IEEE Xplore API Query Construction**:
         - Derive keywords from SMILES (e.g., functional groups, ring types) and combine with terms like "drug design," "cheminformatics," or "molecular modeling."
         - Simulate an API query structure: https://ieeexploreapi.ieee.org/api/v1/search/articles?apikey=[INSERT_YOUR_API_KEY_HERE]&query=KEYWORDS&max_records=3.
         - If API access is unavailable, suggest fallback keywords for manual search.
  
      4. **Output Format**:
         Return the response in a structured JSON format **and nothing else**:
         {
           "summary": "A 2-3 sentence explanation linking the SMILES structure to research areas, noting API or simulation status.",
           "generated_papers": [
             {
               "title": "Unique Paper Title Based on SMILES Analysis",
               "authors": "J. Smith, L. Chen, M. Garcia",
               "year": "2023",
               "abstract": "Abstract detailing SMILES-derived structure and inferred applications (50-100 words).",
               "doi": "10.1109/TBME.2023.1234567",
               "url": "https://ieeexplore.ieee.org/document/1234567",
               "is_simulated": true
             }
           ],
           "ieee_api_papers": [
             {
               "title": "Real Paper Title from API",
               "authors": "Author List from API",
               "year": "Year from API",
               "abstract": "Abstract from API",
               "doi": "DOI from API",
               "url": "https://ieeexplore.ieee.org/document/{DOI_suffix_from_API}",
               "is_simulated": false
             }
           ]
         }
    `;
  
    try {
      const geminiResponse = await axiosInstance.post("/protein/proxy/gemini", {
        prompt: geminiPrompt,
      });
      const geminiContent = geminiResponse.data.content;
      console.log("Gemini raw response:", geminiContent); // Debug log
  
      const jsonMatch = geminiContent.match(/{[\s\S]*}/);
      if (!jsonMatch) {
        throw new Error("No valid JSON found in Gemini response");
      }
  
      const parsedGemini = JSON.parse(jsonMatch[0]);
      setResearchSummary(parsedGemini.summary || "No summary available.");
  
      // Process generated papers with workable URLs
      let papers = (parsedGemini.generated_papers || []).map((paper) => ({
        ...paper,
        is_simulated: true,
        url: constructIeeeUrl(paper.doi || generateUniqueDoi(paper.year || 2025)),
      }));
  
      // Attempt IEEE Xplore API fetch
      const apiKey = "[INSERT_YOUR_API_KEY_HERE]"; // Replace with your actual API key
      if (apiKey && apiKey !== "[INSERT_YOUR_API_KEY_HERE]") {
        const smilesKeywords = selectedMol.newSmiles
          .replace(/[^A-Za-z0-9]/g, " ")
          .split(" ")
          .filter((kw) => kw.length > 2)
          .join(" OR ");
        const ieeeApiUrl = `https://ieeexploreapi.ieee.org/api/v1/search/articles?apikey=${apiKey}&query=${encodeURIComponent(smilesKeywords)}&max_records=3&sort_order=relevance`;
  
        try {
          const ieeeResponse = await axios.get(ieeeApiUrl);
          const ieeePapers = (ieeeResponse.data.articles || []).map((article) => ({
            title: article.title || "Untitled",
            authors: article.authors?.map((a) => a.full_name).join(", ") || "Unknown Authors",
            year: article.publication_year || "Unknown Year",
            abstract: article.abstract || "Abstract not available",
            doi: article.doi || "No DOI available",
            url: constructIeeeUrl(article.doi),
            is_simulated: false,
          }));
  
          papers = [...papers, ...ieeePapers].slice(0, 3); // Combine and limit to 3 papers
          console.log("IEEE API papers fetched:", ieeePapers);
        } catch (apiErr) {
          console.error("IEEE API error:", apiErr.message);
          // Fallback with simulated IEEE-like papers if API fails
          const fallbackPapers = generateFallbackPapers(selectedMol.newSmiles);
          papers = [...papers, ...fallbackPapers].slice(0, 3);
          toast.warning("IEEE API failed; using simulated fallback papers.");
        }
      } else {
        console.warn("No valid API key provided. Using Gemini-generated and fallback papers.");
        const fallbackPapers = generateFallbackPapers(selectedMol.newSmiles);
        papers = [...papers, ...fallbackPapers].slice(0, 3);
      }
  
      // Ensure all URLs are workable
      papers = papers.map((paper) => {
        let validUrl = paper.url;
        if (!validUrl || !validUrl.startsWith("https://ieeexplore.ieee.org/document/")) {
          const doi = paper.doi || generateUniqueDoi(paper.year || 2025);
          validUrl = constructIeeeUrl(doi);
        }
        return {
          ...paper,
          doi: paper.doi || validUrl.split("/document/")[1],
          url: validUrl,
          is_valid: true, // Assume constructed URLs are structurally valid
        };
      });
  
      setResearchPapers(papers);
      await savePapers(selectedMol, papers);
      toast.success("Related research papers fetched and saved successfully!");
    } catch (err) {
      console.error("Error fetching research:", err);
      setError(err.response?.data?.message || "Failed to fetch research papers");
      toast.error("Failed to fetch research papers");
    } finally {
      setLoading(false);
    }
  };
  
  // Helper functions
  const constructIeeeUrl = (doi) => {
    if (!doi) return `https://ieeexplore.ieee.org/document/${generateUniqueDoi()}`;
    const doiSuffix = doi.includes("10.") ? doi.split("/").pop() : doi.replace(/[^0-9]/g, "");
    return `https://ieeexplore.ieee.org/document/${doiSuffix}`;
  };
  
  const generateUniqueDoi = (year = 2025, index = 0) => {
    const randomId = Math.floor(Math.random() * 1000000) + index;
    return `10.1109/TBME.${year}.${randomId}`;
  };
  
  const generateFallbackPapers = (smiles) => {
    // Simulate realistic IEEE papers based on SMILES features
    const inferredFeatures = smiles.includes("c") ? "aromatic rings" : "aliphatic chains";
    return [
      {
        title: `Computational Analysis of ${inferredFeatures} in Drug Design`,
        authors: "R. Johnson, S. Lee, T. Brown",
        year: "2022",
        abstract: `This paper explores the role of ${inferredFeatures} derived from SMILES structures like ${smiles.substring(0, 10)}... in predicting drug efficacy, focusing on cheminformatics techniques.`,
        doi: "10.1109/TBME.2022.987654",
        url: "https://ieeexplore.ieee.org/document/987654",
        is_simulated: true,
      },
      {
        title: `Synthesis and Properties of Novel ${inferredFeatures} Compounds`,
        authors: "M. Davis, P. Kim, Q. Zhang",
        year: "2020",
        abstract: `Investigates synthesis pathways for compounds with ${inferredFeatures}, leveraging SMILES-based modeling for therapeutic applications.`,
        doi: "10.1109/TBME.2020.876543",
        url: "https://ieeexplore.ieee.org/document/876543",
        is_simulated: true,
      },
    ];
  };
  
  const handleGeneratePaperClick = async () => {
    if (!selectedTitle || !selectedSmiles) {
      toast.error("Please select both a title and SMILES string");
      return;
    }

    const paperExists = await checkIfGeneratedPaperExists(selectedTitle, selectedSmiles);
    if (paperExists) {
      toast("A generated research paper for this molecule already exists. Redirecting to Saved Generated Research Papers.", {
        type: "info",
      });
      setActiveTab("savedGenerated");
      return;
    }

    setLoading(true);
    setError(null);
    setGeneratedPaper(null);

    const selectedMol = molecules.find(
      (mol) => mol.newmoleculetitle === selectedTitle && mol.newSmiles === selectedSmiles
    );
    if (!selectedMol) {
      setError("Selected molecule not found");
      setLoading(false);
      return;
    }

    const geminiPrompt = `
    You are an expert in cheminformatics, molecular modeling, and academic writing. The molecule provided is newly generated, so exact research papers about it do not exist. However, you will infer potential applications, synthesis, and significance based on its structure, properties, and related literature. Use the following molecule details to generate a research paper in IEEE format:

    - Molecule Title: "${selectedMol.newmoleculetitle}"
    - SMILES String: "${selectedMol.newSmiles}"
    - IUPAC Name: "${selectedMol.newIupacName}"
    - Conversion Details (Synthesis Pathway): "${selectedMol.conversionDetails}"
    - Potential Diseases (Therapeutic Targets): "${selectedMol.potentialDiseases}"
    - Additional Information (Pharmacokinetics, Toxicity, or Mechanism of Action): "${selectedMol.information}"

    ### Instructions:
    1. **Generate a Research Paper in IEEE Format**:
       - **Title**: Refine the molecule title to reflect its scientific importance.
       - **Authors**: Use the placeholder "Author Name" (to be replaced by the user's name later).
       - **Abstract**: Write a 150-200 word abstract summarizing the molecule’s chemical structure, synthesis, pharmacological relevance, and connection to existing research.
       - **Keywords**: Provide 4-5 relevant keywords such as "cheminformatics," "drug discovery," "bioactivity prediction," "synthetic pathway."
       - **Introduction**: Explain the motivation behind this molecule’s study, its significance in drug discovery, and how it differs from known compounds (1-2 paragraphs).
       - **Methodology**: Describe the synthesis process (use the conversion details) and structural analysis techniques (e.g., NMR, Mass Spectrometry, Molecular Docking).
       - **Results and Discussion**: Discuss the molecule's properties (SMILES, IUPAC name), potential biological activity (based on potential diseases), and compare with existing drugs (using structure-activity relationships).
       - **Conclusion**: Summarize key findings, highlight the molecule’s novelty, and suggest future research directions.
       - **References**: Generate 3-5 **authentic references** based on related molecules or therapeutic areas from sources like PubMed, Google Scholar, or ACS Publications. Format references in IEEE style (e.g., [1] Author(s), "Title," Journal, vol., no., pp., year.).

    2. **Return the Paper in a Structured JSON Format**:
       Provide the following in a structured JSON format **and nothing else**:
       {
         "title": "Refined Paper Title",
         "authors": "Author Name",
         "abstract": "Abstract text summarizing the molecule’s significance.",
         "keywords": ["keyword1", "keyword2", "keyword3", "keyword4"],
         "introduction": "Introduction text providing scientific context.",
         "methodology": "Synthesis pathway and structural analysis methods.",
         "resultsAndDiscussion": "Molecular properties, bioactivity, and comparison with existing drugs.",
         "conclusion": "Summary of key insights and future directions.",
         "references": [
           "[1] Author(s), 'Title,' Journal, vol., no., pp., year.",
           "[2] Author(s), 'Title,' Journal, vol., no., pp., year."
         ]
       }
  `;

    try {
      const geminiResponse = await axiosInstance.post("/protein/proxy/gemini", {
        prompt: geminiPrompt,
      });
      const geminiContent = geminiResponse.data.content;
      const jsonMatch = geminiContent.match(/{[\s\S]*}/);
      if (!jsonMatch) {
        throw new Error("No valid JSON found in Gemini response");
      }

      const parsedPaper = JSON.parse(jsonMatch[0]);
      parsedPaper.authors = user.name || "Author Name";
      setGeneratedPaper(parsedPaper);
      await saveGeneratedPaper(selectedMol, parsedPaper);
      toast.success("Research paper generated and saved successfully!");
    } catch (err) {
      console.error("Error generating research paper:", err);
      setError(err.response?.data?.message || "Failed to generate research paper");
      toast.error("Failed to generate research paper");
    } finally {
      setLoading(false);
    }
  };

  const exportToPDF = (paper) => {
    if (!paper) {
      toast.error("No paper to download");
      return;
    }

    const pdf = new jsPDF("p", "mm", "a4");
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 15;
    let yPosition = margin;

    // Helper function to check if we need a new page
    const checkPageBreak = (requiredHeight) => {
      if (yPosition + requiredHeight > pageHeight - margin) {
        pdf.addPage();
        yPosition = margin;
      }
    };

    // Helper function to add text with pagination
    const addTextWithPagination = (text, fontSize, x, y, maxWidth) => {
      pdf.setFontSize(fontSize);
      pdf.setFont("times", "normal");
      const lines = pdf.splitTextToSize(text, maxWidth);
      const lineHeight = fontSize * 0.4;

      lines.forEach((line) => {
        checkPageBreak(lineHeight);
        pdf.text(line, x, yPosition);
        yPosition += lineHeight;
      });

      return yPosition;
    };

    // Add Title (Centered, with wrapping for long titles)
    pdf.setFontSize(16);
    pdf.setFont("times", "bold");
    const titleLines = pdf.splitTextToSize(paper.title, pageWidth - 2 * margin);
    const titleLineHeight = 16 * 0.4; // Line height for title
    titleLines.forEach((line) => {
      checkPageBreak(titleLineHeight);
      const lineWidth = pdf.getTextWidth(line);
      pdf.text(line, (pageWidth - lineWidth) / 2, yPosition); // Center each line
      yPosition += titleLineHeight;
    });
    yPosition += 10; // Extra spacing after title

    // Add Authors (Centered)
    pdf.setFontSize(12);
    pdf.setFont("times", "normal");
    const authorsWidth = pdf.getTextWidth(paper.authors);
    checkPageBreak(8);
    pdf.text(paper.authors, (pageWidth - authorsWidth) / 2, yPosition);
    yPosition += 15;

    // Add Abstract
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("Abstract", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.abstract, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 10;

    // Add Keywords
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("Keywords", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.keywords.join(", "), 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    // Add Introduction
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("I. Introduction", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.introduction, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    // Add Methodology
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("II. Methodology", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.methodology, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    // Add Results and Discussion
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("III. Results and Discussion", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.resultsAndDiscussion, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    // Add Conclusion
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("IV. Conclusion", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.conclusion, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    // Add References
    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("References", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.setFont("times", "normal");
    paper.references.forEach((ref) => {
      checkPageBreak(6);
      const lines = pdf.splitTextToSize(ref, pageWidth - 2 * margin);
      lines.forEach((line) => {
        checkPageBreak(5);
        pdf.text(line, margin, yPosition);
        yPosition += 5;
      });
      yPosition += 2;
    });

    // Save the PDF
    pdf.save(`${paper.title.replace(/\s+/g, "_")}.pdf`);
    toast.success("PDF downloaded successfully!");
  };

  const handleTabChange = (tab) => {
    setActiveTab(tab);
    setResearchPapers([]);
    setResearchSummary("");
    setGeneratedPaper(null);
    setError(null);
  };

  if (checkingAuth) {
    return <div className="text-center py-10 text-gray-600">Verifying authentication...</div>;
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-100">
        <div className="text-center p-6 bg-white rounded-lg shadow-lg">
          <p className="text-gray-600 mb-4">Please log in to access AI Research Generator</p>
          <button
            className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
            onClick={() => (window.location.href = "/login")}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-100 py-12 px-4 sm:px-6 lg:px-8">
      <div className="max-w-7xl mx-auto">
        <h1 className="text-4xl font-bold text-blue-700 mb-10 text-center">
          AI Research Generator
          <p className="text-xs  p-1 text-blue-700 font-semibold">(Powered by Gemini)</p>

        </h1>

        <div className="flex justify-center mb-8 space-x-4">
          {["related", "saved", "generate", "savedGenerated"].map((tab) => (
            <button
              key={tab}
              className={`px-6 py-2 rounded-lg font-medium transition-all duration-300 ${
                activeTab === tab
                  ? "bg-blue-600 text-white shadow-md"
                  : "bg-gray-200 text-gray-700 hover:bg-gray-300"
              }`}
              onClick={() => handleTabChange(tab)}
            >
              {tab === "related" && "Related Research Papers"}
              {tab === "saved" && "Saved Research Papers"}
              {tab === "generate" && "Generate Research Paper"}
              {tab === "savedGenerated" && "Saved Generated Papers"}
            </button>
          ))}
        </div>

        <div className="bg-white p-8 rounded-xl shadow-lg border border-gray-200">
          {activeTab === "related" && (
            <>
              <h2 className="text-2xl font-semibold text-blue-700 mb-6">
                Related Research Papers
              </h2>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    Select Molecule Title
                  </label>
                  <select
                    value={selectedTitle}
                    onChange={(e) => setSelectedTitle(e.target.value)}
                    className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No titles available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newmoleculetitle))].map((title) => (
                        <option key={title} value={title}>
                          {title}
                        </option>
                      ))
                    )}
                  </select>
                </div>
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    Select SMILES String
                  </label>
                  <select
                    value={selectedSmiles}
                    onChange={(e) => setSelectedSmiles(e.target.value)}
                    className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newSmiles))].map((smiles) => (
                        <option key={smiles} value={smiles}>
                          {smiles}
                        </option>
                      ))
                    )}
                  </select>
                </div>
              </div>

              <button
                onClick={handleResearchClick}
                disabled={loading || !selectedTitle || !selectedSmiles}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300"
              >
                {loading ? "Fetching Research..." : "Fetch Related Research"}
              </button>

              {/* {error && (
                <div className="mt-6 bg-red-50 border border-red-300 text-red-700 px-4 py-3 rounded-lg flex justify-between items-center">
                  <p>{error}</p>
                  <button
                    className="text-red-700 underline hover:text-red-900"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </div>
              )} */}

              {(researchSummary || researchPapers.length > 0) && (
                <div className="mt-8 bg-blue-50 p-6 rounded-xl border border-blue-200">
                  <h3 className="text-lg font-semibold text-blue-700 mb-4">
                    Related Research Information
                  </h3>

                  {researchSummary && (
                    <div className="mb-8">
                      <h4 className="text-md font-semibold text-gray-800 mb-2">Research Context</h4>
                      <p className="text-gray-700 text-sm leading-relaxed">{researchSummary}</p>
                    </div>
                  )}

                  {researchPapers.length > 0 && (
                    <div className="mb-8">
                      <h4 className="text-md font-semibold text-gray-800 mb-4">
                        Newly Fetched Research Papers
                      </h4>
                      <div className="space-y-8">
                        {researchPapers.map((paper, index) => (
                          <div key={index} className="border-l-4 border-blue-500 pl-4">
                            <h5 className="text-lg font-bold text-gray-900 mb-2 uppercase">
                              {paper.title}
                            </h5>
                            <p className="text-gray-700 text-sm mb-1">
                              <span className="font-bold">Authors:</span> {paper.authors}
                            </p>
                            <p className="text-gray-700 text-sm mb-1">
                              <span className="font-bold">Published:</span> {paper.year}
                            </p>
                            <p className="text-gray-700 text-sm mb-1">
                              <span className="font-bold">Abstract:</span> {paper.abstract}
                            </p>
                            {paper.doi && paper.url !== "No URL available" && (
                              <p className="text-blue-600 text-sm">
                                <span className="font-bold">DOI:</span>{" "}
                                <a
                                  href={paper.url}
                                  target="_blank"
                                  rel="noopener noreferrer"
                                  className="underline hover:text-blue-500"
                                >
                                  {paper.doi}
                                </a>
                              </p>
                            )}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </>
          )}

          {activeTab === "saved" && (
            <>
              <h2 className="text-2xl font-semibold text-blue-700 mb-6">
                Saved Research Papers
              </h2>
              {savedPapers.length > 0 ? (
                <div className="space-y-10">
                  {savedPapers.map((entry, index) => {
                    if (!entry.molecule || !entry.molecule.title || !entry.molecule.smiles) {
                      return null;
                    }
                    return (
                      <div key={index} className="bg-blue-50 p-6 rounded-xl border border-blue-200">
                        <h3 className="text-lg font-semibold text-gray-800 mb-4">
                          Molecule: {entry.molecule.title} (SMILES: {entry.molecule.smiles})
                        </h3>
                        <div className="space-y-8">
                          {entry.papers.map((paper, paperIndex) => (
                            <div key={paperIndex} className="border-l-4 border-blue-500 pl-4">
                              <h5 className="text-lg font-bold text-gray-900 mb-2 uppercase">
                                {paper.title}
                              </h5>
                              <p className="text-gray-700 text-sm mb-1">
                                <span className="font-bold">Authors:</span> {paper.authors}
                              </p>
                              <p className="text-gray-700 text-sm mb-1">
                                <span className="font-bold">Published:</span> {paper.year}
                              </p>
                              <p className="text-gray-700 text-sm mb-1">
                                <span className="font-bold">Abstract:</span> {paper.abstract}
                              </p>
                              {paper.doi && paper.url !== "No URL available" && (
                                <p className="text-blue-600 text-sm">
                                  <span className="font-bold">DOI:</span>{" "}
                                  <a
                                    href={paper.url}
                                    target="_blank"
                                    rel="noopener noreferrer"
                                    className="underline hover:text-blue-500"
                                  >
                                    {paper.doi}
                                  </a>
                                </p>
                              )}
                            </div>
                          ))}
                        </div>
                      </div>
                    );
                  })}
                </div>
              ) : (
                <p className="text-gray-600 text-center">No saved research papers found.</p>
              )}
            </>
          )}

          {activeTab === "generate" && (
            <>
              <h2 className="text-2xl font-semibold text-blue-700 mb-6">
                Generate Research Paper
              </h2>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    Select Molecule Title
                  </label>
                  <select
                    value={selectedTitle}
                    onChange={(e) => setSelectedTitle(e.target.value)}
                    className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No titles available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newmoleculetitle))].map((title) => (
                        <option key={title} value={title}>
                          {title}
                        </option>
                      ))
                    )}
                  </select>
                </div>
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    Select SMILES String
                  </label>
                  <select
                    value={selectedSmiles}
                    onChange={(e) => setSelectedSmiles(e.target.value)}
                    className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
                    disabled={loading || molecules.length === 0}
                  >
                    {molecules.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      [...new Set(molecules.map((m) => m.newSmiles))].map((smiles) => (
                        <option key={smiles} value={smiles}>
                          {smiles}
                        </option>
                      ))
                    )}
                  </select>
                </div>
              </div>

              <button
                onClick={handleGeneratePaperClick}
                disabled={loading || !selectedTitle || !selectedSmiles}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300"
              >
                {loading ? "Generating Paper..." : "Generate Research Paper"}
              </button>

              {error && (
                <div className="mt-6 bg-red-50 border border-red-300 text-red-700 px-4 py-3 rounded-lg flex justify-between items-center">
                  <p>{error}</p>
                  <button
                    className="text-red-700 underline hover:text-red-900"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </div>
              )}

              {generatedPaper && (
                <div className="mt-8 p-6 rounded-xl border border-blue-200 bg-blue-50">
                  <h3 className="text-lg font-semibold text-blue-700 mb-4">
                    Generated Research Paper (IEEE Format)
                  </h3>

                  <div className="flex justify-end mb-4">
                    <button
                      onClick={() => exportToPDF(generatedPaper)}
                      className="px-4 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-all duration-300"
                    >
                      Download as PDF
                    </button>
                  </div>

                  <div className="space-y-6">
                    <h4 className="text-2xl font-bold text-center text-gray-900">
                      {generatedPaper.title}
                    </h4>
                    <p className="text-center text-gray-700">{generatedPaper.authors}</p>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">Abstract</h5>
                      <p className="text-sm leading-relaxed text-gray-700">{generatedPaper.abstract}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">Keywords</h5>
                      <p className="text-sm text-gray-700">{generatedPaper.keywords.join(", ")}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">I. Introduction</h5>
                      <p className="text-sm leading-relaxed text-gray-700">{generatedPaper.introduction}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">II. Methodology</h5>
                      <p className="text-sm leading-relaxed text-gray-700">{generatedPaper.methodology}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">III. Results and Discussion</h5>
                      <p className="text-sm leading-relaxed text-gray-700">{generatedPaper.resultsAndDiscussion}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">IV. Conclusion</h5>
                      <p className="text-sm leading-relaxed text-gray-700">{generatedPaper.conclusion}</p>
                    </div>

                    <div>
                      <h5 className="text-lg font-semibold mb-2 text-gray-800">References</h5>
                      <ul className="list-none text-sm space-y-2 text-gray-700">
                        {generatedPaper.references.map((ref, index) => (
                          <li key={index}>{ref}</li>
                        ))}
                      </ul>
                    </div>
                  </div>
                </div>
              )}
            </>
          )}

          {activeTab === "savedGenerated" && (
            <>
              <h2 className="text-2xl font-semibold text-blue-700 mb-6">
                Saved Generated Research Papers
              </h2>
              {savedGeneratedPapers.length > 0 ? (
                <div className="space-y-10">
                  {savedGeneratedPapers.map((entry, index) => {
                    if (!entry.molecule || !entry.molecule.title || !entry.molecule.smiles || !entry.paper) {
                      return null;
                    }
                    return (
                      <div key={index} className="p-6 rounded-xl border border-blue-200 bg-blue-50">
                        <h3 className="text-lg font-semibold text-gray-800 mb-4">
                          Molecule: {entry.molecule.title} (SMILES: {entry.molecule.smiles})
                        </h3>
                        <div className="flex justify-end mb-4">
                          <button
                            onClick={() => exportToPDF(entry.paper)}
                            className="px-4 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-all duration-300"
                          >
                            Download as PDF
                          </button>
                        </div>
                        <div className="space-y-6">
                          <h4 className="text-2xl font-bold text-center text-gray-900">
                            {entry.paper.title}
                          </h4>
                          <p className="text-center text-gray-700">{entry.paper.authors}</p>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">Abstract</h5>
                            <p className="text-sm leading-relaxed text-gray-700">{entry.paper.abstract}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">Keywords</h5>
                            <p className="text-sm text-gray-700">{entry.paper.keywords.join(", ")}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">I. Introduction</h5>
                            <p className="text-sm leading-relaxed text-gray-700">{entry.paper.introduction}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">II. Methodology</h5>
                            <p className="text-sm leading-relaxed text-gray-700">{entry.paper.methodology}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">III. Results and Discussion</h5>
                            <p className="text-sm leading-relaxed text-gray-700">{entry.paper.resultsAndDiscussion}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">IV. Conclusion</h5>
                            <p className="text-sm leading-relaxed text-gray-700">{entry.paper.conclusion}</p>
                          </div>

                          <div>
                            <h5 className="text-lg font-semibold mb-2 text-gray-800">References</h5>
                            <ul className="list-none text-sm space-y-2 text-gray-700">
                              {entry.paper.references.map((ref, refIndex) => (
                                <li key={refIndex}>{ref}</li>
                              ))}
                            </ul>
                          </div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              ) : (
                <p className="text-gray-600 text-center">No saved generated research papers found.</p>
              )}
            </>
          )}
        </div>
      </div>
    </div>
  );
}

export default Airesearchgenerator;