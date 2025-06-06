import { useState, useEffect, useRef } from "react";
import { Bot, Microscope, Atom, Brain, X, Mic, MicOff, Volume2, VolumeX, Dna, Pill } from "lucide-react";
import { useAuthStore } from "../../Store/auth.store.js";
import { useNavigate, useLocation } from "react-router-dom";

export default function Chatbot() {
  const [isListening, setIsListening] = useState(false);
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [transcript, setTranscript] = useState("");
  const [response, setResponse] = useState("");
  const [isSpeaking, setIsSpeaking] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [isRecognizing, setIsRecognizing] = useState(false);
  const [conversationHistory, setConversationHistory] = useState([]);
  const [showTooltip, setShowTooltip] = useState("");
  const [mode, setMode] = useState(null); // null, "doubt", or "beginner"
  const [subMode, setSubMode] = useState(null); // null, "waiting_for_question", "resolving_doubt", "tour_complete"
  const [targetRoute, setTargetRoute] = useState(null); // Track intended route
  const recognitionRef = useRef(null);
  const isCleaningUp = useRef(false);
  const conversationRef = useRef(null);
  const hasPromptedQuestion = useRef(false);
  const recognitionState = useRef("idle"); // idle, starting, running, stopping
  const lastStartAttempt = useRef(0); // For debouncing
  const lastProcessedRoute = useRef(null); // Prevent duplicate useEffect triggers
  const intendedMode = useRef(null); // Track intended mode
  const modeLock = useRef(false); // Prevent mode resets
  const { user } = useAuthStore();
  const navigate = useNavigate();
  const location = useLocation();

  // Normalize path for route comparison
  const normalizePath = (path) => {
    console.log(`Normalizing path: ${path}`);
    if (!path || typeof path !== "string") {
      console.warn(`Invalid path: ${path}, returning default '/dashboard'`);
      return "/dashboard";
    }
    return path.replace(/\/+$/, "").toLowerCase();
  };

  // Log state changes for debugging
  useEffect(() => {
    console.log(`State update: isListening=${isListening}, isRecognizing=${isRecognizing}, isSpeaking=${isSpeaking}, mode=${mode}, subMode=${subMode}, path=${location.pathname}, recognitionState=${recognitionState.current}, targetRoute=${targetRoute}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
  }, [isListening, isRecognizing, isSpeaking, mode, subMode, location.pathname, targetRoute]);

  // Forcefully set subMode when hasPromptedQuestion is true
  useEffect(() => {
    if (hasPromptedQuestion.current && subMode !== "waiting_for_question" && mode === "beginner" && !modeLock.current) {
      console.warn(`SubMode mismatch: subMode=${subMode}, forcing to waiting_for_question`);
      setSubMode("waiting_for_question");
    }
  }, [subMode, mode]);

  // Check for BrowserRouter
  useEffect(() => {
    if (!navigate || typeof navigate !== "function") {
      console.error("Navigate function is undefined. Ensure Chatbot is wrapped in BrowserRouter.");
    }
    if (!location || !location.pathname) {
      console.error("Location or pathname is undefined. Check React Router setup.");
    }
  }, [navigate, location]);

  // Initialize Web Speech API
  const createRecognition = () => {
    const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition;
    if (!SpeechRecognition) {
      console.error("SpeechRecognition API not supported.");
      return null;
    }
    const recognition = new SpeechRecognition();
    recognition.lang = "en-US";
    recognition.interimResults = false;
    recognition.maxAlternatives = 1;
    recognition.continuous = true;
    return recognition;
  };

  // Configure recognition handlers
  const setupRecognition = (recognition) => {
    recognition.onresult = async (event) => {
      const speechResult = event.results[event.results.length - 1][0].transcript.trim().toLowerCase();
      console.log(`Captured transcript: ${speechResult}, mode=${mode}, subMode=${subMode}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
      setTranscript(speechResult);
      setConversationHistory((prev) => [...prev, { type: "user", text: speechResult }]);

      if (speechResult.includes("stop") || speechResult.includes("thank you")) {
        console.log("Stop command detected, stopping conversation");
        await stopConversation();
        return;
      }

      // Restore mode synchronously
      let currentMode = mode;
      if (!currentMode && intendedMode.current) {
        console.warn(`Mode is null, restoring from intendedMode: ${intendedMode.current}`);
        currentMode = intendedMode.current;
        setMode(currentMode);
        setSubMode(currentMode === "beginner" ? "waiting_for_question" : null);
        await new Promise((resolve) => setTimeout(resolve, 0));
      }

      if (!currentMode || currentMode !== intendedMode.current || modeLock.current) {
        console.warn(`Invalid state: mode=${currentMode}, subMode=${subMode}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}, prompting restart`);
        modeLock.current = true;
        await speakResponse("Please start a new conversation in beginner or question mode.");
        modeLock.current = false;
        setIsListening(true);
        await startRecognition();
        return;
      }

      modeLock.current = true;
      try {
        if (currentMode === "beginner") {
          console.log(`Processing in beginner mode with subMode=${subMode}`);
          await handleBeginnerMode(speechResult, currentMode);
        } else if (currentMode === "doubt") {
          console.log("Processing in doubt mode");
          await fetchResponse(speechResult);
        } else {
          console.warn(`Unexpected mode: ${currentMode}, resetting`);
          setMode(null);
          setSubMode(null);
          intendedMode.current = null;
          await speakResponse("Please start a new conversation in beginner or question mode.");
        }
      } catch (error) {
        console.error(`Error in onresult: ${error}`);
        await speakResponse("Sorry, an error occurred. Please try again.");
      } finally {
        modeLock.current = false;
      }
    };

    recognition.onend = async () => {
      console.log(`Recognition ended, isListening: ${isListening}, isCleaningUp: ${isCleaningUp.current}, isSpeaking: ${isSpeaking}, isRecognizing: ${isRecognizing}, recognitionState: ${recognitionState.current}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";
      if (isListening && !isCleaningUp.current && !isSpeaking && mode && mode === intendedMode.current) {
        await startRecognition();
      }
    };

    recognition.onerror = async (event) => {
      console.error(`Speech recognition error: ${event.error}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";
      if (event.error === "no-speech" && isListening && !isCleaningUp.current && mode) {
        console.log("No speech detected, retrying recognition");
        const message = "I didn't hear anything.";
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }]);
        await speakResponse(message);
        if (isListening && !isCleaningUp.current && !isSpeaking && mode) {
          await startRecognition();
        }
      } else if (isListening && !isCleaningUp.current && !isSpeaking && mode) {
        console.log("Retrying recognition after error");
        await startRecognition();
      }
    };
  };

  // Start speech recognition with robust state management
  const startRecognition = async () => {
    const now = Date.now();
    if (now - lastStartAttempt.current < 1000) {
      console.log("Debouncing startRecognition attempt");
      return;
    }
    lastStartAttempt.current = now;

    console.log(`Attempting to start recognition, isListening: ${isListening}, isCleaningUp: ${isCleaningUp.current}, isSpeaking: ${isSpeaking}, isRecognizing: ${isRecognizing}, recognitionState=${recognitionState.current}, recognitionRef.current=${!!recognitionRef.current}`);
    
    // Ensure recognition is stopped before starting
    if (recognitionRef.current && (isRecognizing || recognitionState.current !== "idle")) {
      try {
        recognitionState.current = "stopping";
        recognitionRef.current.stop();
        await new Promise((resolve) => setTimeout(resolve, 200));
        recognitionState.current = "idle";
        setIsRecognizing(false);
      } catch (error) {
        console.error(`Error stopping recognition before start: ${error}`);
      }
    }

    if (!recognitionRef.current || isCleaningUp.current || isSpeaking) {
      console.log("Cannot start recognition: not initialized, cleaning up, or speaking");
      recognitionState.current = "idle";
      setIsRecognizing(false);
      return;
    }

    try {
      recognitionState.current = "starting";
      recognitionRef.current.start();
      recognitionState.current = "running";
      setIsRecognizing(true);
      setIsListening(true);
      console.log("Speech recognition started successfully");
    } catch (error) {
      console.error(`Error starting recognition: ${error}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";
      if (!isCleaningUp.current && !isSpeaking && mode) {
        setTimeout(() => startRecognition(), 1000);
      }
    }
  };

  // Start conversation
  const startConversation = async (selectedMode) => {
    console.log(`Starting conversation with mode: ${selectedMode}`);
    isCleaningUp.current = false;
    setIsPanelOpen(true);
    setIsListening(false);
    setIsRecognizing(false);
    setTranscript("");
    setResponse("");
    setConversationHistory([]);
    setMode(selectedMode);
    setSubMode(null);
    setTargetRoute(null);
    intendedMode.current = selectedMode;
    modeLock.current = false;
    hasPromptedQuestion.current = false;
    recognitionState.current = "idle";
    lastProcessedRoute.current = null;

    if (recognitionRef.current) {
      try {
        recognitionState.current = "stopping";
        recognitionRef.current.stop();
        await new Promise((resolve) => setTimeout(resolve, 200));
        recognitionState.current = "idle";
      } catch (error) {
        console.error(`Error stopping existing recognition: ${error}`);
      }
      recognitionRef.current = null;
    }

    recognitionRef.current = createRecognition();
    if (!recognitionRef.current) {
      const errorMessage = "Speech recognition is not supported in this browser.";
      setResponse(errorMessage);
      setConversationHistory((prev) => [{ type: "jarvis", text: errorMessage }]);
      await speakResponse(errorMessage);
      setIsPanelOpen(false);
      return;
    }
    setupRecognition(recognitionRef.current);

    let welcomeMessage = "";
    if (selectedMode === "doubt") {
      welcomeMessage = "I am Jarvis, here to assist with your drug discovery questions. Ask away!";
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }]);
      await speakResponse(welcomeMessage);
      setIsListening(true);
      await startRecognition();
    } else if (selectedMode === "beginner") {
      welcomeMessage =
        "Welcome to your drug discovery journey! I'm Jarvis, and I'll guide you through the process. " +
        "Let's start with the Protein Structure Mutation tool.";
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }]);
      await speakResponse(welcomeMessage);
      const target = "/dashboard/protein-structure-mutation";
      setTargetRoute(target);
      console.log(`Navigating to ${target}`);
      navigate(target, { replace: true, state: { fromJarvis: true } });
    }
  };

  // Stop conversation
  const stopConversation = async () => {
    console.log("Stopping conversation");
    isCleaningUp.current = true;
    setIsListening(false);
    setIsRecognizing(false);
    setMode(null);
    setSubMode(null);
    setTargetRoute(null);
    intendedMode.current = null;
    modeLock.current = false;
    hasPromptedQuestion.current = false;
    recognitionState.current = "idle";
    lastProcessedRoute.current = null;

    const goodbyeMessage = "Goodbye. Thank you for using Jarvis.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: goodbyeMessage }]);
    await speakResponse(goodbyeMessage);

    if (recognitionRef.current) {
      try {
        recognitionState.current = "stopping";
        recognitionRef.current.stop();
        await new Promise((resolve) => setTimeout(resolve, 200));
        recognitionState.current = "idle";
      } catch (error) {
        console.error(`Error stopping recognition: ${error}`);
      }
      recognitionRef.current = null;
    }

    setTimeout(() => {
      setIsPanelOpen(false);
      isCleaningUp.current = false;
    }, 200);
  };

  // Toggle speech
  const toggleSpeech = async () => {
    if (isSpeaking) {
      setIsSpeaking(false);
      setIsListening(true);
      await startRecognition();
    } else if (response) {
      await speakResponse(response);
    }
  };

  // Fetch response from Gemini API
  const fetchResponse = async (query) => {
    console.log(`Fetching response for query: ${query}`);
    setIsLoading(true);

    const problemStatement = `The process of drug discovery is time-consuming, expensive, and often inefficient, with a high rate of failure in clinical trials. Traditional methods rely heavily on trial and error, requiring years of research and significant financial investment. Additionally, the complexity of biological systems and the vast chemical space make it challenging to identify promising drug candidates efficiently. Generative AI, with its ability to analyze large datasets, predict molecular interactions, and generate novel compounds, has the potential to revolutionize this process. However, there is a lack of accessible, user-friendly tools that leverage generative AI to assist researchers in accelerating drug discovery while reducing costs and improving success rates.`;

    const prompt = `You are Jarvis, an AI assistant specialized in drug discovery. You are trained on the following problem statement: "${problemStatement}". Answer the user's query related to drug discovery, from basic molecule or drug names to complex information about the process, generative AI applications, or challenges. Provide a concise, accurate response suitable for researchers. Query: ${query}`;

    try {
      const response = await fetch(
        "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=AIzaSyDyujm50dHMYdn1Vv50dDDqcAhgUqCOuUGU",
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            contents: [{ parts: [{ text: prompt }] }],
          }),
        }
      );

      if (!response.ok) {
        throw new Error(`HTTP error! Status: ${response.status}`);
      }

      const data = await response.json();
      const text = data.candidates[0].content.parts[0].text;
      console.log(`API response: ${text}`);
      setResponse(text);
      setConversationHistory((prev) => [...prev, { type: "jarvis", text }]);
      await speakResponse(text);

      const jarvisEntry = {
        title: `Query: ${query.slice(0, 50)}`,
        questions: query,
        answers: text,
        userId: user?._id || "66f172e9c7e124f5c4b3c2d1",
        createdAt: new Date().toISOString(),
      };
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]");
      storedEntries.push(jarvisEntry);
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries));
      console.log("Stored entry in localStorage:", jarvisEntry);
    } catch (error) {
      console.error(`Error fetching Gemini API: ${error}`);
      const errorMessage = "Sorry, I encountered an error. Please try again.";
      setResponse(errorMessage);
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }]);
      await speakResponse(errorMessage);

      const jarvisEntry = {
        title: `Query error: ${query.slice(0, 50)}`,
        questions: query,
        answers: `Error: ${error.message}`,
        userId: user?._id || "66f172e9c7e124f5c4b3c2d1",
        createdAt: new Date().toISOString(),
      };
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]");
      storedEntries.push(jarvisEntry);
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries));
      console.log("Stored error entry in localStorage:", jarvisEntry);
    } finally {
      setIsLoading(false);
      if (mode === "beginner" && subMode === "resolving_doubt") {
        console.log("Setting subMode to waiting_for_question after fetch");
        setSubMode("waiting_for_question");
        hasPromptedQuestion.current = true;
        const questionPrompt = "Do you have any more questions about this tool? Say 'yes' or 'no'.";
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: questionPrompt }]);
        await speakResponse(questionPrompt);
      }
    }
  };

  // Speak response using Web Speech API
  const speakResponse = async (text) => {
    return new Promise((resolve) => {
      console.log(`Speaking response: ${text}, mode=${mode}, subMode=${subMode}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
      if (recognitionRef.current && isListening && !isCleaningUp.current) {
        try {
          recognitionState.current = "stopping";
          recognitionRef.current.stop();
          console.log("Recognition paused for speech");
          setTimeout(() => {
            recognitionState.current = "idle";
            setIsRecognizing(false);
          }, 200);
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`);
          recognitionState.current = "idle";
          setIsRecognizing(false);
        }
      }

      modeLock.current = true;
      const utterance = new SpeechSynthesisUtterance(text);
      utterance.lang = "en-US";
      utterance.onstart = () => {
        setIsSpeaking(true);
        console.log(`Speech started, isListening=${isListening}, mode=${mode}, subMode=${subMode}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
      };
      utterance.onend = async () => {
        setIsSpeaking(false);
        console.log(`Speech ended, mode=${mode}, subMode=${subMode}, hasPrompted=${hasPromptedQuestion.current}, isListening=${isListening}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
        if (text.includes("Goodbye.")) {
          setIsListening(false);
          setIsRecognizing(false);
          recognitionState.current = "idle";
          setIsPanelOpen(false);
          isCleaningUp.current = true;
        } else if (mode && mode === intendedMode.current) {
          setIsListening(true);
          if (!recognitionRef.current && !isCleaningUp.current) {
            recognitionRef.current = createRecognition();
            if (recognitionRef.current) {
              setupRecognition(recognitionRef.current);
            }
          }
          await startRecognition();
        }
        modeLock.current = false;
        resolve();
      };
      utterance.onerror = (error) => {
        console.error(`Speech synthesis error: ${error}`);
        setIsSpeaking(false);
        setIsListening(mode && mode === intendedMode.current ? true : false);
        setIsRecognizing(false);
        recognitionState.current = "idle";
        if (mode && mode === intendedMode.current) {
          startRecognition();
        }
        modeLock.current = false;
        resolve();
      };
      window.speechSynthesis.speak(utterance);
    });
  };

  // Handle beginner mode logic
const handleBeginnerMode = async (speechResult, currentMode) => {
  console.log(`Beginner mode handler - Transcript: "${speechResult}" | Mode: ${currentMode} | SubMode: ${subMode} | Path: ${location.pathname} | IntendedMode: ${intendedMode.current} | ModeLock: ${modeLock.current}`);

  // Validate mode
  if (currentMode !== "beginner" || currentMode !== intendedMode.current) {
    console.warn(`Invalid mode in handleBeginnerMode: ${currentMode}, resetting to beginner`);
    setMode("beginner");
    intendedMode.current = "beginner";
    setSubMode("waiting_for_question");
    hasPromptedQuestion.current = true;
    modeLock.current = false;
    await speakResponse("Do you have any questions about this tool? Say 'yes' or 'no'.");
    await startRecognition();
    return;
  }

  const routes = [
    "/dashboard/protein-structure-mutation",
    "/dashboard/cost-estimation",
    "/dashboard/ai-naming",
  ];
  const normalizedCurrentPath = normalizePath(location.pathname);
  const currentRouteIndex = routes.findIndex((route) => normalizePath(route) === normalizedCurrentPath);
  console.log(`Current path: ${normalizedCurrentPath}, Current route index: ${currentRouteIndex}`);

  const yesVariants = ["yes", "yeah", "yep", "sure", "okay"];
  const noVariants = ["no", "nope", "nah", "not really"];

  try {
    // Stop recognition before processing
    if (recognitionRef.current && isListening) {
      recognitionState.current = "stopping";
      recognitionRef.current.stop();
      await new Promise((resolve) => setTimeout(resolve, 200));
      recognitionState.current = "idle";
      setIsRecognizing(false);
      setIsListening(false);
      console.log("Recognition stopped for processing");
    }

    if (subMode === "waiting_for_question" || !subMode) {
      if (yesVariants.some((v) => speechResult.toLowerCase().includes(v))) {
        console.log("Detected 'yes', switching to resolving_doubt");
        setSubMode("resolving_doubt");
        modeLock.current = false;
        await speakResponse("Please ask your question about the current tool.");
        await startRecognition();
      } else if (noVariants.some((v) => speechResult.toLowerCase().includes(v))) {
        console.log("Detected 'no', proceeding to next route");
        let nextRoute;
        if (currentRouteIndex === -1) {
          console.warn(`Invalid current route: ${location.pathname}, navigating to cost-estimation`);
          nextRoute = "/dashboard/cost-estimation";
        } else if (currentRouteIndex === 1) {
          console.warn(`Current route is cost-estimation, navigating to ai-naming`);
          nextRoute = "/dashboard/ai-naming";
        } else {
          const nextRouteIndex = currentRouteIndex + 1;
          console.log(`Next route index: ${nextRouteIndex}, Routes length: ${routes.length}`);
          if (nextRouteIndex < routes.length) {
            nextRoute = routes[nextRouteIndex];
          } else {
            console.log("Tour complete, prompting user choice");
            const completeMsg =
              "You've completed the drug discovery process tour! Would you like to start over or switch to asking questions? Say 'start over' or 'ask questions'.";
            setConversationHistory((prev) => [...prev, { type: "jarvis", text: completeMsg }]);
            setSubMode("tour_complete");
            hasPromptedQuestion.current = false;
            modeLock.current = false;
            await speakResponse(completeMsg);
            await startRecognition();
            return;
          }
        }

        // Perform navigation
        console.log(`Navigating to: ${nextRoute}`);
        setConversationHistory((prev) => [
          ...prev,
          { type: "jarvis", text: `Navigating to ${getRouteName(nextRoute)}` },
        ]);
        setTargetRoute(nextRoute);
        navigate(nextRoute, { replace: true, state: { fromJarvis: true } });
        modeLock.current = false;

        // Defer state reset until navigation is confirmed
        setSubMode("waiting_for_question");
        hasPromptedQuestion.current = true;

        await startRecognition();
      } else {
        console.log("Unrecognized response, prompting for 'yes' or 'no'");
        setSubMode("waiting_for_question");
        hasPromptedQuestion.current = true;
        modeLock.current = false;
        await speakResponse("Do you have any questions about this tool? Say 'yes' or 'no'.");
        await startRecognition();
      }
    } else if (subMode === "resolving_doubt") {
      console.log("Processing doubt query in beginner mode");
      await fetchResponse(speechResult);
      await startRecognition();
    } else if (subMode === "tour_complete") {
      if (speechResult.toLowerCase().includes("start over")) {
        console.log("Restarting tour, navigating to first route");
        setConversationHistory((prev) => [
          ...prev,
          { type: "jarvis", text: `Navigating to ${getRouteName(routes[0])}` },
        ]);
        setTargetRoute(routes[0]);
        navigate(routes[0], { replace: true, state: { fromJarvis: true } });
        setSubMode("waiting_for_question");
        hasPromptedQuestion.current = true;
        modeLock.current = false;
        await startRecognition();
      } else if (speechResult.toLowerCase().includes("ask questions")) {
        console.log("Switching to doubt mode");
        setMode("doubt");
        intendedMode.current = "doubt";
        setSubMode(null);
        hasPromptedQuestion.current = false;
        modeLock.current = false;
        await speakResponse("Switching to question mode. Ask your drug discovery question.");
        await startRecognition();
      } else {
        console.log("Unrecognized response, prompting for 'start over' or 'ask questions'");
        modeLock.current = false;
        await speakResponse("Please say 'start over' to begin again or 'ask questions' to switch modes.");
        await startRecognition();
      }
    } else {
      console.warn(`Unexpected subMode: ${subMode}, resetting to waiting_for_question`);
      setSubMode("waiting_for_question");
      hasPromptedQuestion.current = true;
      modeLock.current = false;
      await speakResponse("Do you have any questions about this tool? Say 'yes' or 'no'.");
      await startRecognition();
    }
  } catch (error) {
    console.error(`Beginner mode error: ${error}`);
    setSubMode("waiting_for_question");
    hasPromptedQuestion.current = true;
    modeLock.current = false;
    await speakResponse("Sorry, I encountered an error. Please say 'yes' or 'no' to continue.");
    await startRecognition();
  }
};

  // Speak form fields and tab info based on route
  const speakFieldsForRoute = async (route) => {
    console.log(`Speaking fields for route: ${route}, mode=${mode}, subMode=${subMode}, intendedMode=${intendedMode.current}, modeLock=${modeLock.current}`);
    if (!route || typeof route !== "string") {
      console.error(`Invalid route: ${route}, defaulting to /dashboard/protein-structure-mutation`);
      route = "/dashboard/protein-structure-mutation";
      setTargetRoute(route);
      navigate(route, { replace: true, state: { fromJarvis: true } });
    }
    const normalizedRoute = normalizePath(route);
    let message = "";
    switch (normalizedRoute) {
      case normalizePath("/dashboard/protein-structure-mutation"):
        message =
          "You are in the Protein Structure Mutation tool. This tool generates new drug molecules by combining two existing ones. " +
          "Please fill in: New Molecule Title, a name for your molecule; " +
          "Search First Molecule, enter the SMILES string for the first molecule; " +
          "Search Second Molecule, enter the SMILES string for the second molecule. " +
          "Click 'Generate Molecule' to create a new compound.";
        break;
      case normalizePath("/dashboard/cost-estimation"):
        message =
          "You are in the Cost Estimation tool. This tool estimates synthesis and production costs. " +
          "Select a SMILES string from your generated molecules in the dropdown, then click 'Estimate Cost'.";
        break;
      case normalizePath("/dashboard/ai-naming"):
        message =
          "You are in the AI Naming Suggestion tool. This tool generates drug names. " +
          "Select a Molecule Title and a SMILES string from the dropdowns, then click 'Start Prediction'.";
        break;
      default:
        message = "It seems you're on an unrecognized page. Let's return to the drug discovery tour.";
        setConversationHistory((prev) => [
          ...prev,
          { type: "jarvis", text: message },
        ]);
        setTargetRoute("/dashboard/protein-structure-mutation");
        navigate("/dashboard/protein-structure-mutation", { replace: true, state: { fromJarvis: true } });
        await speakResponse(message);
        return;
    }
    setConversationHistory((prev) => [
      ...prev,
      { type: "jarvis", text: `Now at ${getRouteName(route)}. ${message}` },
    ]);
    await speakResponse(`Now at ${getRouteName(route)}. ${message}`);

    setSubMode("waiting_for_question");
    hasPromptedQuestion.current = true;
    const questionPrompt = "Do you have any questions about this tool? Say 'yes' or 'no'.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: questionPrompt }]);
    await speakResponse(questionPrompt);
  };

  // Get route name for user-friendly display
  const getRouteName = (route) => {
    if (!route) {
      console.warn(`Invalid route for getRouteName: ${route}, returning 'Unknown'`);
      return "Unknown";
    }
    switch (normalizePath(route)) {
      case normalizePath("/dashboard/protein-structure-mutation"):
        return "Protein Structure Mutation";
      case normalizePath("/dashboard/cost-estimation"):
        return "Cost Estimation";
      case normalizePath("/dashboard/ai-naming"):
        return "AI Naming Suggestion";
      default:
        return "Unknown";
    }
  };

  // Handle route changes and validate navigation
  useEffect(() => {
    console.log(`useEffect triggered, mode: ${mode}, subMode: ${subMode}, path: ${location.pathname}, fromJarvis: ${!!location.state?.fromJarvis}, targetRoute: ${targetRoute}, intendedMode: ${intendedMode.current}, modeLock: ${modeLock.current}`);
    if (!location || !location.pathname) {
      console.error("Location or pathname is undefined, redirecting to default route");
      const defaultRoute = "/dashboard/protein-structure-mutation";
      setTargetRoute(defaultRoute);
      navigate(defaultRoute, { replace: true, state: { fromJarvis: true } });
      return;
    }
    if (!mode && location.state?.fromJarvis && !modeLock.current) {
      console.warn(`Mode is null with fromJarvis, resetting to beginner, intendedMode: ${intendedMode.current}`);
      setMode("beginner");
      intendedMode.current = "beginner";
      setSubMode(null);
      hasPromptedQuestion.current = false;
    }
    if ((mode === "beginner" || mode === null) && location.state?.fromJarvis && location.pathname && !modeLock.current) {
      const normalizedPath = normalizePath(location.pathname);
      const validRoutes = [
        "/dashboard/protein-structure-mutation",
        "/dashboard/cost-estimation",
        "/dashboard/ai-naming",
      ].map(normalizePath);
      console.log(`Route changed to: ${normalizedPath}, validRoutes: ${validRoutes}, lastProcessedRoute: ${lastProcessedRoute.current}, targetRoute: ${targetRoute}`);
      if (validRoutes.includes(normalizedPath) && normalizedPath !== normalizePath(lastProcessedRoute.current)) {
        if (targetRoute && normalizePath(targetRoute) !== normalizedPath) {
          console.warn(`Route mismatch: current=${normalizedPath}, target=${targetRoute}, retrying navigation`);
          navigate(targetRoute, { replace: true, state: { fromJarvis: true } });
          return;
        }
        console.log(`Triggering speakFieldsForRoute for: ${location.pathname}`);
        lastProcessedRoute.current = location.pathname;
        speakFieldsForRoute(location.pathname);
      } else if (!validRoutes.includes(normalizedPath)) {
        console.warn(`Invalid route for beginner mode: ${location.pathname}, redirecting to first route`);
        const firstRoute = "/dashboard/protein-structure-mutation";
        setTargetRoute(firstRoute);
        navigate(firstRoute, { replace: true, state: { fromJarvis: true } });
        setConversationHistory((prev) => [
          ...prev,
          { type: "jarvis", text: `Navigating to ${getRouteName(firstRoute)}` },
        ]);
        setMode("beginner");
        intendedMode.current = "beginner";
        setSubMode(null);
        hasPromptedQuestion.current = false;
        lastProcessedRoute.current = firstRoute;
        speakFieldsForRoute(firstRoute);
      }
    }
  }, [location, mode, subMode, navigate, targetRoute]);

  // Monitor targetRoute for navigation confirmation
  useEffect(() => {
    if (targetRoute && normalizePath(location.pathname) !== normalizePath(targetRoute)) {
      console.warn(`Target route mismatch: current=${location.pathname}, target=${targetRoute}, retrying navigation`);
      navigate(targetRoute, { replace: true, state: { fromJarvis: true } });
    }
  }, [location.pathname, targetRoute, navigate]);

  // Scroll to bottom of conversation
  useEffect(() => {
    if (conversationRef.current) {
      conversationRef.current.scrollTop = conversationRef.current.scrollHeight;
    }
  }, [conversationHistory]);

  // Clean up speech synthesis and recognition on unmount
  useEffect(() => {
    return () => {
      console.log("Cleaning up speech synthesis and recognition");
      isCleaningUp.current = true;
      window.speechSynthesis.cancel();
      if (recognitionRef.current) {
        try {
          recognitionState.current = "stopping";
          recognitionRef.current.stop();
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`);
        }
        recognitionRef.current = null;
      }
      setIsRecognizing(false);
      recognitionState.current = "idle";
      modeLock.current = false;
    };
  }, []);

  return (
    <div className="font-sans">
      {/* Trigger Button */}
      <div className="relative">
        <button
          onClick={() => setIsPanelOpen(true)}
          disabled={isPanelOpen}
          onMouseEnter={() => setShowTooltip("trigger")}
          onMouseLeave={() => setShowTooltip("")}
          className={`fixed bottom-6 right-6 h-16 w-16 rounded-full shadow-lg transition-all transform hover:scale-105 flex items-center justify-center ${
            isPanelOpen
              ? "bg-gradient-to-r from-red-500 to-red-600 hover:from-red-600 hover:to-red-700"
              : "bg-gradient-to-br from-blue-500 to-blue-600 hover:from-blue-600 to-blue-700"
          }`}
        >
          {isPanelOpen ? <MicOff className="h-6 w-6 text-white" /> : <Mic className="h-6 w-6 text-white" />}
        </button>
        {showTooltip === "trigger" && (
          <div className="absolute bottom-full right-0 mb-2 px-2 py-1 bg-gray-800 text-white text-xs rounded shadow-lg whitespace-nowrap">
            <div>{isPanelOpen ? "Close Jarvis" : "Ask Jarvis"}</div>
          </div>
        )}
      </div>

      {/* Pointing Message */}
      {!isPanelOpen && (
        <div className="fixed bottom-24 right-7">
          <div
            className="relative bg-white rounded-lg shadow-lg p-3 mb-2"
            style={{
              animation: "pulse-scale 2s infinite ease-in-out",
              transformOrigin: "center",
            }}
          >
            <style>
              {`
                @keyframes pulse-scale {
                  0%, 100% { transform: scale(1); }
                  50% { transform: scale(1.05); }
                }
              `}
            </style>
            <div className="absolute right-4 bottom-2 transform translate-y-full w-3 h-3 bg-white rotate-45"></div>
            <div className="flex items-center space-x-2">
              <div className="bg-blue-500/20 p-1.5 rounded-full">
                <Mic className="h-4 w-4 text-blue-600" />
              </div>
              <div>
                <div className="text-sm font-medium text-gray-800">Drug Discovery Voice Assistant</div>
                <div className="text-xs text-gray-500">Click to speak with Jarvis</div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Right Side Panel */}
      {isPanelOpen && (
        <div className="fixed top-0 right-0 h-full w-full md:w-96 bg-white shadow-2xl z-50 flex flex-col border-l border-gray-200 animate-in slide-in-from-right duration-300">
          <div className="flex items-center justify-between p-4 border-b border-gray-200">
            <div className="flex items-center space-x-2">
              <Dna className="h-6 w-6 text-blue-600" />
              <h2 className="text-xl font-semibold text-gray-800">Jarvis</h2>
            </div>
            <button
              onClick={stopConversation}
              className="rounded-full h-8 w-8 bg-gray-100 hover:bg-gray-200 text-gray-600 flex items-center justify-center"
            >
              <X className="h-4 w-4" />
              <span className="sr-only">Close</span>
            </button>
          </div>

          {!mode ? (
            <div className="flex-1 p-6 flex flex-col items-center justify-center space-y-4">
              <div className="text-center">
                <h3 className="text-lg font-semibold text-gray-800 mb-2">Welcome to Jarvis</h3>
                <p className="text-sm text-gray-600">
                  Choose an option to start your drug discovery journey.
                </p>
              </div>
              <button
                onClick={() => startConversation("beginner")}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors text-sm font-medium"
              >
                Are you a beginner?
              </button>
              <button
                onClick={() => startConversation("doubt")}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors text-sm font-medium"
              >
                Ask a question
              </button>
            </div>
          ) : (
            <div className="flex flex-col h-full">
              {/* Jarvis Visualization */}
              <div className="p-6 bg-gradient-to-br from-blue-50 to-cyan-50 border-b border-gray-200">
                <div className="flex flex-col items-center text-center space-y-2">
                  <div className="relative">
                    <div
                      className="absolute inset-0 rounded-full bg-blue-500/10 animate-pulse"
                      style={{ animationDuration: "3s" }}
                    ></div>
                    <div className="h-24 w-24 rounded-full bg-gradient-to-br from-blue-400 to-blue-500 flex items-center justify-center">
                      <Dna className="h-12 w-12 text-white" />
                    </div>
                  </div>
                  <h1 className="text-2xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-cyan-600">
                    JARVIS
                  </h1>
                  <p className="text-blue-700 text-sm">Drug Discovery Assistant</p>
                  <div className="flex items-center mt-2 space-x-1">
                    <div
                      className={`h-2 w-2 rounded-full ${isListening ? "bg-green-500 animate-pulse" : "bg-red-500"}`}
                    ></div>
                    <span className="text-xs text-gray-600">{isListening ? "Listening" : "Idle"}</span>
                  </div>
                </div>
                <div className="w-64 mt-4">
                  <img
                    src="/jarvis.jpg"
                    alt="Molecular AI Visualization"
                    className="rounded-xl shadow-md max-w-full h-auto object-cover border border-gray-200"
                    onError={(e) => {
                      e.target.onerror = null;
                      e.target.src = "/jarvis.jpg";
                    }}
                  />
                </div>
                <div className="w-full space-y-4 mt-4">
                  <div className="grid grid-cols-4 gap-2">
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Atom className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">Molecules</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Brain className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">AI</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Microscope className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">Research</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Pill className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">Drugs</span>
                    </div>
                  </div>
                  <div className="text-xs text-gray-500 text-center">Powered by Gemini</div>
                </div>
              </div>

              {/* Conversation History */}
              <div ref={conversationRef} className="flex-1 overflow-y-auto p-6 space-y-4 custom-scrollbar">
                {conversationHistory.map((message, index) => (
                  <div key={index} className={`flex ${message.type === "user" ? "justify-end" : "justify-start"}`}>
                    <div
                      className={`max-w-[80%] rounded-2xl px-4 py-2 ${
                        message.type === "user"
                          ? "bg-blue-500 text-white rounded-tr-none"
                          : "bg-white text-gray-800 rounded-tl-none border border-gray-200 shadow-sm"
                      }`}
                    >
                      {message.type === "jarvis" && (
                        <div className="flex items-center space-x-2 mb-1 pb-1 border-b border-gray-200">
                          <Bot className="h-3 w-3 text-blue-600" />
                          <span className="text-xs font-medium text-blue-600">JARVIS</span>
                        </div>
                      )}
                      <p className="text-sm whitespace-pre-wrap">{message.text}</p>
                    </div>
                  </div>
                ))}
                {isLoading && (
                  <div className="flex justify-start">
                    <div className="max-w-[80%] rounded-2xl px-4 py-2 bg-white border border-gray-200 shadow-sm">
                      <div className="flex items-center space-x-2 mb-1 pb-1 border-b border-gray-200">
                        <Bot className="h-3 w-3 text-blue-600" />
                        <span className="text-xs font-medium text-blue-600">JARVIS</span>
                      </div>
                      <div className="flex space-x-1 items-center">
                        <div
                          className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
                          style={{ animationDelay: "0ms" }}
                        ></div>
                        <div
                          className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
                          style={{ animationDelay: "150ms" }}
                        ></div>
                        <div
                          className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
                          style={{ animationDelay: "300ms" }}
                        ></div>
                      </div>
                    </div>
                  </div>
                )}
              </div>

              {/* Input Area */}
              <div className="p-6 border-t border-gray-200">
                <div className="flex items-center justify-between p-4 rounded-xl bg-white border border-gray-200 shadow-sm">
                  <div className="flex items-center space-x-2 text-sm text-gray-600">
                    {isListening ? (
                      <>
                        <div className="relative">
                          <div className="absolute inset-0 bg-green-500/20 rounded-full animate-ping"></div>
                          <Mic className="h-5 w-5 text-green-600" />
                        </div>
                        <span>Listening... Say "stop" or "thank you" to end</span>
                      </>
                    ) : (
                      <>
                        <MicOff className="h-5 w-5 text-red-500" />
                        <span>Voice recognition paused</span>
                      </>
                    )}
                  </div>
                  <div className="flex space-x-2">
                    <button
                      onClick={toggleSpeech}
                      onMouseEnter={() => setShowTooltip("speech")}
                      onMouseLeave={() => setShowTooltip("")}
                      className="h-8 w-8 rounded-full bg-gray-100 hover:bg-gray-200 flex items-center justify-center text-gray-600"
                    >
                      {isSpeaking ? <VolumeX className="h-4 w-4" /> : <Volume2 className="h-4 w-4" />}
                    </button>
                    <button
                      onClick={isListening ? stopConversation : () => startConversation(mode || "doubt")}
                      className={`px-4 py-2 rounded-lg text-sm font-medium text-white ${
                        isListening ? "bg-red-500 hover:bg-red-600" : "bg-blue-600 hover:bg-blue-700"
                      } transition-colors`}
                    >
                      {isListening ? "Stop" : "Start"}
                    </button>
                  </div>
                </div>
                <div className="mt-2 text-xs text-center text-gray-500">
                  Jarvis can assist with molecular structures, drug discovery processes, and AI applications in
                  pharmaceutical research
                </div>
              </div>
            </div>
          )}
          <style jsx>{`
            .custom-scrollbar::-webkit-scrollbar {
              width: 6px;
            }
            .custom-scrollbar::-webkit-scrollbar-track {
              background: #f1f1f1;
              border-radius: 10px;
            }
            .custom-scrollbar::-webkit-scrollbar-thumb {
              background: #d1d5db;
              border-radius: 10px;
            }
            .custom-scrollbar::-webkit-scrollbar-thumb:hover {
              background: #9ca3af;
            }
          `}</style>
        </div>
      )}
    </div>
  );
}