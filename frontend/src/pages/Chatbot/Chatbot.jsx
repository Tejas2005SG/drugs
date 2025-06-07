import { useState, useEffect, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import {
  Bot,
  Microscope,
  Atom,
  Brain,
  X,
  Mic,
  MicOff,
  Dna,
  Pill,
  ChevronRight,
  ChevronLeft,
  CheckCircle,
  HelpCircle,
  ArrowLeft,
} from "lucide-react";
import { useAuthStore } from "../../Store/auth.store.js";

export default function Chatbot() {
  const [isPanelOpen, setIsPanelOpen] = useState(false);
  const [response, setResponse] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [isSpeaking, setIsSpeaking] = useState(false);
  const [isListening, setIsListening] = useState(false);
  const [isRecognizing, setIsRecognizing] = useState(false);
  const [transcript, setTranscript] = useState("");
  const [conversationHistory, setConversationHistory] = useState([]);
  const [tooltipVisible, setTooltipVisible] = useState("");
  const [mode, setMode] = useState(null); // null, "doubt", or "beginner"
  const [subMode, setSubMode] = useState(null); // null, "waiting_for_selection", "waiting_for_tool_selection", "waiting_for_step_selection", "explaining_tab", "waiting_for_question", "resolving_doubt", "tour_complete", "process_ended"
  const [targetRoute, setTargetRoute] = useState(null);
  const [currentToolIndex, setCurrentToolIndex] = useState(0);
  const [selectedOption, setSelectedOption] = useState(null); // "dashboard" or "drugDiscovery"
  const [selectedTool, setSelectedTool] = useState(null); // Track manually selected tool or step
  const [completedSteps, setCompletedSteps] = useState(
    JSON.parse(localStorage.getItem("drugDiscoveryProgress") || "[]")
  ); // Persist progress
  const recognitionRef = useRef(null);
  const isCleaningUp = useRef(false);
  const recognitionState = useRef("idle");
  const conversationRef = useRef(null);
  const lastProcessedRoute = useRef(null);
  const recognitionLock = useRef(false);
  const { user } = useAuthStore();
  const navigate = useNavigate();
  const location = useLocation();

  // Define dashboard routes and drug discovery steps
  const dashboardRoutes = [
    { path: "/dashboard/protein-structure", name: "Protein Structure", icon: "Atom" },
    { path: "/dashboard/protein-structure-mutation", name: "Protein Structure Mutation", icon: "Dna" },
    { path: "/dashboard/cost-estimation", name: "Cost Estimation", icon: "Pill" },
    { path: "/dashboard/ai-research-paper-generator", name: "AI Research Paper Generator", icon: "Brain" },
    { path: "/dashboard/ai-driven-target-prediction", name: "AI-Driven Target Prediction", icon: "Microscope" },
    { path: "/dashboard/ai-naming", name: "AI Naming Suggestion", icon: "Bot" },
  ];

  const drugDiscoverySteps = [
    {
      path: "/dashboard/protein-structure-mutation",
      name: "Step 1: Protein Mutation",
      description: "Combine molecules to create new compounds",
      icon: "Dna",
    },
    {
      path: "/dashboard/cost-estimation",
      name: "Step 2: Cost Estimation",
      description: "Calculate synthesis and production costs",
      icon: "Pill",
    },
    {
      path: "/dashboard/ai-naming",
      name: "Step 3: AI Naming",
      description: "Generate potential drug names",
      icon: "Bot",
    },
    {
      path: "/dashboard/ai-research-paper-generator",
      name: "Step 4: Research Paper",
      description: "Create research paper drafts",
      icon: "Brain",
    },
    {
      path: "/dashboard/summary",
      name: "Step 5: Summary",
      description: "Review your drug discovery project",
      icon: "Microscope",
    },
  ];

  // Suggested questions per tool/step
  const suggestedQuestions = {
    "/dashboard/protein-structure": ["How do I visualize a protein in 3D?", "What file formats are supported?"],
    "/dashboard/protein-structure-mutation": ["How does mutation affect drug design?", "Can I specify mutation sites?"],
    "/dashboard/cost-estimation": ["What factors are included in cost estimation?", "How accurate is the cost prediction?"],
    "/dashboard/ai-research-paper-generator": ["Can I customize the paper format?", "What data is required for the draft?"],
    "/dashboard/ai-driven-target-prediction": ["How does AI predict drug targets?", "What is the accuracy of predictions?"],
    "/dashboard/ai-naming": ["How are drug names generated?", "Can I filter names by criteria?"],
    "/dashboard/summary": ["What is included in the project summary?", "Can I export the summary?"],
  };

  // Map icon names to components
  const iconMap = {
    Atom: <Atom className="h-4 w-4" />,
    Dna: <Dna className="h-4 w-4" />,
    Pill: <Pill className="h-4 w-4" />,
    Brain: <Brain className="h-4 w-4" />,
    Microscope: <Microscope className="h-4 w-4" />,
    Bot: <Bot className="h-4 w-4" />,
  };

  // Normalize path for route comparison
  const normalizePath = (path) => {
    if (!path || typeof path !== "string") {
      console.warn(`Invalid path: ${path}, returning default '/dashboard'`);
      return "/dashboard";
    }
    return path.replace(/\/+$/, "").toLowerCase();
  };

  // Initialize Web Speech API for recognition (doubt mode only)
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

  // Configure recognition handlers (doubt mode only)
  const setupRecognition = (recognition) => {
    recognition.onstart = () => {
      console.log("Recognition started");
      setIsRecognizing(true);
      recognitionState.current = "active";
    };

    recognition.onresult = async (event) => {
      const speechResult = event.results[event.results.length - 1][0].transcript.trim();
      console.log(`Captured transcript: "${speechResult}"`);
      setTranscript(speechResult);
      setConversationHistory((prev) => [...prev, { type: "user", text: speechResult }]);

      if (speechResult.toLowerCase().includes("stop") || speechResult.toLowerCase().includes("thank you")) {
        console.log("Stop command detected, stopping conversation");
        await stopConversation();
        return;
      }

      await fetchResponse(speechResult);
    };

    recognition.onend = async () => {
      console.log(`Recognition ended - isListening: ${isListening}, isSpeaking: ${isSpeaking}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";
      if (!isCleaningUp.current && !isSpeaking && isListening && mode === "doubt") {
        console.log("Restarting recognition");
        await startRecognition();
      }
    };

    recognition.onerror = async (event) => {
      console.error(`Speech recognition error: ${event.error}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";
      if (event.error === "no-speech" && !isCleaningUp.current && isListening && mode === "doubt") {
        const message = "I didn't hear anything. Please try again.";
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }]);
        await speakResponse(message);
        await startRecognition();
      } else {
        recognitionRef.current = null;
        if (!isCleaningUp.current && isListening && mode === "doubt") {
          await startRecognition();
        }
      }
    };
  };

  // Start speech recognition (doubt mode only)
  const startRecognition = async (retryCount = 0) => {
    if (recognitionLock.current || isCleaningUp.current || isSpeaking) {
      console.log("Recognition blocked");
      if (isSpeaking && retryCount < 3) {
        await new Promise((resolve) => setTimeout(resolve, 500));
        await startRecognition(retryCount + 1);
      }
      return;
    }

    recognitionLock.current = true;

    try {
      if (recognitionRef.current && (isRecognizing || recognitionState.current === "active")) {
        recognitionRef.current.stop();
        await new Promise((resolve) => setTimeout(resolve, 200));
      }

      if (!recognitionRef.current) {
        recognitionRef.current = createRecognition();
        if (!recognitionRef.current) {
          recognitionLock.current = false;
          const errorMessage = "Speech recognition is not supported on this device.";
          setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }]);
          await speakResponse(errorMessage);
          setIsListening(false);
          return;
        }
        setupRecognition(recognitionRef.current);
      }

      recognitionRef.current.start();
      setIsListening(true);
      setIsRecognizing(true);
    } catch (error) {
      console.error(`Error starting recognition: ${error}`);
      setIsRecognizing(false);
      recognitionState.current = "idle";

      if (retryCount < 2) {
        recognitionRef.current = null;
        await new Promise((resolve) => setTimeout(resolve, 200));
        await startRecognition(retryCount + 1);
      } else {
        const message = "Sorry, I couldn't start speech recognition. Please try again.";
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }]);
        await speakResponse(message);
        setIsListening(false);
      }
    } finally {
      recognitionLock.current = false;
    }
  };

  // Speak response using Web Speech API
  const speakResponse = async (text) => {
    return new Promise((resolve) => {
      if (!text || typeof text !== "string") {
        console.error("Invalid text for speech synthesis:", text);
        setIsSpeaking(false);
        resolve();
        return;
      }

      if (recognitionRef.current && isRecognizing) {
        try {
          recognitionRef.current.stop();
          setIsRecognizing(false);
          recognitionState.current = "idle";
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`);
        }
      }

      window.speechSynthesis.cancel();

      const utterance = new SpeechSynthesisUtterance(text);
      utterance.lang = "en-US";
      utterance.rate = 0.9;
      utterance.onstart = () => {
        setIsSpeaking(true);
      };
      utterance.onend = async () => {
        setIsSpeaking(false);
        if (mode === "doubt" && isListening && !isCleaningUp.current) {
          await startRecognition();
        }
        resolve();
      };
      utterance.onerror = (event) => {
        console.error(`Speech synthesis error: ${event.error}`);
        setIsSpeaking(false);
        resolve();
      };
      window.speechSynthesis.speak(utterance);
    });
  };

  // Start conversation
  const startConversation = async (selectedMode) => {
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
    setCurrentToolIndex(0);
    setSelectedTool(null);
    setSelectedOption(null);
    setCompletedSteps(JSON.parse(localStorage.getItem("drugDiscoveryProgress") || "[]"));

    if (recognitionRef.current) {
      try {
        recognitionRef.current.stop();
      } catch (error) {
        console.error(`Error stopping recognition: ${error}`);
      }
      recognitionRef.current = null;
    }

    let welcomeMessage = "";
    if (selectedMode === "doubt") {
      welcomeMessage = "I am Jarvis, here to assist with your medical and pharmaceutical questions. Please speak your question.";
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }]);
      await speakResponse(welcomeMessage);
      setIsListening(true);
      await startRecognition();
    } else if (selectedMode === "beginner") {
      welcomeMessage = "Welcome to your drug discovery journey! I'm Jarvis, and I'll guide you through the dashboard tools or the drug discovery process. Please select an option.";
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }]);
      await speakResponse(welcomeMessage);
      setSubMode("waiting_for_selection");
    }
  };

  // Stop conversation
  const stopConversation = async () => {
    isCleaningUp.current = true;
    setIsListening(false);
    setIsRecognizing(false);
    setMode(null);
    setSubMode(null);
    setTargetRoute(null);
    setCurrentToolIndex(0);
    setSelectedTool(null);
    setSelectedOption(null);
    setCompletedSteps([]);

    const goodbyeMessage = "Goodbye. Thank you for using Jarvis.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: goodbyeMessage }]);
    await speakResponse(goodbyeMessage);

    if (recognitionRef.current) {
      try {
        recognitionRef.current.stop();
      } catch (error) {
        console.error(`Error stopping recognition: ${error}`);
      }
      recognitionRef.current = null;
    }

    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]));
    setTimeout(() => {
      setIsPanelOpen(false);
      isCleaningUp.current = false;
    }, 200);
  };

  // Reset to initial selection screen
  const resetToInitialSelection = async () => {
    const prompt = "Returning to the initial selection screen.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
    setMode(null);
    setSubMode(null);
    setTargetRoute(null);
    setCurrentToolIndex(0);
    setSelectedTool(null);
    setSelectedOption(null);
    setCompletedSteps([]);
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]));
  };

  // Fetch response from Gemini API
  const fetchResponse = async (query) => {
    setIsLoading(true);
    const contextInfo = selectedTool ? `Current tool context: ${selectedTool.name}. ` : "";
    const problemStatement = `You are Jarvis, an AI assistant specialized in medical and pharmaceutical domains. ${contextInfo}Answer the user's query comprehensively and accurately. Keep your response concise but informative. Query: ${query}`;

    try {
      const response = await fetch(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${import.meta.env.VITE_API_KEY}`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            contents: [{ parts: [{ text: problemStatement }] }],
          }),
        }
      );

      if (!response.ok) {
        throw new Error(`HTTP error! Status: ${response.status}`);
      }

      const data = await response.json();
      const text = data.candidates[0].content.parts[0].text;
      setResponse(text);
      setConversationHistory((prev) => [...prev, { type: "jarvis", text }]);
      await speakResponse(text);

      // Store the Q&A in local storage
      const jarvisEntry = {
        title: `Query: ${query.slice(0, 50)}`,
        questions: query,
        answers: text,
        toolContext: selectedTool?.name || "General",
        userId: user?._id || "66f172e9c7e124f5c4b3c2d1",
        createdAt: new Date().toISOString(),
      };
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]");
      storedEntries.push(jarvisEntry);
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries));

      if (mode === "beginner" && (subMode === "resolving_doubt" || subMode === "waiting_for_question")) {
        setSubMode("waiting_for_question");
        const followUpPrompt = "Would you like to ask another question or select a different tool/step?";
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: followUpPrompt }]);
        await speakResponse(followUpPrompt);
      }
    } catch (error) {
      console.error(`Error fetching Gemini API: ${error}`);
      const errorMessage = "Sorry, I encountered an error while processing your question. Would you like to try again or select another tool/step?";
      setResponse(errorMessage);
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }]);
      await speakResponse(errorMessage);

      if (mode === "beginner") {
        setSubMode("waiting_for_question");
      }
    } finally {
      setIsLoading(false);
    }
  };

  // Handle dashboard tools selection
  const handleDashboardSelection = async () => {
    setSelectedOption("dashboard");
    setSubMode("waiting_for_tool_selection");
    const prompt = "You chose to explore dashboard tools. Please select a tool to learn about.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
  };

  // Handle drug discovery selection
  const handleDrugDiscoverySelection = async () => {
    setSelectedOption("drugDiscovery");
    setSubMode("waiting_for_step_selection");
    const prompt = "You chose the drug discovery process. Please select the current step to begin.";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
  };

  // Handle manual tool selection with delay
  const handleToolSelection = async (tool) => {
    setSelectedTool({ ...tool, icon: tool.icon });
    const prompt = `You selected ${tool.name}. Now navigating to that page for more information.`;
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
    await new Promise((resolve) => setTimeout(resolve, 1000)); // 1-second delay
    navigate(tool.path, { replace: true, state: { fromJarvis: true } });
    setTargetRoute(tool.path);
    setSubMode("explaining_tab");
    setCurrentToolIndex(dashboardRoutes.findIndex((t) => t.path === tool.path));
    lastProcessedRoute.current = tool.path;
    await speakFieldsForRoute(tool.path);
  };

  // Handle manual step selection with delay
  const handleStepSelection = async (step) => {
    setSelectedTool({ ...step, icon: step.icon });
    const prompt = `You selected ${step.name}. Now navigating to that page for more information.`;
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
    await new Promise((resolve) => setTimeout(resolve, 1000)); // 1-second delay
    navigate(step.path, { replace: true, state: { fromJarvis: true } });
    setTargetRoute(step.path);
    setSubMode("explaining_tab");
    setCurrentToolIndex(drugDiscoverySteps.findIndex((s) => s.path === step.path));
    setCompletedSteps((prev) => {
      const newSteps = [...prev, step.path];
      localStorage.setItem("drugDiscoveryProgress", JSON.stringify(newSteps));
      return newSteps;
    });
    await speakFieldsForRoute(step.path);
  };

  // Proceed to next step in the process
  const proceedToNextStep = async () => {
    setCompletedSteps((prev) => {
      const newSteps = [...prev, drugDiscoverySteps[currentToolIndex].path];
      localStorage.setItem("drugDiscoveryProgress", JSON.stringify(newSteps));
      return newSteps;
    });
    const completionPrompt = `Step ${currentToolIndex + 1} completed!`;
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: completionPrompt }]);
    await speakResponse(completionPrompt);

    if (currentToolIndex < drugDiscoverySteps.length - 1) {
      const nextStep = drugDiscoverySteps[currentToolIndex + 1];
      const prompt = `Moving to ${nextStep.name}.`;
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
      await speakResponse(prompt);
      await new Promise((resolve) => setTimeout(resolve, 1000)); // 1-second delay
      navigate(nextStep.path, { replace: true, state: { fromJarvis: true } });
      setTargetRoute(nextStep.path);
      setSubMode("explaining_tab");
      setCurrentToolIndex(currentToolIndex + 1);
      lastProcessedRoute.current = nextStep.path;
      await speakFieldsForRoute(nextStep.path);
    } else {
      const prompt = "You've reached the end of the drug discovery process. Click 'End Process' to complete or explore other options.";
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
      await speakResponse(prompt);
      setSubMode("tour_complete");
    }
  };

  // End drug discovery process
  const endProcess = async () => {
    const congratsPrompt = "Congratulations! You've successfully completed the drug discovery process! Would you like to restart the process or explore dashboard tools?";
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: congratsPrompt }]);
    await speakResponse(congratsPrompt);
    setSubMode("process_ended");
    setCurrentToolIndex(0);
    setCompletedSteps([]);
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]));
  };

  // Restart the current process
  const restartProcess = async () => {
    setCompletedSteps([]);
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]));
    const routes = selectedOption === "drugDiscovery" ? drugDiscoverySteps : dashboardRoutes;
    const prompt = `Starting over with ${routes[0].name}.`;
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
    await speakResponse(prompt);
    await new Promise((resolve) => setTimeout(resolve, 1000)); // 1-second delay
    navigate(routes[0].path, { replace: true, state: { fromJarvis: true } });
    setTargetRoute(routes[0].path);
    setSubMode("explaining_tab");
    setCurrentToolIndex(0);
    lastProcessedRoute.current = routes[0].path;
    await speakFieldsForRoute(routes[0].path);
  };

  // Speak form fields and tab info based on route
  const speakFieldsForRoute = async (route) => {
    const normalizedRoute = normalizePath(route);
    let message = "";

    switch (normalizedRoute) {
      case normalizePath("/dashboard/protein-structure"):
        message = "You are in the Protein Structure tool. This tool visualizes protein structures in 3D, helping you analyze molecular configurations.";
        break;
      case normalizePath("/dashboard/protein-structure-mutation"):
        message = "You are in the Protein Structure Mutation tool. This tool generates new drug molecules by modifying protein structures.";
        break;
      case normalizePath("/dashboard/cost-estimation"):
        message = "You are in the Cost Estimation tool. This tool estimates synthesis and production costs for drug development.";
        break;
      case normalizePath("/dashboard/ai-research-paper-generator"):
        message = "You are in the AI Research Paper Generator tool. This tool generates research paper drafts based on your drug discovery data.";
        break;
      case normalizePath("/dashboard/ai-driven-target-prediction"):
        message = "You are in the AI-Driven Target Prediction tool. This tool predicts potential drug targets using artificial intelligence.";
        break;
      case normalizePath("/dashboard/ai-naming"):
        message = "You are in the AI Naming Suggestion tool. This tool generates creative and relevant drug names.";
        break;
      case normalizePath("/dashboard/summary"):
        message = "You are in the Summary tool. This tool provides a comprehensive review of your drug discovery project.";
        break;
      default:
        message = "It seems you're on an unrecognized page. Let's return to the first tool.";
        navigate(dashboardRoutes[0].path, { replace: true, state: { fromJarvis: true } });
        lastProcessedRoute.current = dashboardRoutes[0].path;
        break;
    }

    const fullMessage = `Now at ${getRouteName(route)}. ${message}`;
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: fullMessage }]);
    await speakResponse(fullMessage);

    if (mode === "beginner") {
      setSubMode("waiting_for_question");
      const questionPrompt = "Do you have any questions about this tool? If not, you can select another step or proceed to the next step.";
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: questionPrompt }]);
      await speakResponse(questionPrompt);
    }
  };

  // Get route name for display
  const getRouteName = (route) => {
    const foundRoute =
      dashboardRoutes.find((r) => normalizePath(r.path) === normalizePath(route)) ||
      drugDiscoverySteps.find((s) => normalizePath(s.path) === normalizePath(route));
    return foundRoute ? foundRoute.name : "Unknown";
  };

  // Scroll to bottom of conversation
  useEffect(() => {
    if (conversationRef.current) {
      conversationRef.current.scrollTop = conversationRef.current.scrollHeight;
    }
  }, [conversationHistory]);

  // Clean up on unmount
  useEffect(() => {
    return () => {
      isCleaningUp.current = true;
      window.speechSynthesis.cancel();
      if (recognitionRef.current) {
        try {
          recognitionRef.current.stop();
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`);
        }
      }
      recognitionRef.current = null;
    };
  }, []);

  // Save progress to localStorage when completedSteps changes
  useEffect(() => {
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify(completedSteps));
  }, [completedSteps]);

  return (
    <div className="font-sans">
      {/* Trigger Button */}
      <button
        onClick={() => setIsPanelOpen(true)}
        disabled={isPanelOpen}
        onMouseEnter={() => setTooltipVisible("trigger")}
        onMouseLeave={() => setTooltipVisible("")}
        aria-label={isPanelOpen ? "Close Jarvis Assistant" : "Open Jarvis Assistant"}
        className={`fixed bottom-6 right-6 h-16 w-16 rounded-full shadow-lg transition-all duration-200 transform hover:scale-105 flex items-center justify-center ${
          isPanelOpen
            ? "bg-gradient-to-r from-blue-500 to-blue-600 hover:from-blue-600 hover:to-blue-700"
            : "bg-gradient-to-br from-blue-500 to-blue-600 hover:from-blue-600 hover:to-blue-700"
        }`}
      >
        <Bot className="h-6 w-6 text-white" />
      </button>
      {tooltipVisible === "trigger" && (
        <div className="absolute bottom-20 right-4 px-2 py-1 bg-gray-800 text-white text-xs rounded-lg shadow-sm">
          {isPanelOpen ? "Close Jarvis" : "Ask Jarvis"}
        </div>
      )}

      {/* Pointing Message */}
      {!isPanelOpen && (
        <div className="fixed bottom-24 right-6">
          <div className="relative bg-white rounded-lg shadow-lg p-3 mb-2 animate-pulse-scale">
            <div className="absolute right-4 bottom-2 transform translate-y-full w-3 h-3 bg-white rotate-45"></div>
            <div className="flex items-center space-x-2">
              <div className="bg-blue-500/20 p-1.5 rounded-full">
                <Bot className="h-4 w-4 text-blue-600" />
              </div>
              <div>
                <div className="text-sm font-semibold text-gray-800">Drug Discovery Assistant</div>
                <div className="text-xs text-gray-500">Click to interact with Jarvis</div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Right Side Panel */}
      {isPanelOpen && (
        <div className="fixed top-0 right-0 h-full md:w-96 bg-white shadow-2xl z-50 flex flex-col border-l border-gray-200 animate-in slide-in-from-right duration-300">
          <div className="flex items-center justify-between p-4 border-b border-gray-200">
            <div className="flex items-center space-x-2">
              <button
                onClick={resetToInitialSelection}
                aria-label="Back to initial selection"
                className="rounded-full h-8 w-8 bg-gray-100 hover:bg-gray-200 text-gray-600 flex items-center justify-center"
              >
                <ArrowLeft className="h-4 w-4" />
              </button>
              <Dna className="h-6 w-6 text-blue-600" />
              <h2 className="text-xl font-semibold text-gray-800">Jarvis</h2>
            </div>
            <div className="flex items-center space-x-2">
              <button
                onClick={() => {
                  const helpPrompt =
                    "In beginner mode, you can type questions or follow guided steps. In doubt mode, use voice input to ask questions. Use the buttons to navigate or say 'stop' to end the conversation.";
                  setConversationHistory((prev) => [...prev, { type: "jarvis", text: helpPrompt }]);
                  speakResponse(helpPrompt);
                }}
                aria-label="Help"
                className="rounded-full h-8 w-8 bg-gray-100 hover:bg-gray-200 text-gray-600 flex items-center justify-center"
              >
                <HelpCircle className="h-4 w-4" />
              </button>
              <button
                onClick={stopConversation}
                aria-label="Close panel"
                className="rounded-full h-8 w-8 bg-gray-100 hover:bg-gray-200 text-gray-600 flex items-center justify-center"
              >
                <X className="h-4 w-4" />
              </button>
            </div>
          </div>

          {!mode && (
            <div className="flex-1 p-6 flex flex-col items-center justify-center space-y-6">
              <div className="text-center">
                <h3 className="text-lg font-semibold text-gray-800 mb-2">Welcome to Jarvis</h3>
                <p className="text-sm text-gray-600">Choose an option to start your drug discovery journey.</p>
              </div>
              <button
                onClick={() => startConversation("beginner")}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-all duration-200 text-sm font-semibold"
                aria-label="Start as a beginner"
              >
                Are you a beginner?
              </button>
              <button
                onClick={() => startConversation("doubt")}
                className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-all duration-200 text-sm font-semibold"
                aria-label="Ask a question"
              >
                Ask a question
              </button>
            </div>
          )}

          {mode && (
            <div className="flex flex-col h-full">
              {/* Header with progress bar for drug discovery */}
              {mode === "beginner" &&
                selectedOption === "drugDiscovery" &&
                subMode !== "waiting_for_selection" &&
                subMode !== "waiting_for_step_selection" &&
                subMode !== "process_ended" && (
                  <div className="p-4 bg-blue-50 border-b border-gray-200">
                    <div className="flex items-center justify-between mb-2">
                      <h3 className="text-sm font-medium text-blue-600">Drug Discovery Process</h3>
                      <div className="text-xs text-blue-600">
                        Step {currentToolIndex + 1} of {drugDiscoverySteps.length}
                      </div>
                    </div>
                    <div className="w-full bg-gray-200 border rounded-full h-2">
                      <div
                        className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                        style={{ width: `${((currentToolIndex + 1) / drugDiscoverySteps.length) * 100}%` }}
                      ></div>
                    </div>
                    <div className="flex items-center space-x-4 overflow-x-auto pb-2 mt-2">
                      {drugDiscoverySteps.map((step, index) => (
                        <div
                          key={index}
                          className={`flex-shrink-0 flex items-center ${
                            completedSteps.includes(step.path)
                              ? "text-green-600"
                              : index === currentToolIndex
                              ? "text-blue-600"
                              : "text-gray-400"
                          }`}
                        >
                          <div
                            className={`h-8 w-8 rounded-full flex items-center justify-center ${
                              completedSteps.includes(step.path)
                                ? "bg-green-100 border border-green-300"
                                : index === currentToolIndex
                                ? "bg-blue-100 border border-blue-500"
                                : "bg-gray-100"
                            }`}
                          >
                            {completedSteps.includes(step.path) ? <CheckCircle className="h-4 w-4" /> : index + 1}
                          </div>
                          {index < drugDiscoverySteps.length - 1 && <ChevronRight className="h-4 w-4 mx-1" />}
                        </div>
                      ))}
                    </div>
                  </div>
                )}

              {/* Current step details */}
              {mode === "beginner" &&
                selectedOption === "drugDiscovery" &&
                subMode !== "waiting_for_selection" &&
                subMode !== "waiting_for_step_selection" &&
                subMode !== "process_ended" && (
                  <div className="p-4 bg-blue-50 border-b border-gray-200">
                    <div className="flex items-start space-x-3">
                      <div className="bg-blue-100 p-2 rounded-full">
                        {iconMap[drugDiscoverySteps[currentToolIndex]?.icon]}
                      </div>
                      <div>
                        <h4 className="text-sm font-medium text-gray-800">
                          {drugDiscoverySteps[currentToolIndex]?.name}
                        </h4>
                        <p className="text-xs text-gray-600">{drugDiscoverySteps[currentToolIndex]?.description}</p>
                      </div>
                    </div>
                  </div>
                )}

              {/* Manual navigation buttons */}
              {mode === "beginner" &&
                selectedOption === "drugDiscovery" &&
                subMode !== "waiting_for_selection" &&
                subMode !== "waiting_for_step_selection" &&
                subMode !== "process_ended" && (
                  <div className="flex justify-between p-4 bg-white border-b border-gray-200">
                    <button
                      disabled={true} // Disable Previous button
                      className="flex items-center space-x-1 px-3 py-2 rounded-md text-sm font-medium text-gray-400 cursor-not-allowed"
                      aria-label="Previous step (disabled)"
                    >
                      <ChevronLeft className="h-4 w-4" />
                      <span>Previous</span>
                    </button>
                    {currentToolIndex === drugDiscoverySteps.length - 1 ? (
                      <button
                        onClick={endProcess}
                        disabled={isSpeaking}
                        className={`flex items-center space-x-1 px-3 py-2 rounded-md text-sm font-medium ${
                          isSpeaking ? "text-gray-400 cursor-not-allowed" : "text-blue-600 hover:bg-blue-50"
                        }`}
                        aria-label="End drug discovery process"
                      >
                        <span>End Process</span>
                        <ChevronRight className="h-4 w-4" />
                      </button>
                    ) : (
                      <button
                        onClick={proceedToNextStep}
                        disabled={isSpeaking}
                        className={`flex items-center space-x-1 px-3 py-2 rounded-md text-sm font-medium ${
                          isSpeaking ? "text-gray-400 cursor-not-allowed" : "text-blue-600 hover:bg-blue-50"
                        }`}
                        aria-label="Next step"
                      >
                        <span>Next</span>
                        <ChevronRight className="h-4 w-4" />
                      </button>
                    )}
                  </div>
                )}

              {/* Manual selection for beginner mode */}
              {mode === "beginner" && subMode === "waiting_for_selection" && (
                <div className="p-6 bg-blue-50 border-b border-gray-200">
                  <h3 className="text-sm font-medium text-gray-800 mb-3">Select an option to continue:</h3>
                  <div className="grid grid-cols-1 gap-3">
                    <button
                      onClick={handleDashboardSelection}
                      className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                      aria-label="Explore dashboard tools"
                    >
                      <div className="bg-blue-100 p-2 rounded-full">{iconMap["Atom"]}</div>
                      <div className="text-left">
                        <div className="text-sm font-medium text-gray-800">Dashboard Tools</div>
                        <div className="text-xs text-gray-500">Explore all available tools</div>
                      </div>
                    </button>
                    <button
                      onClick={handleDrugDiscoverySelection}
                      className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                      aria-label="Start drug discovery process"
                    >
                      <div className="bg-blue-100 p-2 rounded-full">{iconMap["Pill"]}</div>
                      <div className="text-left">
                        <div className="text-sm font-medium text-gray-800">Drug Discovery Process</div>
                        <div className="text-xs text-gray-500">Guided step-by-step journey</div>
                      </div>
                    </button>
                  </div>
                  <div className="flex justify-between gap-2 mt-4">
                    <button
                      onClick={stopConversation}
                      className="flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white bg-red-600 hover:bg-red-700 transition-all duration-200"
                      aria-label="End conversation"
                    >
                      End Conversation
                    </button>
                  </div>
                </div>
              )}

              {/* Tool selection for dashboard */}
              {mode === "beginner" && subMode === "waiting_for_tool_selection" && (
                <div className="p-6 bg-blue-50 border-b border-gray-200">
                  <h3 className="text-sm font-medium text-gray-800 mb-3">Select a tool to explore:</h3>
                  <div className="grid grid-cols-1 gap-3">
                    {dashboardRoutes.map((tool, index) => (
                      <button
                        key={index}
                        onClick={() => handleToolSelection(tool)}
                        className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                        aria-label={`Select ${tool.name} tool`}
                      >
                        <div className="bg-blue-100 p-2 rounded-full">{iconMap[tool.icon]}</div>
                        <div className="text-left">
                          <div className="text-sm font-medium text-gray-800">{tool.name}</div>
                        </div>
                      </button>
                    ))}
                  </div>
                  <div className="flex justify-between gap-2 mt-4">
                    <button
                      onClick={async () => {
                        const prompt = "Returning to option selection.";
                        setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
                        await speakResponse(prompt);
                        setSubMode("waiting_for_selection");
                        setSelectedOption(null);
                      }}
                      className="flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white bg-blue-600 hover:bg-blue-700 transition-all duration-200"
                      aria-label="Go back to option selection"
                    >
                      Back
                    </button>
                    <button
                      onClick={stopConversation}
                      className="flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white bg-red-600 hover:bg-red-700 transition-all duration-200"
                      aria-label="End conversation"
                    >
                      End Conversation
                    </button>
                  </div>
                </div>
              )}

              {/* Step selection for drug discovery */}
              {mode === "beginner" && subMode === "waiting_for_step_selection" && (
                <div className="p-6 bg-blue-50 border-b border-gray-200">
                  <h3 className="text-sm font-medium text-gray-800 mb-3">Select a step to explore:</h3>
                  <div className="grid grid-cols-1 gap-3">
                    {drugDiscoverySteps.map((step, index) => {
                      const isCurrentStep = index === currentToolIndex;
                      const isCompleted = completedSteps.includes(step.path);
                      const isDisabled = !isCurrentStep || isCompleted;
                      return (
                        <button
                          key={index}
                          onClick={() => handleStepSelection(step)}
                          disabled={isDisabled}
                          className={`flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 transition-all duration-200 ${
                            isDisabled ? "opacity-50 cursor-not-allowed" : "hover:border-blue-300 hover:bg-blue-50"
                          }`}
                          aria-label={isDisabled ? `${step.name} (locked)` : `Select ${step.name}`}
                        >
                          <div className="bg-blue-100 p-2 rounded-full">
                            {isCompleted ? <CheckCircle className="h-4 w-4 text-green-600" /> : iconMap[step.icon]}
                          </div>
                          <div className="text-left">
                            <div className="text-sm font-medium text-gray-800">{step.name}</div>
                            <div className="text-xs text-gray-500">{step.description}</div>
                          </div>
                        </button>
                      );
                    })}
                  </div>
                  <div className="flex justify-between gap-2 mt-4">
                    <button
                      onClick={async () => {
                        const prompt = "Returning to option selection.";
                        setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
                        await speakResponse(prompt);
                        setSubMode("waiting_for_selection");
                        setSelectedOption(null);
                      }}
                      className="flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white bg-blue-600 hover:bg-blue-700 transition-all duration-200"
                      aria-label="Go back to option selection"
                    >
                      Back
                    </button>
                    <button
                      onClick={stopConversation}
                      className="flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white bg-red-600 hover:bg-red-700 transition-all duration-200"
                      aria-label="End conversation"
                    >
                      End Conversation
                    </button>
                  </div>
                </div>
              )}

              {/* Process ended options */}
              {mode === "beginner" && subMode === "process_ended" && (
                <div className="p-6 bg-blue-50 border-b border-gray-200">
                  <h3 className="text-sm font-medium text-gray-800 mb-3">Drug Discovery Completed!</h3>
                  <div className="grid grid-cols-1 gap-3">
                    <button
                      onClick={restartProcess}
                      className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                      aria-label="Restart drug discovery process"
                    >
                      <div className="bg-blue-100 p-2 rounded-full">{iconMap["Pill"]}</div>
                      <div className="text-left">
                        <div className="text-sm font-medium text-gray-800">Restart Drug Discovery</div>
                        <div className="text-xs text-gray-500">Start the process again</div>
                      </div>
                    </button>
                    <button
                      onClick={handleDashboardSelection}
                      className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                      aria-label="Explore dashboard tools"
                    >
                      <div className="bg-blue-100 p-2 rounded-full">{iconMap["Atom"]}</div>
                      <div className="text-left">
                        <div className="text-sm font-medium text-gray-800">Explore Dashboard Tools</div>
                        <div className="text-xs text-gray-500">Try other tools</div>
                      </div>
                    </button>
                    <button
                      onClick={() => {
                        setSubMode("waiting_for_selection");
                        const prompt = "Returning to option selection.";
                        setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
                        speakResponse(prompt);
                      }}
                      className="flex items-center space-x-3 p-3 bg-white rounded-lg border border-gray-200 hover:border-blue-300 hover:bg-blue-50 transition-all duration-200"
                      aria-label="Return to option selection"
                    >
                      <div className="bg-blue-100 p-2 rounded-full">{iconMap["Bot"]}</div>
                      <div className="text-left">
                        <div className="text-sm font-medium text-gray-800">Back to Options</div>
                        <div className="text-xs text-gray-500">Choose between dashboard or drug discovery</div>
                      </div>
                    </button>
                  </div>
                </div>
              )}

              {/* Conversation History */}
              <div ref={conversationRef} className="flex-auto overflow-y-auto p-6 space-y-4">
                {conversationHistory.map((msg, index) => (
                  <div key={index} className={`flex ${msg.type === "user" ? "justify-end" : "justify-start"}`}>
                    <div
                      className={`max-w-[80%] ${
                        msg.type === "user"
                          ? "bg-blue-600 text-white rounded-tr-none"
                          : "bg-white text-gray-800 rounded-tl-none border border-gray-200 shadow-sm"
                      } rounded-2xl px-4 py-2 transition-all duration-200`}
                    >
                      {msg.type === "jarvis" && (
                        <div className="flex items-center space-x-2 mb-1 pb-1 border-b border-gray-200">
                          <Bot className="h-3 w-3 text-blue-600" />
                          <span className="text-xs font-medium text-blue-600">JARVIS</span>
                        </div>
                      )}
                      <p className="text-sm whitespace-pre-wrap">{msg.text}</p>
                    </div>
                  </div>
                ))}
                {(isLoading || isSpeaking) && (
                  <div className="flex justify-start">
                    <div className="max-w-[80%] w-full py-2 px-4 bg-white border rounded-2xl border-gray-200 shadow-sm">
                      <div className="flex items-center space-x-2 pb-0">
                        <Bot className="h-3 w-3 text-blue-600" />
                        <span className="text-xs font-semibold text-blue-600">JARVIS</span>
                      </div>
                      <div className="flex space-x-1 items-center">
                        <div className="h-3 w-1 bg-blue-400 animate-wave"></div>
                        <div className="h-4 w-1 bg-blue-400 animate-wave" style={{ animationDelay: "100ms" }}></div>
                        <div className="h-5 w-1 bg-blue-400 animate-wave" style={{ animationDelay: "200ms" }}></div>
                        <div className="h-4 w-1 bg-blue-400 animate-wave" style={{ animationDelay: "300ms" }}></div>
                        <div className="h-3 w-1 bg-blue-400 animate-wave" style={{ animationDelay: "400ms" }}></div>
                      </div>
                    </div>
                  </div>
                )}
              </div>

              {/* Input Area */}
              <div className="p-6 border-t border-gray-200">
                {mode === "doubt" && (
                  <div className="flex items-center justify-between px-4 py-2 rounded-full bg-white border border-gray-200 shadow-sm">
                    <div className="flex items-center space-x-2 text-sm text-gray-600">
                      {isListening && mode === "doubt" ? (
                        <>
                          <div className="relative">
                            <div className="absolute inset-0 bg-green-100 rounded-full animate-ping"></div>
                            <Mic className="h-5 w-5 text-green-600" />
                          </div>
                          <span className="text-sm">Listening...</span>
                        </>
                      ) : (
                        <>
                          <MicOff className="h-5 w-5 text-red-600" />
                          <span className="text-sm">Voice recognition paused</span>
                        </>
                      )}
                    </div>
                    <button
                      onClick={() => (isListening ? stopConversation() : startConversation("doubt"))}
                      className={`px-4 py-2 rounded-lg text-sm font-semibold text-white ${
                        isListening ? "bg-red-600 hover:bg-red-700" : "bg-blue-600 hover:bg-blue-700"
                      } transition-all duration-200`}
                      aria-label={isListening ? "Stop voice conversation" : "Start voice conversation"}
                    >
                      {isListening ? "Stop" : "Start"}
                    </button>
                  </div>
                )}
                {mode === "beginner" && subMode === "waiting_for_question" && (
                  <div className="flex flex-col space-y-3">
                    <input
                      type="text"
                      placeholder="Ask a question about this tool..."
                      className="w-full px-4 py-2 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-blue-500 text-sm"
                      onKeyDown={async (e) => {
                        if (e.key === "Enter" && e.target.value.trim()) {
                          setSubMode("resolving_doubt");
                          const query = e.target.value.trim();
                          setConversationHistory((prev) => [...prev, { type: "user", text: query }]);
                          await fetchResponse(query);
                          e.target.value = "";
                        }
                      }}
                      aria-label="Ask a question about the current tool"
                    />
                    {suggestedQuestions[targetRoute]?.length > 0 && (
                      <div className="flex flex-wrap gap-2">
                        {suggestedQuestions[targetRoute].map((q, index) => (
                          <button
                            key={index}
                            onClick={() => {
                              setSubMode("resolving_doubt");
                              setConversationHistory((prev) => [...prev, { type: "user", text: q }]);
                              fetchResponse(q);
                            }}
                            className="px-3 py-1 bg-blue-100 text-blue-600 rounded-full text-xs hover:bg-blue-200 transition-all duration-200"
                            aria-label={`Ask suggested question: ${q}`}
                          >
                            {q}
                          </button>
                        ))}
                      </div>
                    )}
                    <div className="flex justify-end gap-2">
                      {selectedOption === "drugDiscovery" && currentToolIndex < drugDiscoverySteps.length - 1 && (
                        <button
                          onClick={proceedToNextStep}
                          disabled={isSpeaking}
                          className={`flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white ${
                            isSpeaking ? "bg-gray-400 cursor-not-allowed" : "bg-blue-600 hover:bg-blue-700"
                          } transition-all duration-200`}
                          aria-label="Proceed to next step"
                        >
                          Next Step
                        </button>
                      )}
                      {/* <button
                        onClick={async () => {
                          const prompt = `Returning to ${selectedOption === "drugDiscovery" ? "step" : "tool"} selection.`;
                          setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }]);
                          await speakResponse(prompt);
                          setSubMode(selectedOption === "drugDiscovery" ? "waiting_for_step_selection" : "waiting_for_tool_selection");
                        }}
                        disabled={isSpeaking}
                        className={`flex-1 px-4 py-2 rounded-lg text-sm font-semibold text-white ${
                          isSpeaking ? "bg-gray-400 cursor-not-allowed" : "bg-blue-600 hover:bg-blue-700"
                        } transition-all duration-200`}
                        aria-label={`Select another ${selectedOption === "drugDiscovery" ? "step" : "tool"}`}
                      >
                        Select Another {selectedOption === "drugDiscovery" ? "Step" : "Tool"}
                      </button> */}
                    </div>
                  </div>
                )}
                <div className="mt-2 text-xs text-center text-gray-500">
                  Jarvis can assist with molecular structures, drug discovery processes, and medical queries
                </div>
              </div>
            </div>
          )}
        </div>
      )}
      <style>
        {`
          @keyframes pulse-scale {
            0%, 100% { transform: scale(1); }
            50% { transform: scale(1.05); }
          }
          @keyframes wave {
            0%, 100% { transform: translateY(0); }
            50% { transform: translateY(-4px); }
          }
          .animate-wave {
            animation: wave 0.8s infinite;
          }
          .animate-ping {
            animation: ping 1s cubic-bezier(0, 0, 0.2, 1) infinite;
          }
          @keyframes ping {
            75%, 100% { transform: scale(2); opacity: 0; }
          }
        `}
      </style>
    </div>
  );
}