import { useState, useEffect, useRef } from "react"
import { Bot, Microscope, Atom, Brain, X, Mic, MicOff, Volume2, VolumeX, Dna, Pill, FlaskConical } from "lucide-react"
import { useAuthStore } from "../../Store/auth.store.js"
export default function JarvisVoiceAssistant() {
  const [isListening, setIsListening] = useState(false)
  const [isPopupOpen, setIsPopupOpen] = useState(false)
  const [transcript, setTranscript] = useState("")
  const [response, setResponse] = useState("")
  const [isSpeaking, setIsSpeaking] = useState(false)
  const [isLoading, setIsLoading] = useState(false)
  const [conversationHistory, setConversationHistory] = useState([])
  const [showTooltip, setShowTooltip] = useState("")
  const recognitionRef = useRef(null)
  const isCleaningUp = useRef(false)
  const conversationRef = useRef(null)
  const { user } = useAuthStore()

  // Initialize Web Speech API
  const createRecognition = () => {
    const SpeechRecognition = window.SpeechRecognition || (window).webkitSpeechRecognition
    const recognition = new SpeechRecognition()
    recognition.lang = "en-US"
    recognition.interimResults = false
    recognition.maxAlternatives = 1
    recognition.continuous = true
    return recognition
  }

  // Configure recognition handlers
  const setupRecognition = (recognition) => {
    recognition.onresult = (event) => {
      const speechResult = event.results[event.results.length - 1][0].transcript.trim().toLowerCase()
      console.log("Captured transcript:", speechResult)
      setTranscript(speechResult)

      // Add to conversation history
      setConversationHistory((prev) => [...prev, { type: "user", text: speechResult }])

      // Check for stop commands
      if (speechResult.includes("stop") || speechResult.includes("thank you")) {
        stopConversation()
        return
      }

      fetchResponse(speechResult)
    }

    recognition.onend = () => {
      console.log("Recognition ended, isListening:", isListening, "isCleaningUp:", isCleaningUp.current)
      if (isListening && !isCleaningUp.current) {
        setTimeout(() => {
          try {
            recognition.start()
            console.log("Speech recognition restarted")
          } catch (error) {
            console.error("Error restarting recognition:", error)
            stopConversation()
          }
        }, 500)
      }
    }

    recognition.onerror = (event) => {
      console.error("Speech recognition error:", event.error)
      if (isListening && !isCleaningUp.current) {
        setTimeout(() => {
          try {
            recognition.start()
            console.log("Retrying speech recognition")
          } catch (error) {
            console.error("Retry failed:", error)
            stopConversation()
          }
        }, 500)
      }
    }
  }

  // Start conversation
  const startConversation = () => {
    if (isListening || isPopupOpen) {
      console.log("Conversation already active, ignoring start")
      return
    }
    console.log("Starting conversation")
    isCleaningUp.current = false
    setIsPopupOpen(true)
    setIsListening(true)
    setTranscript("")
    setResponse("")
    setConversationHistory([])

    // Create new recognition instance
    recognitionRef.current = createRecognition()
    setupRecognition(recognitionRef.current)

    const welcomeMessage =
      "I am Jarvis, I am here to assist you in your new drug discovery journey. How can I help you?"
    setConversationHistory([{ type: "jarvis", text: welcomeMessage }])
    speakResponse(welcomeMessage)

    setTimeout(() => {
      try {
        recognitionRef.current.start()
        console.log("Speech recognition started")
      } catch (error) {
        console.error("Error starting recognition:", error)
        stopConversation()
      }
    }, 100)
  }

  // Stop conversation
  const stopConversation = () => {
    console.log("Stopping conversation")
    isCleaningUp.current = true
    setIsListening(false)

    const goodbyeMessage = "Goodbye. Thank you for using Jarvis."
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: goodbyeMessage }])
    speakResponse(goodbyeMessage)

    if (recognitionRef.current) {
      recognitionRef.current.stop()
      recognitionRef.current = null
    }

    // Close the popup after the goodbye message is spoken
    setTimeout(() => {
      setIsPopupOpen(false)
    }, 2000)
  }

  // Toggle speech
  const toggleSpeech = () => {
    if (isSpeaking) {
      window.speechSynthesis.cancel()
      setIsSpeaking(false)
    } else if (response) {
      speakResponse(response)
    }
  }

  // Fetch response from Gemini API
  const fetchResponse = async (query) => {
    console.log("Fetching response for query:", query)
    setIsLoading(true)

    const problemStatement = `The process of drug discovery is time-consuming, expensive, and often inefficient, with a high rate of failure in clinical trials. Traditional methods rely heavily on trial and error, requiring years of research and significant financial investment. Additionally, the complexity of biological systems and the vast chemical space make it challenging to identify promising drug candidates efficiently. Generative AI, with its ability to analyze large datasets, predict molecular interactions, and generate novel compounds, has the potential to revolutionize this process. However, there is a lack of accessible, user-friendly tools that leverage generative AI to assist researchers in accelerating drug discovery while reducing costs and improving success rates.`

    const prompt = `You are Jarvis, an AI assistant specialized in drug discovery. You are trained on the following problem statement: "${problemStatement}". Answer the user's query related to drug discovery, from basic molecule or drug names to complex information about the process, generative AI applications, or challenges. Provide a concise, accurate response suitable for researchers. Query: ${query}`

    try {
      const response = await fetch(
        "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=AIzaSyDyujm50dHMYvn1V50dDDqcAhgUqCOuUGU",
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            contents: [{ parts: [{ text: prompt }] }],
          }),
        },
      )

      if (!response.ok) {
        throw new Error(`HTTP error! Status: ${response.status}`)
      }

      const data = await response.json()
      const text = data.candidates[0].content.parts[0].text
      console.log("API response:", text)
      setResponse(text)
      setConversationHistory((prev) => [...prev, { type: "jarvis", text }])
      speakResponse(text)

      // Store in localStorage
      const jarvisEntry = {
        title: `Drug Discovery Query: ${query.slice(0, 50)}`,
        questions: query,
        answers: text,
        userId: user._id, // Hardcoded; replace with auth in production
        createdAt: new Date().toISOString(),
      }
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]")
      storedEntries.push(jarvisEntry)
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries))
      console.log("Stored entry in localStorage:", jarvisEntry)
    } catch (error) {
      console.error("Error fetching Gemini API:", error)
      const errorMessage = "Sorry, I encountered an API error. Please try again."
      setResponse(errorMessage)
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }])
      speakResponse(errorMessage)

      // Store error in localStorage
      const jarvisEntry = {
        title: `Drug Discovery Query: ${query.slice(0, 50)}`,
        questions: query,
        answers: "Error: Unable to fetch response",
        userId: "66f8b1e4c7e1b2f5a4b3c2d1",
        createdAt: new Date().toISOString(),
      }
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]")
      storedEntries.push(jarvisEntry)
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries))
      console.log("Stored error entry in localStorage:", jarvisEntry)
    } finally {
      setIsLoading(false)
    }
  }

  // Speak the response using Web Speech API
  const speakResponse = (text) => {
    console.log("Speaking response:", text)
    window.speechSynthesis.cancel()
    const utterance = new SpeechSynthesisUtterance(text)
    utterance.lang = "en-US"
    utterance.onstart = () => setIsSpeaking(true)
    utterance.onend = () => {
      setIsSpeaking(false)
      if (text === "Goodbye.") {
        setIsListening(false)
        setIsPopupOpen(false)
        isCleaningUp.current = true
      }
    }
    window.speechSynthesis.speak(utterance)
  }

  // Scroll to bottom of conversation
  useEffect(() => {
    if (conversationRef.current) {
      conversationRef.current.scrollTop = conversationRef.current.scrollHeight
    }
  }, [conversationHistory])

  // Clean up speech synthesis and recognition on component unmount
  useEffect(() => {
    return () => {
      console.log("Cleaning up speech synthesis and recognition")
      isCleaningUp.current = true
      window.speechSynthesis.cancel()
      if (recognitionRef.current) {
        recognitionRef.current.stop()
        recognitionRef.current = null
      }
    }
  }, [])

  // return (
  //   <div className="font-sans">
  //     {/* Trigger Button */}
  //     <div className="relative">
  //       <button
  //         onClick={startConversation}
  //         disabled={isListening}
  //         onMouseEnter={() => setShowTooltip("trigger")}
  //         onMouseLeave={() => setShowTooltip("")}
  //         className={`fixed bottom-6 right-6 h-16 w-16 rounded-full shadow-lg transition-all transform hover:scale-105 flex items-center justify-center ${
  //           isListening
  //             ? "bg-gradient-to-r from-red-500 to-red-600 hover:from-red-600 hover:to-red-700"
  //             : "bg-gradient-to-r from-cyan-600 to-blue-700 hover:from-cyan-700 hover:to-blue-800"
  //         }`}
  //       >
  //         {isListening ? <MicOff className="h-6 w-6 text-white" /> : <Mic className="h-6 w-6 text-white" />}
  //       </button>
  //       {showTooltip === "trigger" && (
  //         <div className="absolute bottom-full right-0 mb-2 px-2 py-1 bg-gray-800 text-white text-xs rounded shadow-lg whitespace-nowrap">
  //           {isListening ? "Stop Listening" : "Ask Jarvis"}
  //         </div>
  //       )}
  //     </div>

  //     {/* Popup Modal */}
  //     {isPopupOpen && (
  //       <div className="fixed inset-0 bg-black/70 flex items-center justify-center z-50 p-4 backdrop-blur-sm animate-in fade-in duration-300">
  //         <div className="w-full max-w-5xl overflow-hidden rounded-xl border-0 shadow-2xl bg-gradient-to-br from-slate-900 to-slate-800 text-white">
  //           <div className="absolute top-4 right-4">
  //             <button
  //               onClick={stopConversation}
  //               className="rounded-full h-8 w-8 bg-white/10 hover:bg-white/20 text-white flex items-center justify-center"
  //             >
  //               <X className="h-4 w-4" />
  //               <span className="sr-only">Close</span>
  //             </button>
  //           </div>

  //           <div className="flex flex-col md:flex-row h-[80vh] max-h-[800px]">
  //             {/* Left side: Jarvis visualization */}
  //             <div className="w-full md:w-2/5 bg-gradient-to-br from-blue-900/50 to-indigo-900/50 flex flex-col items-center justify-between p-6 border-r border-white/10">
  //               <div className="flex flex-col items-center text-center space-y-2">
  //                 <div className="relative">
  //                   <div
  //                     className="absolute inset-0 rounded-full bg-cyan-500/20 animate-pulse"
  //                     style={{ animationDuration: "3s" }}
  //                   ></div>
  //                   <div className="h-32 w-32 rounded-full bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center">
  //                     <Dna className="h-16 w-16 text-white" />
  //                   </div>
  //                 </div>

  //                 <h1 className="text-3xl font-bold mt-4 bg-clip-text text-transparent bg-gradient-to-r from-cyan-400 to-blue-500">
  //                   JARVIS
  //                 </h1>
  //                 <p className="text-blue-200 text-sm">Drug Discovery Assistant</p>

  //                 <div className="flex items-center mt-2 space-x-1">
  //                   <div
  //                     className={`h-2 w-2 rounded-full ${isListening ? "bg-green-500 animate-pulse" : "bg-red-500"}`}
  //                   ></div>
  //                   <span className="text-xs text-blue-200">{isListening ? "Listening" : "Idle"}</span>
  //                 </div>
  //               </div>

  //               <div className="w-full ">
  //                 <img
  //                   src="../../../public/jarvis.jpg"
  //                   alt="Molecular AI Visualization"
  //                   className="rounded-xl shadow-lg max-w-full h-auto mt-1 object-cover border border-white/20"
  //                   onError={(e) => {
  //                     const target = e.target
  //                     target.onerror = null
  //                     target.src = "../../../public/jarvis.jpg"
  //                   }}
  //                 />
  //               </div>

  //               <div className="w-full space-y-4 mt-4">
  //                 <div className="grid grid-cols-4 gap-2">
  //                   <div className="flex items-center justify-center gap-1 py-1 bg-white/5 border border-white/10 rounded-md px-2">
  //                     <Atom className="h-3 w-3" /> <span className="text-xs">Molecules</span>
  //                   </div>
  //                   <div className="flex items-center justify-center gap-1 py-1 bg-white/5 border border-white/10 rounded-md px-2">
  //                     <Brain className="h-3 w-3" /> <span className="text-xs">AI</span>
  //                   </div>
  //                   <div className="flex items-center justify-center gap-1 py-1 bg-white/5 border border-white/10 rounded-md px-2">
  //                     <Microscope className="h-3 w-3" /> <span className="text-xs">Research</span>
  //                   </div>
  //                   <div className="flex items-center justify-center gap-1 py-1 bg-white/5 border border-white/10 rounded-md px-2">
  //                     <Pill className="h-3 w-3" /> <span className="text-xs">Drugs</span>
  //                   </div>
  //                 </div>

  //                 <div className="text-xs text-blue-200/70 text-center">Powered by Gemini 1.5 Flash</div>
  //               </div>
  //             </div>

  //             {/* Right side: Conversation */}
  //             <div className="w-full md:w-3/5 flex flex-col p-6">
  //               <div className="flex items-center justify-between mb-4">
  //                 <div className="flex items-center space-x-2">
  //                   <FlaskConical className="h-5 w-5 text-cyan-400" />
  //                   <h2 className="text-xl font-semibold">Drug Discovery Assistant</h2>
  //                 </div>

  //                 <div className="flex space-x-2">
  //                   <div className="relative">
  //                     <button
  //                       onClick={toggleSpeech}
  //                       onMouseEnter={() => setShowTooltip("speech")}
  //                       onMouseLeave={() => setShowTooltip("")}
  //                       className="h-8 w-8 rounded-full bg-white/10 hover:bg-white/20 flex items-center justify-center"
  //                     >
  //                       {isSpeaking ? <VolumeX className="h-4 w-4" /> : <Volume2 className="h-4 w-4" />}
  //                     </button>
  //                     {showTooltip === "speech" && (
  //                       <div className="absolute bottom-full right-0 mb-2 px-2 py-1 bg-gray-800 text-white text-xs rounded shadow-lg whitespace-nowrap">
  //                         {isSpeaking ? "Mute" : "Speak"}
  //                       </div>
  //                     )}
  //                   </div>
  //                 </div>
  //               </div>

  //               {/* Conversation history */}
  //               <div ref={conversationRef} className="flex-1 overflow-y-auto pr-2 space-y-4 mb-4 custom-scrollbar">
  //                 {conversationHistory.map((message, index) => (
  //                   <div key={index} className={`flex ${message.type === "user" ? "justify-end" : "justify-start"}`}>
  //                     <div
  //                       className={`max-w-[80%] rounded-2xl px-4 py-2 ${
  //                         message.type === "user"
  //                           ? "bg-blue-600 text-white rounded-tr-none"
  //                           : "bg-slate-700 text-white rounded-tl-none"
  //                       }`}
  //                     >
  //                       {message.type === "jarvis" && (
  //                         <div className="flex items-center space-x-2 mb-1 pb-1 border-b border-white/10">
  //                           <Bot className="h-3 w-3 text-cyan-400" />
  //                           <span className="text-xs font-medium text-cyan-400">JARVIS</span>
  //                         </div>
  //                       )}
  //                       <p className="text-sm whitespace-pre-wrap">{message.text}</p>
  //                     </div>
  //                   </div>
  //                 ))}

  //                 {isLoading && (
  //                   <div className="flex justify-start">
  //                     <div className="max-w-[80%] rounded-2xl px-4 py-2 bg-slate-700 text-white rounded-tl-none">
  //                       <div className="flex items-center space-x-2 mb-1 pb-1 border-b border-white/10">
  //                         <Bot className="h-3 w-3 text-cyan-400" />
  //                         <span className="text-xs font-medium text-cyan-400">JARVIS</span>
  //                       </div>
  //                       <div className="flex space-x-1 items-center">
  //                         <div
  //                           className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
  //                           style={{ animationDelay: "0ms" }}
  //                         ></div>
  //                         <div
  //                           className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
  //                           style={{ animationDelay: "150ms" }}
  //                         ></div>
  //                         <div
  //                           className="h-2 w-2 bg-blue-400 rounded-full animate-bounce"
  //                           style={{ animationDelay: "300ms" }}
  //                         ></div>
  //                       </div>
  //                     </div>
  //                   </div>
  //                 )}
  //               </div>

  //               {/* Input area */}
  //               <div className="mt-auto">
  //                 <div className="relative">
  //                   <div className="flex items-center justify-between p-4 rounded-xl bg-slate-700/50 border border-white/10">
  //                     <div className="flex items-center space-x-2 text-sm text-blue-200">
  //                       {isListening ? (
  //                         <>
  //                           <div className="relative">
  //                             <div className="absolute inset-0 bg-green-500/20 rounded-full animate-ping"></div>
  //                             <Mic className="h-5 w-5 text-green-400" />
  //                           </div>
  //                           <span>Listening... Say "stop" or "thank you" to end</span>
  //                         </>
  //                       ) : (
  //                         <>
  //                           <MicOff className="h-5 w-5 text-red-400" />
  //                           <span>Voice recognition paused</span>
  //                         </>
  //                       )}
  //                     </div>

  //                     <button
  //                       onClick={isListening ? stopConversation : startConversation}
  //                       className={`px-4 py-2 rounded-lg text-sm font-medium text-white ${
  //                         isListening ? "bg-red-600 hover:bg-red-700" : "bg-cyan-600 hover:bg-cyan-700"
  //                       } transition-colors`}
  //                     >
  //                       {isListening ? "Stop" : "Start"}
  //                     </button>
  //                   </div>
  //                 </div>

  //                 <div className="mt-2 text-xs text-center text-blue-200/70">
  //                   Jarvis can assist with molecular structures, drug discovery processes, and AI applications in
  //                   pharmaceutical research
  //                 </div>
  //               </div>
  //             </div>
  //           </div>
  //         </div>
  //       </div>
  //     )}


  //   </div>
  // )


  return (
    <div className="font-sans">
      {/* Trigger Button */}
      <div className="relative">
        <button
          onClick={startConversation}
          disabled={isListening}
          onMouseEnter={() => setShowTooltip("trigger")}
          onMouseLeave={() => setShowTooltip("")}
          className={`fixed bottom-6 right-6 h-16 w-16 rounded-full shadow-lg transition-all transform hover:scale-105 flex items-center justify-center ${isListening
              ? "bg-gradient-to-r from-red-500 to-red-600 hover:from-red-600 hover:to-red-700"
              : "bg-gradient-to-r from-blue-500 to-cyan-500 hover:from-blue-600 hover:to-cyan-600"
            }`}
        >
          {isListening ? <MicOff className="h-6 w-6 text-white" /> : <Mic className="h-6 w-6 text-white" />}
        </button>
        {showTooltip === "trigger" && (
          <div className="absolute bottom-full right-0 mb-2 px-2 py-1 bg-gray-800 text-white text-xs rounded shadow-lg whitespace-nowrap">
            {isListening ? "Stop Listening" : "Ask Jarvis"}
          </div>
        )}
      </div>
      {/* Pointing Message - hidden when popup is open */}
      {!isPopupOpen && (
        <div className="fixed bottom-24 right-7">
          <div
            className="relative bg-white rounded-lg shadow-lg p-3 mb-2"
            style={{
              animation: 'pulse-scale 2s infinite ease-in-out',
              transformOrigin: 'center'
            }}
          >
            <style>
              {`
          @keyframes pulse-scale {
            0%, 100% {
              transform: scale(1);
            }
            50% {
              transform: scale(1.05);
            }
          }
        `}
            </style>
            <div className="absolute right-4 bottom-2 transform translate-y-full w-3 h-3 bg-white rotate-45"></div>
            {/* <div className="absolute -rig transform -translate-x-1/2 translate-y-1/2 w-3 h-3 bg-white rotate-45"></div> */}
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

      {/* Popup Modal */}
      {isPopupOpen && (
        <div className="fixed inset-0 bg-black/30 flex items-center justify-center z-50 p-4 backdrop-blur-sm animate-in fade-in duration-300">
          <div className="w-full max-w-5xl overflow-hidden rounded-xl border border-gray-200 shadow-2xl bg-white text-gray-800">
            <div className="absolute top-4 right-4">
              <button
                onClick={stopConversation}
                className="rounded-full h-8 w-8 bg-gray-100 hover:bg-gray-200 text-gray-600 flex items-center justify-center"
              >
                <X className="h-4 w-4" />
                <span className="sr-only">Close</span>
              </button>
            </div>

            <div className="flex flex-col md:flex-row h-[80vh] max-h-[800px]">
              {/* Left side: Jarvis visualization */}
              <div className="w-full md:w-2/5 bg-gradient-to-br from-blue-50 to-cyan-50 flex flex-col items-center justify-between p-6 border-r border-gray-200">
                <div className="flex flex-col items-center text-center space-y-2">
                  <div className="relative">
                    <div
                      className="absolute inset-0 rounded-full bg-blue-500/10 animate-pulse"
                      style={{ animationDuration: "3s" }}
                    ></div>
                    <div className="h-32 w-32 rounded-full bg-gradient-to-br from-blue-400 to-cyan-500 flex items-center justify-center">
                      <Dna className="h-16 w-16 text-white" />
                    </div>
                  </div>

                  <h1 className="text-3xl font-bold mt-4 bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-cyan-600">
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

                <div className="w-64">
                  <img
                    src="../../../public/jarvis.jpg"
                    alt="Molecular AI Visualization"
                    className="rounded-xl shadow-md max-w-full h-auto object-cover border border-gray-200"
                    onError={(e) => {
                      const target = e.target
                      target.onerror = null
                      target.src = "../../../public/jarvis.jpg"
                    }}
                  />
                </div>

                <div className="w-full space-y-4 mt-4">
                  <div className="grid grid-cols-4 gap-2">
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Atom className="h-3 w-3 text-blue-600" />{" "}
                      <span className="text-xs text-gray-700">Molecules</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Brain className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">AI</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Microscope className="h-3 w-3 text-blue-600" />{" "}
                      <span className="text-xs text-gray-700">Research</span>
                    </div>
                    <div className="flex items-center justify-center gap-1 py-1 bg-white border border-gray-200 rounded-md px-2 shadow-sm">
                      <Pill className="h-3 w-3 text-blue-600" /> <span className="text-xs text-gray-700">Drugs</span>
                    </div>
                  </div>

                  <div className="text-xs text-gray-500 text-center">Powered by Gemini</div>
                </div>
              </div>

              {/* Right side: Conversation */}
              <div className="w-full md:w-3/5 flex flex-col p-6 bg-gray-50">
                <div className="flex items-center justify-between mb-4">
                  <div className="flex items-center space-x-2">
                    <FlaskConical className="h-5 w-5 text-blue-600" />
                    <h2 className="text-xl font-semibold text-gray-800">Drug Discovery Assistant</h2>
                  </div>

                  <div className="flex space-x-2">
                    <div className="relative">
                      <button
                        onClick={toggleSpeech}
                        onMouseEnter={() => setShowTooltip("speech")}
                        onMouseLeave={() => setShowTooltip("")}
                        className="h-8 w-8 rounded-full bg-gray-100 hover:bg-gray-200 flex items-center justify-center text-gray-600"
                      >
                        {isSpeaking ? <VolumeX className="h-4 w-4" /> : <Volume2 className="h-4 w-4" />}
                      </button>
                      {showTooltip === "speech" && (
                        <div className="absolute bottom-full right-0 mb-2 px-2 py-1 bg-gray-800 text-white text-xs rounded shadow-lg whitespace-nowrap">
                          {isSpeaking ? "Mute" : "Speak"}
                        </div>
                      )}
                    </div>
                  </div>
                </div>

                {/* Conversation history */}
                <div ref={conversationRef} className="flex-1 overflow-y-auto pr-2 space-y-4 mb-4 custom-scrollbar">
                  {conversationHistory.map((message, index) => (
                    <div key={index} className={`flex ${message.type === "user" ? "justify-end" : "justify-start"}`}>
                      <div
                        className={`max-w-[80%] rounded-2xl px-4 py-2 ${message.type === "user"
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
                      <div className="max-w-[80%] rounded-2xl px-4 py-2 bg-white text-gray-800 rounded-tl-none border border-gray-200 shadow-sm">
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

                {/* Input area */}
                <div className="mt-auto">
                  <div className="relative">
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

                      <button
                        onClick={isListening ? stopConversation : startConversation}
                        className={`px-4 py-2 rounded-lg text-sm font-medium text-white ${isListening ? "bg-red-500 hover:bg-red-600" : "bg-blue-500 hover:bg-blue-600"
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
            </div>
          </div>
        </div>
      )}


    </div>
  )
}
