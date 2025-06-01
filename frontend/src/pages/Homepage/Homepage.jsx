"use client"

import { useState, useEffect, useRef } from "react"
import {
  ArrowRight,
  Atom,
  Beaker,
  Brain,
  Database,
  FileText,
  Microscope,
  Pill,
  Search,
  Shield,
  Mic,
  Users,
  MessageCircle,
  BarChart3,
  Star,
  Mail,
  Phone,
  MapPin,
  Zap,
  Target,
  TrendingUp,
} from "lucide-react"
import { motion, useScroll, useTransform, useInView, useMotionValue, useSpring } from "framer-motion"
import Navbar from "../../components/Navbar"
import { Link } from 'react-router-dom';

function Homepage() {
  const [activeFeature, setActiveFeature] = useState(null)
  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 })
  const [counters, setCounters] = useState({ stat1: 0, stat2: 0, stat3: 0, stat4: 0 })

  // Refs for scroll animations
  const statsRef = useRef(null)
  const featuresRef = useRef(null)
  const testimonialsRef = useRef(null)
  const ctaRef = useRef(null)

  // InView hooks for sections
  const statsInView = useInView(statsRef, { once: true, amount: 0.3 })
  const featuresInView = useInView(featuresRef, { once: true, amount: 0.1 })
  const testimonialsInView = useInView(testimonialsRef, { once: true, amount: 0.2 })
  const ctaInView = useInView(ctaRef, { once: true, amount: 0.5 })

  // Scroll progress for parallax effects
  const { scrollYProgress } = useScroll()
  const heroOpacity = useTransform(scrollYProgress, [0, 0.2], [1, 0])
  const heroY = useTransform(scrollYProgress, [0, 0.2], [0, 100])

  // Mouse tracking for hero parallax effect
  useEffect(() => {
    const handleMouseMove = (e) => {
      setMousePosition({
        x: (e.clientX / window.innerWidth) * 100,
        y: (e.clientY / window.innerHeight) * 100,
      })
    }

    window.addEventListener("mousemove", handleMouseMove)
    return () => window.removeEventListener("mousemove", handleMouseMove)
  }, [])

  // Counter animation for stats using Framer Motion
  const value1 = useMotionValue(0)
  const value2 = useMotionValue(0)
  const value3 = useMotionValue(0)
  const value4 = useMotionValue(0)

  const spring1 = useSpring(value1, { duration: 2000 })
  const spring2 = useSpring(value2, { duration: 2000 })
  const spring3 = useSpring(value3, { duration: 2000 })
  const spring4 = useSpring(value4, { duration: 2000 })

  useEffect(() => {
    spring1.set(85)
  }, [spring1])

  useEffect(() => {
    return spring1.onChange((latest) => {
      setCounters((prev) => ({ ...prev, stat1: Math.floor(latest) }))
    })
  }, [spring1])

  useEffect(() => {
    spring2.set(60)
  }, [spring2])

  useEffect(() => {
    return spring2.onChange((latest) => {
      setCounters((prev) => ({ ...prev, stat2: Math.floor(latest) }))
    })
  }, [spring2])

  useEffect(() => {
    spring3.set(3.5)
  }, [spring3])

  useEffect(() => {
    return spring3.onChange((latest) => {
      setCounters((prev) => ({ ...prev, stat3: Math.floor(latest) }))
    })
  }, [spring3])

  useEffect(() => {
    if (statsInView) {
      setTimeout(() => setCounters((prev) => ({ ...prev, stat4: "24/7" })), 1000)
    }
  }, [statsInView])

  // Animation variants
  const containerVariants = {
    hidden: { opacity: 0 },
    visible: {
      opacity: 1,
      transition: {
        staggerChildren: 0.1,
        delayChildren: 0.2,
      },
    },
  }

  const itemVariants = {
    hidden: { opacity: 0, y: 20 },
    visible: {
      opacity: 1,
      y: 0,
      transition: {
        type: "spring",
        stiffness: 100,
        damping: 10,
      },
    },
  }

  const fadeInVariants = {
    hidden: { opacity: 0 },
    visible: {
      opacity: 1,
      transition: {
        duration: 0.6,
      },
    },
  }

  const scaleVariants = {
    hidden: { scale: 0.8, opacity: 0 },
    visible: {
      scale: 1,
      opacity: 1,
      transition: {
        type: "spring",
        stiffness: 100,
        damping: 10,
      },
    },
  }

  const features = [
    {
      icon: <Brain className="w-6 h-6" />,
      title: "AI-Powered Molecule Generation",
      description: "Generate novel molecules using SMILES notation with advanced generative AI algorithms.",
      color: "bg-blue-50 text-blue-600",
    },
    {
      icon: <Microscope className="w-6 h-6" />,
      title: "Molecule Visualization & Analysis",
      description: "Visualize and analyze molecules in real-time using RDKit for data-driven insights.",
      color: "bg-green-50 text-green-600",
    },
    {
      icon: <Beaker className="w-6 h-6" />,
      title: "Molecular Evolution",
      description: "Iteratively mutate molecules to optimize for stability, solubility, and efficacy.",
      color: "bg-purple-50 text-purple-600",
    },
    {
      icon: <BarChart3 className="w-6 h-6" />,
      title: "Synthesis Cost Estimation",
      description: "Estimate real-world synthesis costs by analyzing lab materials and complexity.",
      color: "bg-orange-50 text-orange-600",
    },
    {
      icon: <FileText className="w-6 h-6" />,
      title: "Research Paper Generator",
      description: "Automatically generate research papers in IEEE, APA, or Nature journal styles.",
      color: "bg-red-50 text-red-600",
    },
    {
      icon: <Search className="w-6 h-6" />,
      title: "Research Summaries",
      description: "Summarize the latest research from PubMed, arXiv, and Google Scholar using Gemini AI.",
      color: "bg-indigo-50 text-indigo-600",
    },
    {
      icon: <Database className="w-6 h-6" />,
      title: "Discovery Recommendations",
      description: "Get AI-driven recommendations for new molecules based on your previous work.",
      color: "bg-teal-50 text-teal-600",
    },
    {
      icon: <Pill className="w-6 h-6" />,
      title: "Drug Naming",
      description: "Get intelligent, systematic names for your novel drug candidates.",
      color: "bg-pink-50 text-pink-600",
    },
    {
      icon: <Shield className="w-6 h-6" />,
      title: "Toxicity Prediction",
      description: "Predict potential toxicity and side effects in real-time.",
      color: "bg-emerald-50 text-emerald-600",
    },
    {
      icon: <Mic className="w-6 h-6" />,
      title: "Voice-to-Text Notes",
      description: "Dictate observations, and AI will auto-transcribe and summarize them for you.",
      color: "bg-yellow-50 text-yellow-600",
    },
    {
      icon: <Users className="w-6 h-6" />,
      title: "Real-Time Collaboration",
      description: "Collaborate seamlessly with your team through shared workspaces.",
      color: "bg-cyan-50 text-cyan-600",
    },
    {
      icon: <MessageCircle className="w-6 h-6" />,
      title: "AI Chatbot",
      description: "A Gemini-powered assistant to explain complex molecular structures.",
      color: "bg-violet-50 text-violet-600",
    },
  ]

  return (
    <>
      <style jsx>{`
        @keyframes rotate-slow {
          from {
            transform: rotate(0deg);
          }
          to {
            transform: rotate(360deg);
          }
        }

        @keyframes pulse-glow {
          0% {
            box-shadow: 0 0 0 0 rgba(79, 70, 229, 0.4);
          }
          50% {
            box-shadow: 0 0 20px 10px rgba(79, 70, 229, 0.6);
          }
          100% {
            box-shadow: 0 0 0 0 rgba(79, 70, 229, 0.4);
          }
        }

        @keyframes bounce-gentle {
          0%,
          100% {
            transform: translateY(0);
          }
          50% {
            transform: translateY(-5px);
          }
        }

        @keyframes float {
          0%,
          100% {
            transform: translateY(0);
          }
          50% {
            transform: translateY(-10px);
          }
        }

        @keyframes slide-up {
          from {
            opacity: 0;
            transform: translateY(30px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }

        @keyframes slide-in-right {
          from {
            opacity: 0;
            transform: translateX(30px);
          }
          to {
            opacity: 1;
            transform: translateX(0);
          }
        }

        @keyframes fade-in {
          from {
            opacity: 0;
          }
          to {
            opacity: 1;
          }
        }

        .gradient-text {
          background: linear-gradient(to right, #2563eb, #9333ea);
          -webkit-background-clip: text;
          background-clip: text;
          -webkit-text-fill-color: transparent;
          color: transparent;
        }

        .animate-rotate-slow {
          animation: rotate-slow 20s linear infinite;
        }

        .animate-pulse-glow {
          animation: pulse-glow 3s infinite;
        }

        .animate-bounce-gentle {
          animation: bounce-gentle 3s infinite ease-in-out;
        }

        .animate-float {
          animation: float 5s infinite ease-in-out;
        }

        .animate-slide-up {
          animation: slide-up 0.6s ease-out forwards;
        }

        .animate-slide-in-right {
          animation: slide-in-right 0.8s ease-out forwards;
        }

        .animate-fade-in {
          animation: fade-in 0.6s ease-out forwards;
        }

        .delay-100 {
          animation-delay: 0.1s;
        }

        .delay-200 {
          animation-delay: 0.2s;
        }

        .delay-300 {
          animation-delay: 0.3s;
        }

        .delay-400 {
          animation-delay: 0.4s;
        }

        .delay-500 {
          animation-delay: 0.5s;
        }

        .delay-600 {
          animation-delay: 0.6s;
        }

        .hover-lift {
          transition: transform 0.3s ease, box-shadow 0.3s ease;
        }

        .hover-lift:hover {
          transform: translateY(-10px);
          box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
        }

        .feature-card {
          transition: all 0.3s ease;
        }

        .feature-card:hover {
          transform: translateY(-10px);
          box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
        }

        .floating-shapes {
          position: absolute;
          inset: 0;
          overflow: hidden;
          pointer-events: none;
        }

        .shape {
          position: absolute;
          opacity: 0.3;
        }

        .shape-1 {
          top: 10%;
          left: 5%;
        }

        .shape-2 {
          top: 60%;
          right: 10%;
        }

        .shape-3 {
          bottom: 15%;
          left: 15%;
        }

        .shape-4 {
          top: 25%;
          right: 20%;
        }

        .parallax-element {
          transition: transform 0.1s ease-out;
        }
      `}</style>

      <div className="min-h-screen bg-white overflow-hidden">
        <Navbar />

        {/* Hero Section */}
        <motion.section
          className="relative bg-gradient-to-br mt-3 from-blue-50 via-white to-purple-50 py-20 overflow-hidden"
          style={{ opacity: heroOpacity, y: heroY }}
        >
          {/* Floating Background Shapes */}
          <div className="floating-shapes">
            <motion.div
              className="shape shape-1"
              animate={{
                x: [0, 10, 0],
                y: [0, -15, 0],
                rotate: [0, 5, 0],
              }}
              transition={{
                repeat: Number.POSITIVE_INFINITY,
                duration: 8,
                ease: "easeInOut",
              }}
            >
              <Atom className="w-24 h-24 text-blue-300" />
            </motion.div>
            <motion.div
              className="shape shape-2"
              animate={{
                x: [0, -15, 0],
                y: [0, 10, 0],
                rotate: [0, -5, 0],
              }}
              transition={{
                repeat: Number.POSITIVE_INFINITY,
                duration: 10,
                ease: "easeInOut",
              }}
            >
              <Brain className="w-32 h-32 text-purple-300" />
            </motion.div>
            <motion.div
              className="shape shape-3"
              animate={{
                x: [0, 20, 0],
                y: [0, 15, 0],
                rotate: [0, 10, 0],
              }}
              transition={{
                repeat: Number.POSITIVE_INFINITY,
                duration: 12,
                ease: "easeInOut",
              }}
            >
              <Beaker className="w-20 h-20 text-green-300" />
            </motion.div>
            <motion.div
              className="shape shape-4"
              animate={{
                x: [0, -10, 0],
                y: [0, -20, 0],
                rotate: [0, -8, 0],
              }}
              transition={{
                repeat: Number.POSITIVE_INFINITY,
                duration: 9,
                ease: "easeInOut",
              }}
            >
              <Microscope className="w-28 h-28 text-orange-300" />
            </motion.div>
          </div>

          {/* Animated Particles */}
          <div className="absolute inset-0 overflow-hidden pointer-events-none">
            {[...Array(20)].map((_, i) => (
              <motion.div
                key={i}
                className="absolute rounded-full bg-blue-400 opacity-20"
                initial={{
                  x: Math.random() * 100 + "%",
                  y: Math.random() * 100 + "%",
                  scale: Math.random() * 0.5 + 0.5,
                }}
                animate={{
                  y: [Math.random() * 100 + "%", Math.random() * 100 + "%"],
                  x: [Math.random() * 100 + "%", Math.random() * 100 + "%"],
                }}
                transition={{
                  duration: Math.random() * 20 + 10,
                  repeat: Number.POSITIVE_INFINITY,
                  repeatType: "reverse",
                  ease: "linear",
                }}
                style={{
                  width: Math.random() * 4 + 2 + "px",
                  height: Math.random() * 4 + 2 + "px",
                }}
              />
            ))}
          </div>

          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 relative">
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
              <motion.div className="space-y-8" initial="hidden" animate="visible" variants={containerVariants}>
                <motion.div variants={itemVariants}>
                  <div className="inline-flex items-center px-4 py-2 rounded-full bg-blue-100 text-blue-800 text-sm font-medium mb-6 hover:bg-blue-200 transition-colors duration-300">
                    <motion.span
                      className="w-2 h-2 bg-blue-600 rounded-full mr-2"
                      animate={{ scale: [1, 1.5, 1] }}
                      transition={{ duration: 2, repeat: Number.POSITIVE_INFINITY }}
                    ></motion.span>
                    Next-Gen AI for Drug Discovery
                  </div>
                </motion.div>

                <motion.h1
                  className="text-5xl lg:text-6xl font-bold text-gray-900 leading-tight"
                  variants={containerVariants}
                >
                  <motion.span variants={itemVariants} className="inline-block">
                    Revolutionizing
                  </motion.span>
                  <motion.span variants={itemVariants} className="gradient-text block">
                    Drug Discovery
                  </motion.span>
                  <motion.span variants={itemVariants} className="inline-block">
                    with AI
                  </motion.span>
                </motion.h1>

                <motion.p className="text-xl text-gray-600 leading-relaxed" variants={itemVariants}>
                  Harness the power of generative AI to accelerate your research, reduce costs, and discover
                  breakthrough medications faster than ever before.
                </motion.p>

                <motion.div className="flex flex-col sm:flex-row gap-4" variants={containerVariants}>
                  <motion.button
                    className="group bg-blue-600 text-white px-8 py-4 rounded-lg hover:bg-blue-700 transition-all duration-300 flex items-center justify-center font-semibold"
                    variants={itemVariants}
                    whileHover={{ scale: 1.05, boxShadow: "0 10px 25px -5px rgba(59, 130, 246, 0.5)" }}
                    whileTap={{ scale: 0.98 }}
                  >
                    <Link to="/dashboard" className="flex items-center">
                      Start Your Discovery Journey
                      <motion.span
                        className="ml-2 inline-block"
                        whileHover={{ x: 5 }}
                        transition={{ type: "spring", stiffness: 400, damping: 10 }}
                      >
                        <ArrowRight className="w-5 h-5" />
                      </motion.span>
                    </Link>
                  </motion.button>
                  <motion.button
                    className="group border-2 border-gray-300 text-gray-700 px-8 py-4 rounded-lg hover:border-blue-600 hover:text-blue-600 transition-all duration-300 font-semibold"
                    variants={itemVariants}
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.98 }}
                  >
                    Explore Features
                    <motion.span
                      className="ml-2 inline-block"
                      whileHover={{ rotate: 12 }}
                      transition={{ type: "spring", stiffness: 400, damping: 10 }}
                    >
                      <Zap className="w-5 h-5" />
                    </motion.span>
                  </motion.button>
                </motion.div>

                <motion.div className="flex items-center space-x-6" variants={itemVariants}>
                  <div className="flex -space-x-2">
                    {[1, 2, 3, 4].map((i) => (
                      <motion.div
                        key={i}
                        className="w-10 h-10 rounded-full bg-gray-300 border-2 border-white flex items-center justify-center text-gray-600 font-semibold cursor-pointer"
                        initial={{ opacity: 0, scale: 0.5 }}
                        animate={{ opacity: 1, scale: 1 }}
                        transition={{ delay: i * 0.1 }}
                        whileHover={{ scale: 1.1, y: -5 }}
                      >
                        {i}
                      </motion.div>
                    ))}
                  </div>
                  <div className="text-sm text-gray-600">
                    <motion.span
                      className="font-semibold text-gray-900"
                      initial={{ opacity: 0 }}
                      animate={{ opacity: 1 }}
                      transition={{ delay: 0.6 }}
                    >
                      500+
                    </motion.span>{" "}
                    researchers already using our platform
                  </div>
                </motion.div>
              </motion.div>

              <motion.div
                className="relative"
                initial={{ opacity: 0, x: 50 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ duration: 0.8, delay: 0.3 }}
              >
                {/* Main Dashboard Mockup */}
                <motion.div
                  className="bg-white rounded-2xl shadow-2xl p-8 border border-gray-100"
                  style={{
                    x: useTransform(useMotionValue(mousePosition.x), [0, 100], [-15, 15]),
                    y: useTransform(useMotionValue(mousePosition.y), [0, 100], [-15, 15]),
                  }}
                  whileHover={{ scale: 1.02 }}
                  transition={{ type: "spring", stiffness: 300, damping: 20 }}
                >
                  {/* Header */}
                  <div className="flex items-center justify-between mb-6">
                    <div className="flex items-center space-x-2">
                      <motion.div className="w-3 h-3 bg-red-400 rounded-full" whileHover={{ scale: 1.5 }}></motion.div>
                      <motion.div
                        className="w-3 h-3 bg-yellow-400 rounded-full"
                        whileHover={{ scale: 1.5 }}
                      ></motion.div>
                      <motion.div
                        className="w-3 h-3 bg-green-400 rounded-full"
                        whileHover={{ scale: 1.5 }}
                      ></motion.div>
                    </div>
                    <div className="text-sm text-gray-500">Drug Discovery Dashboard</div>
                  </div>

                  {/* Molecule Grid */}
                  <div className="grid grid-cols-3 gap-4 mb-6">
                    {[1, 2, 3, 4, 5, 6].map((i) => (
                      <motion.div
                        key={i}
                        className="h-16 bg-gradient-to-br from-blue-50 to-purple-50 rounded-lg flex items-center justify-center"
                        whileHover={{ scale: 1.1, rotate: 5 }}
                        initial={{ opacity: 0, scale: 0.8 }}
                        animate={{ opacity: 1, scale: 1 }}
                        transition={{ delay: i * 0.1 }}
                      >
                        <motion.div
                          animate={{ rotate: 360 }}
                          transition={{ duration: 20, repeat: Number.POSITIVE_INFINITY, ease: "linear" }}
                        >
                          <Beaker className="w-6 h-6 text-blue-500" />
                        </motion.div>
                      </motion.div>
                    ))}
                  </div>

                  {/* Progress Bars */}
                  <div className="space-y-3">
                    <div className="flex items-center justify-between text-sm text-gray-600">
                      <span>AI Analysis Progress</span>
                      <span>87%</span>
                    </div>
                    <motion.div
                      className="h-2 bg-gray-100 rounded-full overflow-hidden"
                      initial={{ width: 0 }}
                      animate={{ width: "100%" }}
                      transition={{ duration: 1 }}
                    >
                      <motion.div
                        className="h-full bg-gradient-to-r from-blue-500 to-purple-500 rounded-full"
                        initial={{ width: 0 }}
                        animate={{ width: "87%" }}
                        transition={{ duration: 1.5, delay: 0.2 }}
                      ></motion.div>
                    </motion.div>

                    <div className="flex items-center justify-between text-sm text-gray-600">
                      <span>Synthesis Optimization</span>
                      <span>64%</span>
                    </div>
                    <motion.div
                      className="h-2 bg-gray-100 rounded-full overflow-hidden"
                      initial={{ width: 0 }}
                      animate={{ width: "100%" }}
                      transition={{ duration: 1, delay: 0.3 }}
                    >
                      <motion.div
                        className="h-full bg-gradient-to-r from-green-500 to-teal-500 rounded-full"
                        initial={{ width: 0 }}
                        animate={{ width: "64%" }}
                        transition={{ duration: 1.5, delay: 0.5 }}
                      ></motion.div>
                    </motion.div>

                    <div className="flex items-center justify-between text-sm text-gray-600">
                      <span>Toxicity Prediction</span>
                      <span>92%</span>
                    </div>
                    <motion.div
                      className="h-2 bg-gray-100 rounded-full overflow-hidden"
                      initial={{ width: 0 }}
                      animate={{ width: "100%" }}
                      transition={{ duration: 1, delay: 0.6 }}
                    >
                      <motion.div
                        className="h-full bg-gradient-to-r from-purple-500 to-pink-500 rounded-full"
                        initial={{ width: 0 }}
                        animate={{ width: "92%" }}
                        transition={{ duration: 1.5, delay: 0.8 }}
                      ></motion.div>
                    </motion.div>
                  </div>
                </motion.div>

                {/* Floating AI Badge */}
                <motion.div
                  className="absolute -top-4 -right-4 w-24 h-24 bg-gradient-to-br from-blue-600 to-purple-600 rounded-full flex items-center justify-center"
                  animate={{
                    boxShadow: [
                      "0 0 0 rgba(79, 70, 229, 0.4)",
                      "0 0 20px rgba(79, 70, 229, 0.6)",
                      "0 0 0 rgba(79, 70, 229, 0.4)",
                    ],
                  }}
                  transition={{
                    duration: 2,
                    repeat: Number.POSITIVE_INFINITY,
                    repeatType: "reverse",
                  }}
                  style={{
                    x: useTransform(useMotionValue(mousePosition.x), [0, 100], [10, -10]),
                    y: useTransform(useMotionValue(mousePosition.y), [0, 100], [10, -10]),
                  }}
                >
                  <motion.div
                    animate={{ rotate: 360 }}
                    transition={{ duration: 15, repeat: Number.POSITIVE_INFINITY, ease: "linear" }}
                  >
                    <Brain className="w-12 h-12 text-white" />
                  </motion.div>
                </motion.div>

                {/* Floating Stats Cards */}
                <motion.div
                  className="absolute -bottom-4 -left-4 bg-white rounded-xl p-4 shadow-lg border border-gray-100"
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: 0.8 }}
                  style={{
                    x: useTransform(useMotionValue(mousePosition.x), [0, 100], [-5, 5]),
                    y: useTransform(useMotionValue(mousePosition.y), [0, 100], [-5, 5]),
                  }}
                >
                  <div className="flex items-center space-x-2">
                    <motion.div
                      animate={{ rotate: [0, 10, 0] }}
                      transition={{ duration: 2, repeat: Number.POSITIVE_INFINITY }}
                    >
                      <TrendingUp className="w-5 h-5 text-green-500" />
                    </motion.div>
                    <div>
                      <motion.div
                        className="text-lg font-bold text-gray-900"
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        transition={{ delay: 1 }}
                      >
                        3.5x
                      </motion.div>
                      <div className="text-xs text-gray-500">Faster Discovery</div>
                    </div>
                  </div>
                </motion.div>

                <motion.div
                  className="absolute top-1/2 -right-8 bg-white rounded-xl p-4 shadow-lg border border-gray-100"
                  initial={{ opacity: 0, y: -20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: 1 }}
                  style={{
                    x: useTransform(useMotionValue(mousePosition.x), [0, 100], [5, -5]),
                    y: useTransform(useMotionValue(mousePosition.y), [0, 100], [5, -5]),
                  }}
                >
                  <div className="flex items-center space-x-2">
                    <motion.div
                      animate={{ scale: [1, 1.2, 1] }}
                      transition={{ duration: 2, repeat: Number.POSITIVE_INFINITY }}
                    >
                      <Target className="w-5 h-5 text-blue-500" />
                    </motion.div>
                    <div>
                      <motion.div
                        className="text-lg font-bold text-gray-900"
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        transition={{ delay: 1.2 }}
                      >
                        85%
                      </motion.div>
                      <div className="text-xs text-gray-500">Accuracy</div>
                    </div>
                  </div>
                </motion.div>
              </motion.div>
            </div>
          </div>
        </motion.section>

        {/* Stats Section */}
        <section className="py-16 bg-white border-t border-gray-100" id="stats" ref={statsRef}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <motion.div
              className="grid grid-cols-2 md:grid-cols-4 gap-8"
              initial="hidden"
              animate={statsInView ? "visible" : "hidden"}
              variants={containerVariants}
            >
              {[
                { key: "stat1", value: counters.stat1, suffix: "%", label: "Faster Discovery", color: "text-blue-600" },
                { key: "stat2", value: counters.stat2, suffix: "%", label: "Cost Reduction", color: "text-green-600" },
                {
                  key: "stat3",
                  value: counters.stat3,
                  suffix: "x",
                  label: "More Candidates",
                  color: "text-purple-600",
                },
                { key: "stat4", value: counters.stat4, suffix: "", label: "AI Assistance", color: "text-orange-600" },
              ].map((stat, index) => (
                <motion.div key={index} className="text-center" variants={itemVariants} custom={index}>
                  <motion.div
                    className={`text-4xl font-bold ${stat.color} mb-2`}
                    initial={{ scale: 0.5, opacity: 0 }}
                    animate={statsInView ? { scale: 1, opacity: 1 } : { scale: 0.5, opacity: 0 }}
                    transition={{ delay: index * 0.1 + 0.2 }}
                  >
                    {stat.value}
                    {stat.suffix}
                  </motion.div>
                  <motion.div
                    className="text-gray-600"
                    initial={{ opacity: 0 }}
                    animate={statsInView ? { opacity: 1 } : { opacity: 0 }}
                    transition={{ delay: index * 0.1 + 0.4 }}
                  >
                    {stat.label}
                  </motion.div>
                </motion.div>
              ))}
            </motion.div>
          </div>
        </section>

        {/* Features Section */}
        <section id="features" className="py-20 bg-gray-50" ref={featuresRef}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <motion.div
              className="text-center mb-16"
              initial="hidden"
              animate={featuresInView ? "visible" : "hidden"}
              variants={fadeInVariants}
            >
              <motion.div
                className="inline-flex items-center px-4 py-2 rounded-full bg-blue-100 text-blue-800 text-sm font-medium mb-6"
                variants={scaleVariants}
              >
                Powerful Capabilities
              </motion.div>
              <motion.h2 className="text-4xl font-bold text-gray-900 mb-4" variants={itemVariants}>
                Accelerate Your Research
              </motion.h2>
              <motion.p className="text-xl text-gray-600 max-w-3xl mx-auto" variants={itemVariants}>
                Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
              </motion.p>
            </motion.div>

            <motion.div
              className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8"
              initial="hidden"
              animate={featuresInView ? "visible" : "hidden"}
              variants={containerVariants}
            >
              {features.map((feature, index) => (
                <motion.div
                  key={index}
                  className="feature-card bg-white rounded-xl p-8 shadow-sm border border-gray-100"
                  variants={itemVariants}
                  custom={index}
                  whileHover={{
                    y: -10,
                    boxShadow: "0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
                  }}
                  onMouseEnter={() => setActiveFeature(index)}
                  onMouseLeave={() => setActiveFeature(null)}
                >
                  <motion.div
                    className={`inline-flex items-center justify-center w-12 h-12 rounded-lg ${feature.color} mb-6`}
                    whileHover={{ scale: 1.1, rotate: 6 }}
                    animate={activeFeature === index ? { scale: 1.1, rotate: 6 } : { scale: 1, rotate: 0 }}
                  >
                    {feature.icon}
                  </motion.div>

                  <motion.h3
                    className="text-xl font-semibold text-gray-900 mb-3"
                    animate={activeFeature === index ? { color: "#3B82F6" } : { color: "#111827" }}
                  >
                    {feature.title}
                  </motion.h3>

                  <motion.p className="text-gray-600 leading-relaxed mb-4">{feature.description}</motion.p>

                  <motion.div
                    className="flex items-center text-blue-600 font-medium"
                    animate={activeFeature === index ? { x: 8 } : { x: 0 }}
                  >
                    <span className="text-sm">Learn more</span>
                    <ArrowRight className="w-4 h-4 ml-2" />
                  </motion.div>
                </motion.div>
              ))}
            </motion.div>
          </div>
        </section>

        {/* Testimonials Section */}
        <section id="testimonials" className="py-20 bg-white" ref={testimonialsRef}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <motion.div
              className="text-center mb-16"
              initial="hidden"
              animate={testimonialsInView ? "visible" : "hidden"}
              variants={fadeInVariants}
            >
              <motion.div
                className="inline-flex items-center px-4 py-2 rounded-full bg-green-100 text-green-800 text-sm font-medium mb-6"
                variants={scaleVariants}
              >
                Success Stories
              </motion.div>
              <motion.h2 className="text-4xl font-bold text-gray-900 mb-4" variants={itemVariants}>
                What Researchers Are Saying
              </motion.h2>
            </motion.div>

            <motion.div
              className="grid grid-cols-1 md:grid-cols-3 gap-8"
              initial="hidden"
              animate={testimonialsInView ? "visible" : "hidden"}
              variants={containerVariants}
            >
              {[
                {
                  quote:
                    "This platform reduced our initial discovery phase from 6 months to just 2 weeks. Unbelievable!",
                  author: "Dr. Sarah Chen",
                  role: "Lead Researcher",
                  affiliation: "Stanford University",
                  rating: 5,
                },
                {
                  quote: "The AI-generated molecules had better binding affinity than our manually designed compounds.",
                  author: "Prof. Raj Patel",
                  role: "Principal Investigator",
                  affiliation: "MIT Bioengineering",
                  rating: 5,
                },
                {
                  quote:
                    "Finally, a tool that bridges the gap between computational chemistry and practical drug development.",
                  author: "Dr. Elena Rodriguez",
                  role: "Senior Scientist",
                  affiliation: "Novartis Pharmaceuticals",
                  rating: 5,
                },
              ].map((testimonial, index) => (
                <motion.div
                  key={index}
                  className="bg-gray-50 rounded-xl p-8 border border-gray-100"
                  variants={itemVariants}
                  custom={index}
                  whileHover={{
                    y: -10,
                    boxShadow: "0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
                  }}
                >
                  <motion.div className="flex mb-4">
                    {[...Array(testimonial.rating)].map((_, i) => (
                      <motion.div
                        key={i}
                        animate={{ scale: [1, 1.2, 1] }}
                        transition={{
                          duration: 2,
                          repeat: Number.POSITIVE_INFINITY,
                          repeatType: "reverse",
                          delay: i * 0.1,
                        }}
                      >
                        <Star className="w-5 h-5 text-yellow-400 fill-current" />
                      </motion.div>
                    ))}
                  </motion.div>

                  <motion.p
                    className="text-gray-700 text-lg mb-6 leading-relaxed"
                    initial={{ opacity: 0 }}
                    animate={testimonialsInView ? { opacity: 1 } : { opacity: 0 }}
                    transition={{ delay: index * 0.2 + 0.3 }}
                  >
                    "{testimonial.quote}"
                  </motion.p>

                  <motion.div
                    className="flex items-center"
                    initial={{ opacity: 0, x: -20 }}
                    animate={testimonialsInView ? { opacity: 1, x: 0 } : { opacity: 0, x: -20 }}
                    transition={{ delay: index * 0.2 + 0.5 }}
                  >
                    <motion.div
                      className="w-12 h-12 rounded-full bg-blue-600 flex items-center justify-center text-white font-semibold mr-4"
                      whileHover={{ scale: 1.1, rotate: 10 }}
                    >
                      {testimonial.author.charAt(0)}
                    </motion.div>
                    <div>
                      <h4 className="font-semibold text-gray-900">{testimonial.author}</h4>
                      <p className="text-gray-600 text-sm">{testimonial.role}</p>
                      <p className="text-gray-500 text-sm">{testimonial.affiliation}</p>
                    </div>
                  </motion.div>
                </motion.div>
              ))}
            </motion.div>
          </div>
        </section>

        {/* CTA Section */}
        <motion.section
          className="py-20 bg-gradient-to-r from-blue-600 to-purple-600 relative overflow-hidden"
          ref={ctaRef}
        >
          <motion.div
            className="absolute inset-0"
            initial={{ opacity: 0 }}
            animate={ctaInView ? { opacity: 0.1 } : { opacity: 0 }}
          >
            {[...Array(30)].map((_, i) => (
              <motion.div
                key={i}
                className="absolute rounded-full bg-white"
                initial={{
                  x: Math.random() * 100 + "%",
                  y: Math.random() * 100 + "%",
                  opacity: Math.random() * 0.5 + 0.1,
                  scale: Math.random() * 0.5 + 0.5,
                }}
                animate={{
                  y: [Math.random() * 100 + "%", Math.random() * 100 + "%"],
                }}
                transition={{
                  duration: Math.random() * 10 + 10,
                  repeat: Number.POSITIVE_INFINITY,
                  repeatType: "reverse",
                  ease: "linear",
                }}
                style={{
                  width: Math.random() * 8 + 2 + "px",
                  height: Math.random() * 8 + 2 + "px",
                }}
              />
            ))}
          </motion.div>

          <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center relative">
            <motion.div
              className="inline-flex items-center px-4 py-2 rounded-full bg-blue-500 text-blue-100 text-sm font-medium mb-6"
              initial={{ opacity: 0, y: 20 }}
              animate={ctaInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 20 }}
              transition={{ delay: 0.2 }}
            >
              <motion.span
                animate={{ scale: [1, 1.2, 1] }}
                transition={{ duration: 2, repeat: Number.POSITIVE_INFINITY }}
              >
                Limited Time Offer
              </motion.span>
            </motion.div>

            <motion.h2
              className="text-4xl font-bold text-white mb-6"
              initial={{ opacity: 0, y: 20 }}
              animate={ctaInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 20 }}
              transition={{ delay: 0.4 }}
            >
              Ready to Transform Drug Discovery?
            </motion.h2>

            <motion.p
              className="text-xl text-blue-100 mb-10 max-w-3xl mx-auto"
              initial={{ opacity: 0, y: 20 }}
              animate={ctaInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 20 }}
              transition={{ delay: 0.6 }}
            >
              Join thousands of researchers using our platform to accelerate their drug discovery process with the power
              of generative AI.
            </motion.p>

            <motion.div
              className="flex flex-col sm:flex-row justify-center gap-4"
              initial={{ opacity: 0, y: 20 }}
              animate={ctaInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 20 }}
              transition={{ delay: 0.8 }}
            >
              <motion.button
                className="group bg-white text-blue-600 px-8 py-4 rounded-lg hover:bg-gray-50 transition-all duration-300 font-semibold flex items-center justify-center"
                whileHover={{ scale: 1.05, boxShadow: "0 10px 25px -5px rgba(255, 255, 255, 0.5)" }}
                whileTap={{ scale: 0.98 }}
              >
                Get Started Now
                <motion.span
                  className="ml-2 inline-block"
                  whileHover={{ x: 5 }}
                  transition={{ type: "spring", stiffness: 400, damping: 10 }}
                >
                  <ArrowRight className="w-5 h-5" />
                </motion.span>
              </motion.button>
              <motion.button
                className="group border-2 border-blue-400 text-white px-8 py-4 rounded-lg hover:bg-blue-500 transition-all duration-300 font-semibold"
                whileHover={{ scale: 1.05 }}
                whileTap={{ scale: 0.98 }}
              >
                Contact Our Team
                <motion.span
                  className="ml-2 inline-block"
                  whileHover={{ rotate: 12 }}
                  transition={{ type: "spring", stiffness: 400, damping: 10 }}
                >
                  <MessageCircle className="w-5 h-5" />
                </motion.span>
              </motion.button>
            </motion.div>
          </div>
        </motion.section>

        {/* Footer */}
        <footer className="bg-gray-900 text-white py-16">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="grid grid-cols-1 md:grid-cols-4 gap-12">
              <div className="space-y-4">
                <motion.div className="flex items-center" whileHover={{ scale: 1.05 }}>
                  <motion.div
                    animate={{ rotate: 360 }}
                    transition={{ duration: 20, repeat: Number.POSITIVE_INFINITY, ease: "linear" }}
                  >
                    <Atom className="w-8 h-8 text-blue-400 mr-3" />
                  </motion.div>
                  <span className="text-xl font-bold">Drug Discovery AI</span>
                </motion.div>
                <p className="text-gray-400 leading-relaxed">
                  Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more
                  efficient.
                </p>
              </div>

              <div>
                <h3 className="text-lg font-semibold mb-4">Features</h3>
                <ul className="space-y-3">
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#features" className="text-gray-400 hover:text-white transition-all duration-300">
                      Molecule Generation
                    </a>
                  </motion.li>
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#features" className="text-gray-400 hover:text-white transition-all duration-300">
                      Visualization Tools
                    </a>
                  </motion.li>
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#features" className="text-gray-400 hover:text-white transition-all duration-300">
                      Research Automation
                    </a>
                  </motion.li>
                </ul>
              </div>

              <div>
                <h3 className="text-lg font-semibold mb-4">Company</h3>
                <ul className="space-y-3">
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#" className="text-gray-400 hover:text-white transition-all duration-300">
                      About Us
                    </a>
                  </motion.li>
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#" className="text-gray-400 hover:text-white transition-all duration-300">
                      Careers
                    </a>
                  </motion.li>
                  <motion.li whileHover={{ x: 5 }} transition={{ type: "spring", stiffness: 400, damping: 10 }}>
                    <a href="#" className="text-gray-400 hover:text-white transition-all duration-300">
                      Blog
                    </a>
                  </motion.li>
                </ul>
              </div>

              <div id="contact">
                <h3 className="text-lg font-semibold mb-4">Contact</h3>
                <div className="space-y-3 text-gray-400">
                  <motion.div
                    className="flex items-center hover:text-white transition-colors duration-300"
                    whileHover={{ x: 5 }}
                  >
                    <Mail className="w-4 h-4 mr-2" />
                    <span>alokchaturvedi190@gmail.com</span>
                  </motion.div>
                  <motion.div
                    className="flex items-center hover:text-white transition-colors duration-300"
                    whileHover={{ x: 5 }}
                  >
                    <Phone className="w-4 h-4 mr-2" />
                    <span>+91 9975175098</span>
                  </motion.div>
                  <motion.div
                    className="flex items-center hover:text-white transition-colors duration-300"
                    whileHover={{ x: 5 }}
                  >
                    <Mail className="w-4 h-4 mr-2" />
                    <span>tbhangale9@gmail.com</span>
                  </motion.div>
                  <motion.div
                    className="flex items-center hover:text-white transition-colors duration-300"
                    whileHover={{ x: 5 }}
                  >
                    <Phone className="w-4 h-4 mr-2" />
                    <span>+91 8766816061</span>
                  </motion.div>
                  <motion.div
                    className="flex items-center hover:text-white transition-colors duration-300"
                    whileHover={{ x: 5 }}
                  >
                    <MapPin className="w-4 h-4 mr-2" />
                    <span>Pune, Maharashtra, India</span>
                  </motion.div>
                </div>
              </div>
            </div>

            <div className="mt-12 pt-8 border-t border-gray-800 flex flex-col md:flex-row justify-between items-center">
              <p className="text-gray-400">© {new Date().getFullYear()} Drug Discovery AI. All rights reserved.</p>
              <div className="mt-4 md:mt-0">
                <ul className="flex space-x-6 text-sm">
                  <motion.li whileHover={{ y: -3 }}>
                    <a href="#" className="text-gray-400 hover:text-white transition-colors duration-300">
                      Terms
                    </a>
                  </motion.li>
                  <motion.li whileHover={{ y: -3 }}>
                    <a href="#" className="text-gray-400 hover:text-white transition-colors duration-300">
                      Privacy
                    </a>
                  </motion.li>
                </ul>
              </div>
            </div>
          </div>
        </footer>
      </div>
    </>
  )
}

export default Homepage
