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
  ChevronDown,
  Play,
  Sparkles,
  Cpu,
  Activity,
  Globe,
} from "lucide-react"

function Homepage() {
  const [activeFeature, setActiveFeature] = useState()
  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 })
  const [counters, setCounters] = useState({ stat1: 0, stat2: 0, stat3: 0, stat4: 0 })
  const [isVisible, setIsVisible] = useState({
    hero: false,
    stats: false,
    features: false,
    testimonials: false,
    cta: false,
  })

  // Refs for scroll animations
  const heroRef = useRef(null)
  const statsRef = useRef(null)
  const featuresRef = useRef(null)
  const testimonialsRef = useRef(null)
  const ctaRef = useRef(null)

  // Intersection Observer for animations
  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const id = entry.target.id
            setIsVisible((prev) => ({ ...prev, [id]: true }))
          }
        })
      },
      { threshold: 0.1 },
    )

    const refs = [heroRef, statsRef, featuresRef, testimonialsRef, ctaRef]
    refs.forEach((ref) => {
      if (ref.current) observer.observe(ref.current)
    })

    return () => observer.disconnect()
  }, [])

  // Mouse tracking for parallax effects
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

  // Counter animations
  useEffect(() => {
    if (isVisible.stats) {
      const animateCounter = (target, key, duration = 2000) => {
        let start = 0
        const increment = target / (duration / 16)
        const timer = setInterval(() => {
          start += increment
          if (start >= target) {
            setCounters((prev) => ({ ...prev, [key]: target }))
            clearInterval(timer)
          } else {
            setCounters((prev) => ({ ...prev, [key]: Math.floor(start) }))
          }
        }, 16)
      }

      animateCounter(85, "stat1")
      animateCounter(60, "stat2")
      animateCounter(3, "stat3")
      setTimeout(() => setCounters((prev) => ({ ...prev, stat4: 24 })), 1000)
    }
  }, [isVisible.stats])

  const features = [
    {
      icon: <Brain className="w-6 h-6" />,
      title: "AI-Powered Molecule Generation",
      description: "Generate novel molecules using SMILES notation with advanced generative AI algorithms.",
      gradient: "from-purple-500 to-blue-500",
    },
    {
      icon: <Microscope className="w-6 h-6" />,
      title: "Molecule Visualization & Analysis",
      description: "Visualize and analyze molecules in real-time using RDKit for data-driven insights.",
      gradient: "from-blue-500 to-cyan-500",
    },
    {
      icon: <Beaker className="w-6 h-6" />,
      title: "Molecular Evolution",
      description: "Iteratively mutate molecules to optimize for stability, solubility, and efficacy.",
      gradient: "from-cyan-500 to-teal-500",
    },
    {
      icon: <BarChart3 className="w-6 h-6" />,
      title: "Synthesis Cost Estimation",
      description: "Estimate real-world synthesis costs by analyzing lab materials and complexity.",
      gradient: "from-teal-500 to-green-500",
    },
    {
      icon: <FileText className="w-6 h-6" />,
      title: "Research Paper Generator",
      description: "Automatically generate research papers in IEEE, APA, or Nature journal styles.",
      gradient: "from-green-500 to-yellow-500",
    },
    {
      icon: <Search className="w-6 h-6" />,
      title: "Research Summaries",
      description: "Summarize the latest research from PubMed, arXiv, and Google Scholar using Gemini AI.",
      gradient: "from-yellow-500 to-orange-500",
    },
    {
      icon: <Database className="w-6 h-6" />,
      title: "Discovery Recommendations",
      description: "Get AI-driven recommendations for new molecules based on your previous work.",
      gradient: "from-orange-500 to-red-500",
    },
    {
      icon: <Pill className="w-6 h-6" />,
      title: "Drug Naming",
      description: "Get intelligent, systematic names for your novel drug candidates.",
      gradient: "from-red-500 to-pink-500",
    },
    {
      icon: <Shield className="w-6 h-6" />,
      title: "Toxicity Prediction",
      description: "Predict potential toxicity and side effects in real-time.",
      gradient: "from-pink-500 to-purple-500",
    },
    {
      icon: <Mic className="w-6 h-6" />,
      title: "Voice-to-Text Notes",
      description: "Dictate observations, and AI will auto-transcribe and summarize them for you.",
      gradient: "from-purple-500 to-indigo-500",
    },
    {
      icon: <Users className="w-6 h-6" />,
      title: "Real-Time Collaboration",
      description: "Collaborate seamlessly with your team through shared workspaces.",
      gradient: "from-indigo-500 to-blue-500",
    },
    {
      icon: <MessageCircle className="w-6 h-6" />,
      title: "AI Chatbot",
      description: "A Gemini-powered assistant to explain complex molecular structures.",
      gradient: "from-blue-500 to-cyan-500",
    },
  ]

  return (
    <>
      <style jsx>{`
        @import url('https://fonts.googleapis.com/css2?family=Exo+2:wght@400;600;700&family=Inter:wght@400;500;600&family=Fira+Code:wght@400;500&display=swap');

        :root {
          --primary: #6A00F4;
          --primary-light: #9A4AFF;
          --primary-dark: #4A00B0;
          --secondary: #00D1FF;
          --secondary-dark: #0085A3;
          --ai-accent: #FF2E88;
          --bg-main: #0A0B0E;
          --bg-card: #161A24;
          --bg-hover: #1E2435;
          --text-primary: #E0E4FF;
          --text-secondary: #A2A8C3;
          --text-disabled: #5A6178;
          --success: #00F0A0;
          --warning: #FF9E3D;
          --error: #FF4D6D;
          --font-heading: 'Exo 2', sans-serif;
          --font-body: 'Inter', sans-serif;
          --font-code: 'Fira Code', monospace;
        }

        * {
          margin: 0;
          padding: 0;
          box-sizing: border-box;
        }

        body {
          font-family: var(--font-body);
          background: var(--bg-main);
          color: var(--text-primary);
          overflow-x: hidden;
        }

        .font-heading {
          font-family: var(--font-heading);
        }

        .font-body {
          font-family: var(--font-body);
        }

        .font-code {
          font-family: var(--font-code);
        }

        @keyframes float {
          0%, 100% { transform: translateY(0px); }
          50% { transform: translateY(-20px); }
        }

        @keyframes pulse-glow {
          0%, 100% { box-shadow: 0 0 20px rgba(106, 0, 244, 0.3); }
          50% { box-shadow: 0 0 40px rgba(106, 0, 244, 0.6), 0 0 60px rgba(106, 0, 244, 0.4); }
        }

        @keyframes rotate-slow {
          from { transform: rotate(0deg); }
          to { transform: rotate(360deg); }
        }

        @keyframes slide-up {
          from { opacity: 0; transform: translateY(50px); }
          to { opacity: 1; transform: translateY(0); }
        }

        @keyframes slide-in-left {
          from { opacity: 0; transform: translateX(-50px); }
          to { opacity: 1; transform: translateX(0); }
        }

        @keyframes slide-in-right {
          from { opacity: 0; transform: translateX(50px); }
          to { opacity: 1; transform: translateX(0); }
        }

        @keyframes scale-in {
          from { opacity: 0; transform: scale(0.8); }
          to { opacity: 1; transform: scale(1); }
        }

        @keyframes gradient-shift {
          0% { background-position: 0% 50%; }
          50% { background-position: 100% 50%; }
          100% { background-position: 0% 50%; }
        }

        @keyframes particle-float {
          0% { transform: translateY(100vh) rotate(0deg); opacity: 0; }
          10% { opacity: 1; }
          90% { opacity: 1; }
          100% { transform: translateY(-100px) rotate(360deg); opacity: 0; }
        }

        .animate-float {
          animation: float 6s ease-in-out infinite;
        }

        .animate-pulse-glow {
          animation: pulse-glow 3s ease-in-out infinite;
        }

        .animate-rotate-slow {
          animation: rotate-slow 20s linear infinite;
        }

        .animate-slide-up {
          animation: slide-up 0.8s ease-out forwards;
        }

        .animate-slide-in-left {
          animation: slide-in-left 0.8s ease-out forwards;
        }

        .animate-slide-in-right {
          animation: slide-in-right 0.8s ease-out forwards;
        }

        .animate-scale-in {
          animation: scale-in 0.6s ease-out forwards;
        }

        .animate-gradient {
          background-size: 200% 200%;
          animation: gradient-shift 4s ease infinite;
        }

        .animate-particle {
          animation: particle-float linear infinite;
        }

        .gradient-text {
          background: linear-gradient(135deg, var(--primary-light), var(--secondary), var(--ai-accent));
          background-clip: text;
          -webkit-background-clip: text;
          -webkit-text-fill-color: transparent;
          background-size: 200% 200%;
          animation: gradient-shift 3s ease infinite;
        }

        .glass-effect {
          background: rgba(22, 26, 36, 0.8);
          backdrop-filter: blur(20px);
          border: 1px solid rgba(160, 168, 195, 0.1);
        }

        .hover-lift {
          transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        .hover-lift:hover {
          transform: translateY(-10px);
          box-shadow: 0 20px 40px rgba(106, 0, 244, 0.2);
        }

        .feature-card {
          transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
          position: relative;
          overflow: hidden;
        }

        .feature-card::before {
          content: '';
          position: absolute;
          top: 0;
          left: -100%;
          width: 100%;
          height: 100%;
          background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.1), transparent);
          transition: left 0.5s;
        }

        .feature-card:hover::before {
          left: 100%;
        }

        .feature-card:hover {
          transform: translateY(-15px) scale(1.02);
          box-shadow: 0 25px 50px rgba(106, 0, 244, 0.3);
        }

        .parallax-bg {
          background: radial-gradient(circle at 20% 80%, rgba(106, 0, 244, 0.3) 0%, transparent 50%),
                      radial-gradient(circle at 80% 20%, rgba(0, 209, 255, 0.3) 0%, transparent 50%),
                      radial-gradient(circle at 40% 40%, rgba(255, 46, 136, 0.2) 0%, transparent 50%);
        }

        .particle {
          position: absolute;
          width: 4px;
          height: 4px;
          background: var(--secondary);
          border-radius: 50%;
          pointer-events: none;
        }

        .delay-100 { animation-delay: 0.1s; }
        .delay-200 { animation-delay: 0.2s; }
        .delay-300 { animation-delay: 0.3s; }
        .delay-400 { animation-delay: 0.4s; }
        .delay-500 { animation-delay: 0.5s; }
        .delay-600 { animation-delay: 0.6s; }
        .delay-700 { animation-delay: 0.7s; }
        .delay-800 { animation-delay: 0.8s; }

        .nav-blur {
          backdrop-filter: blur(20px);
          background: rgba(10, 11, 14, 0.9);
          border-bottom: 1px solid rgba(160, 168, 195, 0.1);
        }

        .btn-primary {
          background: linear-gradient(135deg, var(--primary), var(--primary-light));
          border: none;
          color: white;
          padding: 12px 24px;
          border-radius: 12px;
          font-weight: 600;
          font-family: var(--font-heading);
          cursor: pointer;
          transition: all 0.3s ease;
          position: relative;
          overflow: hidden;
        }

        .btn-primary::before {
          content: '';
          position: absolute;
          top: 0;
          left: -100%;
          width: 100%;
          height: 100%;
          background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.2), transparent);
          transition: left 0.5s;
        }

        .btn-primary:hover::before {
          left: 100%;
        }

        .btn-primary:hover {
          transform: translateY(-2px);
          box-shadow: 0 10px 25px rgba(106, 0, 244, 0.4);
        }

        .btn-secondary {
          background: transparent;
          border: 2px solid var(--secondary);
          color: var(--secondary);
          padding: 10px 22px;
          border-radius: 12px;
          font-weight: 600;
          font-family: var(--font-heading);
          cursor: pointer;
          transition: all 0.3s ease;
        }

        .btn-secondary:hover {
          background: var(--secondary);
          color: var(--bg-main);
          transform: translateY(-2px);
          box-shadow: 0 10px 25px rgba(0, 209, 255, 0.4);
        }

        .text-glow {
          text-shadow: 0 0 20px rgba(106, 0, 244, 0.5);
        }

        .border-glow {
          border: 1px solid rgba(106, 0, 244, 0.3);
          box-shadow: 0 0 20px rgba(106, 0, 244, 0.1);
        }

        .scroll-indicator {
          position: absolute;
          bottom: 30px;
          left: 50%;
          transform: translateX(-50%);
          animation: float 3s ease-in-out infinite;
        }
      `}</style>

      <div className="min-h-screen" style={{ background: "var(--bg-main)" }}>
        {/* Navigation */}
        <nav className="fixed top-0 w-full z-50 nav-blur">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="flex justify-between items-center h-16">
              <div className="flex items-center space-x-3">
                <div className="animate-rotate-slow">
                  <Atom className="w-8 h-8" style={{ color: "var(--primary-light)" }} />
                </div>
                <span className="font-heading font-bold text-xl gradient-text">Drug Discovery AI</span>
              </div>

              <div className="hidden md:flex items-center space-x-8">
                <a
                  href="#features"
                  className="font-body hover:text-white transition-colors duration-300"
                  style={{ color: "var(--text-secondary)" }}
                >
                  Features
                </a>
                <a
                  href="#testimonials"
                  className="font-body hover:text-white transition-colors duration-300"
                  style={{ color: "var(--text-secondary)" }}
                >
                  Testimonials
                </a>
                <a
                  href="#contact"
                  className="font-body hover:text-white transition-colors duration-300"
                  style={{ color: "var(--text-secondary)" }}
                >
                  Contact
                </a>
                <button className="btn-primary">Get Started</button>
              </div>
            </div>
          </div>
        </nav>

        {/* Floating Particles */}
        <div className="fixed inset-0 pointer-events-none overflow-hidden">
          {[...Array(50)].map((_, i) => (
            <div
              key={i}
              className="particle animate-particle"
              style={{
                left: Math.random() * 100 + "%",
                animationDuration: Math.random() * 10 + 10 + "s",
                animationDelay: Math.random() * 10 + "s",
                opacity: Math.random() * 0.5 + 0.2,
              }}
            />
          ))}
        </div>

        {/* Hero Section */}
        <section
          id="hero"
          ref={heroRef}
          className="relative min-h-screen flex items-center justify-center parallax-bg pt-16"
        >
          {/* Floating Background Elements */}
          <div className="absolute inset-0 overflow-hidden">
            <div
              className="absolute top-20 left-10 animate-float delay-100"
              style={{
                transform: `translate(${mousePosition.x * 0.02}px, ${mousePosition.y * 0.02}px)`,
              }}
            >
              <Brain className="w-24 h-24 opacity-20" style={{ color: "var(--primary-light)" }} />
            </div>
            <div
              className="absolute top-40 right-20 animate-float delay-300"
              style={{
                transform: `translate(${mousePosition.x * -0.03}px, ${mousePosition.y * 0.03}px)`,
              }}
            >
              <Beaker className="w-32 h-32 opacity-15" style={{ color: "var(--secondary)" }} />
            </div>
            <div
              className="absolute bottom-40 left-20 animate-float delay-500"
              style={{
                transform: `translate(${mousePosition.x * 0.025}px, ${mousePosition.y * -0.025}px)`,
              }}
            >
              <Microscope className="w-28 h-28 opacity-20" style={{ color: "var(--ai-accent)" }} />
            </div>
          </div>

          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 relative">
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
              {/* Left Content */}
              <div className="space-y-8">
                <div
                  className={`inline-flex items-center px-4 py-2 rounded-full glass-effect ${isVisible.hero ? "animate-scale-in" : "opacity-0"}`}
                >
                  <div
                    className="w-2 h-2 rounded-full animate-pulse-glow mr-3"
                    style={{ background: "var(--success)" }}
                  ></div>
                  <span className="font-code text-sm" style={{ color: "var(--text-secondary)" }}>
                    Next-Gen AI for Drug Discovery
                  </span>
                </div>

                <h1
                  className={`font-heading font-bold text-5xl lg:text-7xl leading-tight ${isVisible.hero ? "animate-slide-up delay-200" : "opacity-0"}`}
                >
                  <span style={{ color: "var(--text-primary)" }}>Revolutionizing</span>
                  <br />
                  <span className="gradient-text text-glow">Drug Discovery</span>
                  <br />
                  <span style={{ color: "var(--text-primary)" }}>with AI</span>
                </h1>

                <p
                  className={`font-body text-xl leading-relaxed ${isVisible.hero ? "animate-slide-up delay-400" : "opacity-0"}`}
                  style={{ color: "var(--text-secondary)" }}
                >
                  Harness the power of generative AI to accelerate your research, reduce costs, and discover
                  breakthrough medications faster than ever before.
                </p>

                <div
                  className={`flex flex-col sm:flex-row gap-4 ${isVisible.hero ? "animate-slide-up delay-600" : "opacity-0"}`}
                >
                  <button className="btn-primary flex items-center justify-center group">
                    <Sparkles className="w-5 h-5 mr-2" />
                    Start Your Discovery Journey
                    <ArrowRight className="w-5 h-5 ml-2 group-hover:translate-x-1 transition-transform" />
                  </button>
                  <button className="btn-secondary flex items-center justify-center group">
                    <Play className="w-5 h-5 mr-2" />
                    Watch Demo
                    <Zap className="w-5 h-5 ml-2 group-hover:rotate-12 transition-transform" />
                  </button>
                </div>

                <div
                  className={`flex items-center space-x-6 ${isVisible.hero ? "animate-slide-up delay-800" : "opacity-0"}`}
                >
                  <div className="flex -space-x-2">
                    {[1, 2, 3, 4].map((i) => (
                      <div
                        key={i}
                        className="w-10 h-10 rounded-full border-2 flex items-center justify-center font-semibold hover-lift cursor-pointer"
                        style={{
                          background: "var(--bg-card)",
                          borderColor: "var(--primary)",
                          color: "var(--text-primary)",
                        }}
                      >
                        {i}
                      </div>
                    ))}
                  </div>
                  <div className="font-body text-sm" style={{ color: "var(--text-secondary)" }}>
                    <span className="font-semibold gradient-text">500+</span> researchers already using our platform
                  </div>
                </div>
              </div>

              {/* Right Content - Dashboard Mockup */}
              <div className={`relative ${isVisible.hero ? "animate-slide-in-right delay-300" : "opacity-0"}`}>
                <div
                  className="glass-effect rounded-2xl p-8 border-glow hover-lift"
                  style={{
                    transform: `translate(${mousePosition.x * 0.01}px, ${mousePosition.y * 0.01}px)`,
                  }}
                >
                  {/* Header */}
                  <div className="flex items-center justify-between mb-6">
                    <div className="flex items-center space-x-2">
                      <div className="w-3 h-3 rounded-full" style={{ background: "var(--error)" }}></div>
                      <div className="w-3 h-3 rounded-full" style={{ background: "var(--warning)" }}></div>
                      <div className="w-3 h-3 rounded-full" style={{ background: "var(--success)" }}></div>
                    </div>
                    <div className="font-code text-sm" style={{ color: "var(--text-secondary)" }}>
                      Drug Discovery Dashboard
                    </div>
                  </div>

                  {/* Molecule Grid */}
                  <div className="grid grid-cols-3 gap-4 mb-6">
                    {[1, 2, 3, 4, 5, 6].map((i) => (
                      <div
                        key={i}
                        className="h-16 rounded-lg flex items-center justify-center hover-lift cursor-pointer"
                        style={{
                          background: `linear-gradient(135deg, rgba(106, 0, 244, 0.2), rgba(0, 209, 255, 0.2))`,
                        }}
                      >
                        <div className="animate-rotate-slow">
                          <Beaker className="w-6 h-6" style={{ color: "var(--secondary)" }} />
                        </div>
                      </div>
                    ))}
                  </div>

                  {/* Progress Bars */}
                  <div className="space-y-4">
                    {[
                      { label: "AI Analysis Progress", value: 87, color: "var(--primary)" },
                      { label: "Synthesis Optimization", value: 64, color: "var(--secondary)" },
                      { label: "Toxicity Prediction", value: 92, color: "var(--success)" },
                    ].map((item, index) => (
                      <div key={index}>
                        <div
                          className="flex items-center justify-between text-sm mb-2"
                          style={{ color: "var(--text-secondary)" }}
                        >
                          <span className="font-body">{item.label}</span>
                          <span className="font-code">{item.value}%</span>
                        </div>
                        <div
                          className="h-2 rounded-full overflow-hidden"
                          style={{ background: "rgba(160, 168, 195, 0.1)" }}
                        >
                          <div
                            className="h-full rounded-full transition-all duration-1000 ease-out"
                            style={{
                              width: isVisible.hero ? `${item.value}%` : "0%",
                              background: item.color,
                              transitionDelay: `${index * 200}ms`,
                            }}
                          ></div>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>

                {/* Floating AI Badge */}
                <div
                  className="absolute -top-4 -right-4 w-24 h-24 rounded-full flex items-center justify-center animate-pulse-glow"
                  style={{
                    background: `linear-gradient(135deg, var(--primary), var(--ai-accent))`,
                    transform: `translate(${mousePosition.x * -0.02}px, ${mousePosition.y * -0.02}px)`,
                  }}
                >
                  <div className="animate-rotate-slow">
                    <Cpu className="w-12 h-12 text-white" />
                  </div>
                </div>

                {/* Floating Stats */}
                <div
                  className="absolute -bottom-4 -left-4 glass-effect rounded-xl p-4 border-glow"
                  style={{
                    transform: `translate(${mousePosition.x * 0.015}px, ${mousePosition.y * 0.015}px)`,
                  }}
                >
                  <div className="flex items-center space-x-3">
                    <div className="animate-float">
                      <TrendingUp className="w-6 h-6" style={{ color: "var(--success)" }} />
                    </div>
                    <div>
                      <div className="font-heading font-bold text-lg gradient-text">3.5x</div>
                      <div className="font-body text-xs" style={{ color: "var(--text-secondary)" }}>
                        Faster Discovery
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </div>

          {/* Scroll Indicator */}
          <div className="scroll-indicator">
            <ChevronDown className="w-6 h-6" style={{ color: "var(--text-secondary)" }} />
          </div>
        </section>

        {/* Stats Section */}
        <section id="stats" ref={statsRef} className="py-20" style={{ background: "var(--bg-card)" }}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="grid grid-cols-2 md:grid-cols-4 gap-8">
              {[
                {
                  key: "stat1",
                  value: counters.stat1,
                  suffix: "%",
                  label: "Faster Discovery",
                  icon: <Zap className="w-8 h-8" />,
                },
                {
                  key: "stat2",
                  value: counters.stat2,
                  suffix: "%",
                  label: "Cost Reduction",
                  icon: <TrendingUp className="w-8 h-8" />,
                },
                {
                  key: "stat3",
                  value: counters.stat3,
                  suffix: "x",
                  label: "More Candidates",
                  icon: <Target className="w-8 h-8" />,
                },
                {
                  key: "stat4",
                  value: counters.stat4,
                  suffix: "/7",
                  label: "AI Assistance",
                  icon: <Activity className="w-8 h-8" />,
                },
              ].map((stat, index) => (
                <div
                  key={index}
                  className={`text-center glass-effect rounded-xl p-6 hover-lift ${isVisible.stats ? `animate-scale-in delay-${index * 100 + 200}` : "opacity-0"}`}
                >
                  <div className="flex justify-center mb-4" style={{ color: "var(--primary-light)" }}>
                    {stat.icon}
                  </div>
                  <div className="font-heading font-bold text-4xl gradient-text mb-2">
                    {stat.value}
                    {stat.suffix}
                  </div>
                  <div className="font-body" style={{ color: "var(--text-secondary)" }}>
                    {stat.label}
                  </div>
                </div>
              ))}
            </div>
          </div>
        </section>

        {/* Features Section */}
        <section id="features" ref={featuresRef} className="py-20">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className={`text-center mb-16 ${isVisible.features ? "animate-slide-up" : "opacity-0"}`}>
              <div className="inline-flex items-center px-4 py-2 rounded-full glass-effect mb-6">
                <Sparkles className="w-4 h-4 mr-2" style={{ color: "var(--ai-accent)" }} />
                <span className="font-code text-sm" style={{ color: "var(--text-secondary)" }}>
                  Powerful Capabilities
                </span>
              </div>
              <h2 className="font-heading font-bold text-4xl mb-4 gradient-text text-glow">Accelerate Your Research</h2>
              <p className="font-body text-xl max-w-3xl mx-auto" style={{ color: "var(--text-secondary)" }}>
                Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
              </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
              {features.map((feature, index) => (
                <div
                  key={index}
                  className={`feature-card glass-effect rounded-xl p-8 border-glow cursor-pointer ${isVisible.features ? `animate-slide-up delay-${index * 100}` : "opacity-0"}`}
                  onMouseEnter={() => setActiveFeature(index)}
                  onMouseLeave={() => setActiveFeature(null)}
                >
                  <div
                    className="inline-flex items-center justify-center w-12 h-12 rounded-lg mb-6 animate-gradient"
                    style={{
                      background: `linear-gradient(135deg, ${feature.gradient.split(" ")[1]}, ${feature.gradient.split(" ")[3]})`,
                    }}
                  >
                    {feature.icon}
                  </div>

                  <h3
                    className="font-heading font-semibold text-xl mb-3"
                    style={{
                      color: activeFeature === index ? "var(--primary-light)" : "var(--text-primary)",
                    }}
                  >
                    {feature.title}
                  </h3>

                  <p className="font-body leading-relaxed mb-4" style={{ color: "var(--text-secondary)" }}>
                    {feature.description}
                  </p>

                  <div
                    className="flex items-center font-body font-medium"
                    style={{
                      color: "var(--secondary)",
                      transform: activeFeature === index ? "translateX(8px)" : "translateX(0)",
                      transition: "transform 0.3s ease",
                    }}
                  >
                    <span className="text-sm">Learn more</span>
                    <ArrowRight className="w-4 h-4 ml-2" />
                  </div>
                </div>
              ))}
            </div>
          </div>
        </section>

        {/* Testimonials Section */}
        <section id="testimonials" ref={testimonialsRef} className="py-20" style={{ background: "var(--bg-card)" }}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className={`text-center mb-16 ${isVisible.testimonials ? "animate-slide-up" : "opacity-0"}`}>
              <div className="inline-flex items-center px-4 py-2 rounded-full glass-effect mb-6">
                <Users className="w-4 h-4 mr-2" style={{ color: "var(--success)" }} />
                <span className="font-code text-sm" style={{ color: "var(--text-secondary)" }}>
                  Success Stories
                </span>
              </div>
              <h2 className="font-heading font-bold text-4xl gradient-text text-glow">What Researchers Are Saying</h2>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
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
                <div
                  key={index}
                  className={`glass-effect rounded-xl p-8 border-glow hover-lift ${isVisible.testimonials ? `animate-scale-in delay-${index * 200}` : "opacity-0"}`}
                >
                  <div className="flex mb-4">
                    {[...Array(testimonial.rating)].map((_, i) => (
                      <Star
                        key={i}
                        className="w-5 h-5 fill-current animate-float"
                        style={{
                          color: "var(--warning)",
                          animationDelay: `${i * 100}ms`,
                        }}
                      />
                    ))}
                  </div>

                  <p className="font-body text-lg mb-6 leading-relaxed" style={{ color: "var(--text-primary)" }}>
                    "{testimonial.quote}"
                  </p>

                  <div className="flex items-center">
                    <div
                      className="w-12 h-12 rounded-full flex items-center justify-center font-heading font-semibold mr-4 hover-lift"
                      style={{ background: `linear-gradient(135deg, var(--primary), var(--secondary))` }}
                    >
                      {testimonial.author.charAt(0)}
                    </div>
                    <div>
                      <h4 className="font-heading font-semibold" style={{ color: "var(--text-primary)" }}>
                        {testimonial.author}
                      </h4>
                      <p className="font-body text-sm" style={{ color: "var(--text-secondary)" }}>
                        {testimonial.role}
                      </p>
                      <p className="font-body text-sm" style={{ color: "var(--text-disabled)" }}>
                        {testimonial.affiliation}
                      </p>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        </section>

        {/* CTA Section */}
        <section
          id="cta"
          ref={ctaRef}
          className="py-20 relative overflow-hidden"
          style={{
            background: `linear-gradient(135deg, var(--primary), var(--ai-accent), var(--secondary))`,
            backgroundSize: "200% 200%",
          }}
        >
          <div className="absolute inset-0 animate-gradient opacity-50"></div>

          <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center relative">
            <div
              className={`inline-flex items-center px-4 py-2 rounded-full glass-effect mb-6 ${isVisible.cta ? "animate-scale-in" : "opacity-0"}`}
            >
              <div
                className="w-2 h-2 rounded-full animate-pulse-glow mr-3"
                style={{ background: "var(--success)" }}
              ></div>
              <span className="font-code text-sm text-white">Limited Time Offer</span>
            </div>

            <h2
              className={`font-heading font-bold text-4xl text-white mb-6 text-glow ${isVisible.cta ? "animate-slide-up delay-200" : "opacity-0"}`}
            >
              Ready to Transform Drug Discovery?
            </h2>

            <p
              className={`font-body text-xl text-white opacity-90 mb-10 max-w-3xl mx-auto ${isVisible.cta ? "animate-slide-up delay-400" : "opacity-0"}`}
            >
              Join thousands of researchers using our platform to accelerate their drug discovery process with the power
              of generative AI.
            </p>

            <div
              className={`flex flex-col sm:flex-row justify-center gap-4 ${isVisible.cta ? "animate-slide-up delay-600" : "opacity-0"}`}
            >
              <button className="bg-white text-black px-8 py-4 rounded-lg font-heading font-semibold hover-lift flex items-center justify-center group">
                <Globe className="w-5 h-5 mr-2" />
                Get Started Now
                <ArrowRight className="w-5 h-5 ml-2 group-hover:translate-x-1 transition-transform" />
              </button>
              <button className="border-2 border-white text-white px-8 py-4 rounded-lg font-heading font-semibold hover-lift flex items-center justify-center group">
                <MessageCircle className="w-5 h-5 mr-2" />
                Contact Our Team
                <Sparkles className="w-5 h-5 ml-2 group-hover:rotate-12 transition-transform" />
              </button>
            </div>
          </div>
        </section>

        {/* Footer */}
        <footer className="py-16" style={{ background: "var(--bg-main)" }}>
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="grid grid-cols-1 md:grid-cols-4 gap-12">
              <div className="space-y-4">
                <div className="flex items-center">
                  <div className="animate-rotate-slow">
                    <Atom className="w-8 h-8 mr-3" style={{ color: "var(--primary-light)" }} />
                  </div>
                  <span className="font-heading font-bold text-xl gradient-text">Drug Discovery AI</span>
                </div>
                <p className="font-body leading-relaxed" style={{ color: "var(--text-secondary)" }}>
                  Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more
                  efficient.
                </p>
              </div>

              <div>
                <h3 className="font-heading font-semibold text-lg mb-4" style={{ color: "var(--text-primary)" }}>
                  Features
                </h3>
                <ul className="space-y-3">
                  {["Molecule Generation", "Visualization Tools", "Research Automation"].map((item, index) => (
                    <li key={index}>
                      <a
                        href="#features"
                        className="font-body hover-lift inline-block transition-all duration-300"
                        style={{ color: "var(--text-secondary)" }}
                        onMouseEnter={(e) => (e.target.style.color = "var(--text-primary)")}
                        onMouseLeave={(e) => (e.target.style.color = "var(--text-secondary)")}
                      >
                        {item}
                      </a>
                    </li>
                  ))}
                </ul>
              </div>

              <div>
                <h3 className="font-heading font-semibold text-lg mb-4" style={{ color: "var(--text-primary)" }}>
                  Company
                </h3>
                <ul className="space-y-3">
                  {["About Us", "Careers", "Blog"].map((item, index) => (
                    <li key={index}>
                      <a
                        href="#"
                        className="font-body hover-lift inline-block transition-all duration-300"
                        style={{ color: "var(--text-secondary)" }}
                        onMouseEnter={(e) => (e.target.style.color = "var(--text-primary)")}
                        onMouseLeave={(e) => (e.target.style.color = "var(--text-secondary)")}
                      >
                        {item}
                      </a>
                    </li>
                  ))}
                </ul>
              </div>

              <div id="contact">
                <h3 className="font-heading font-semibold text-lg mb-4" style={{ color: "var(--text-primary)" }}>
                  Contact
                </h3>
                <div className="space-y-3">
                  {[
                    { icon: <Mail className="w-4 h-4" />, text: "alokchaturvedi190@gmail.com" },
                    { icon: <Phone className="w-4 h-4" />, text: "+91 9975175098" },
                    { icon: <Mail className="w-4 h-4" />, text: "tbhangale9@gmail.com" },
                    { icon: <Phone className="w-4 h-4" />, text: "+91 8766816061" },
                    { icon: <MapPin className="w-4 h-4" />, text: "Pune, Maharashtra, India" },
                  ].map((contact, index) => (
                    <div
                      key={index}
                      className="flex items-center hover-lift cursor-pointer transition-all duration-300"
                      style={{ color: "var(--text-secondary)" }}
                      onMouseEnter={(e) => (e.currentTarget.style.color = "var(--text-primary)")}
                      onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-secondary)")}
                    >
                      <span className="mr-2">{contact.icon}</span>
                      <span className="font-body text-sm">{contact.text}</span>
                    </div>
                  ))}
                </div>
              </div>
            </div>

            <div
              className="mt-12 pt-8 border-t flex flex-col md:flex-row justify-between items-center"
              style={{ borderColor: "rgba(160, 168, 195, 0.1)" }}
            >
              <p className="font-body" style={{ color: "var(--text-secondary)" }}>
                © {new Date().getFullYear()} Drug Discovery AI. All rights reserved.
              </p>
              <div className="mt-4 md:mt-0">
                <ul className="flex space-x-6 text-sm">
                  {["Terms", "Privacy"].map((item, index) => (
                    <li key={index}>
                      <a
                        href="#"
                        className="font-body hover-lift inline-block transition-all duration-300"
                        style={{ color: "var(--text-secondary)" }}
                        onMouseEnter={(e) => (e.target.style.color = "var(--text-primary)")}
                        onMouseLeave={(e) => (e.target.style.color = "var(--text-secondary)")}
                      >
                        {item}
                      </a>
                    </li>
                  ))}
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
