"use client"

import { useState, useEffect } from "react"
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
} from "lucide-react";
import "./homepage.css" // Import custom styles for animations and effects
import Navbar from "../../components/Navbar";
function Homepage() {
  const [activeFeature, setActiveFeature] = useState(null)
  const [isVisible, setIsVisible] = useState({})
  const [counters, setCounters] = useState({ stat1: 0, stat2: 0, stat3: 0, stat4: 0 })
  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 })

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

  // Intersection Observer for scroll animations
  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            setIsVisible((prev) => ({ ...prev, [entry.target.id]: true }))
          }
        })
      },
      { threshold: 0.1 },
    )

    const elements = document.querySelectorAll("[data-animate]")
    elements.forEach((el) => observer.observe(el))

    return () => observer.disconnect()
  }, [])

  // Counter animation for stats
  useEffect(() => {
    const animateCounter = (key, target, duration = 2000) => {
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

    if (isVisible.stats) {
      animateCounter("stat1", 85)
      animateCounter("stat2", 60)
      animateCounter("stat3", 3.5)
      setTimeout(() => setCounters((prev) => ({ ...prev, stat4: "24/7" })), 1000)
    }
  }, [isVisible.stats])

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
    <div className="min-h-screen bg-white overflow-hidden">
     
     
      <Navbar />

      {/* Hero Section */}
      <section className="relative bg-gradient-to-br mt-3 from-blue-50 via-white to-purple-50 py-20 overflow-hidden">
        {/* Floating Background Shapes */}
        <div className="floating-shapes">
          <div className="shape shape-1">
            <Atom className="w-24 h-24 text-blue-300" />
          </div>
          <div className="shape shape-2">
            <Brain className="w-32 h-32 text-purple-300" />
          </div>
          <div className="shape shape-3">
            <Beaker className="w-20 h-20 text-green-300" />
          </div>
          <div className="shape shape-4">
            <Microscope className="w-28 h-28 text-orange-300" />
          </div>
        </div>

        {/* Animated Particles */}
        <div className="absolute inset-0 overflow-hidden pointer-events-none">
          {[...Array(20)].map((_, i) => (
            <div
              key={i}
              className="absolute rounded-full bg-blue-400 opacity-20"
              style={{
                width: Math.random() * 4 + 2 + "px",
                height: Math.random() * 4 + 2 + "px",
                left: Math.random() * 100 + "%",
                top: Math.random() * 100 + "%",
                animationDelay: Math.random() * 5 + "s",
              }}
            />
          ))}
        </div>

        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 relative">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            <div className="space-y-8">
              <div className="opacity-0 animate-slide-up">
                <div className="inline-flex items-center px-4 py-2 rounded-full bg-blue-100 text-blue-800 text-sm font-medium mb-6 hover:bg-blue-200 transition-colors duration-300">
                  <span className="w-2 h-2 bg-blue-600 rounded-full mr-2 animate-pulse"></span>
                  Next-Gen AI for Drug Discovery
                </div>
              </div>

              <h1 className="text-5xl lg:text-6xl font-bold text-gray-900 leading-tight">
                <span className="opacity-0 animate-slide-up delay-100 inline-block">Revolutionizing</span>
                <span className="gradient-text block opacity-0 animate-slide-up delay-200">Drug Discovery</span>
                <span className="opacity-0 animate-slide-up delay-300 inline-block">with AI</span>
              </h1>

              <p className="text-xl text-gray-600 leading-relaxed opacity-0 animate-slide-up delay-400">
                Harness the power of generative AI to accelerate your research, reduce costs, and discover breakthrough
                medications faster than ever before.
              </p>

              <div className="flex flex-col sm:flex-row gap-4 opacity-0 animate-slide-up delay-500">
                <button className="group bg-blue-600 text-white px-8 py-4 rounded-lg hover:bg-blue-700 transition-all duration-300 flex items-center justify-center font-semibold hover:scale-105 hover:shadow-xl">
                  Start Your Discovery Journey
                  <ArrowRight className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-300" />
                </button>
                <button className="group border-2 border-gray-300 text-gray-700 px-8 py-4 rounded-lg hover:border-blue-600 hover:text-blue-600 transition-all duration-300 font-semibold hover:scale-105">
                  Explore Features
                  <Zap className="ml-2 w-5 h-5 group-hover:rotate-12 transition-transform duration-300" />
                </button>
              </div>

              <div className="flex items-center space-x-6 opacity-0 animate-slide-up delay-600">
                <div className="flex -space-x-2">
                  {[1, 2, 3, 4].map((i) => (
                    <div
                      key={i}
                      className="w-10 h-10 rounded-full bg-gray-300 border-2 border-white flex items-center justify-center text-gray-600 font-semibold hover:scale-110 transition-transform duration-300 cursor-pointer"
                      style={{ animationDelay: `${i * 0.1}s` }}
                    >
                      {i}
                    </div>
                  ))}
                </div>
                <div className="text-sm text-gray-600">
                  <span className="font-semibold text-gray-900">500+</span> researchers already using our platform
                </div>
              </div>
            </div>

            <div className="relative opacity-0 animate-slide-in-right delay-300">
              {/* Main Dashboard Mockup */}
              <div
                className="bg-white rounded-2xl shadow-2xl p-8 border border-gray-100 hover-lift parallax-element"
                style={{
                  transform: `translate(${(mousePosition.x - 50) * 0.02}px, ${(mousePosition.y - 50) * 0.02}px)`,
                }}
              >
                {/* Header */}
                <div className="flex items-center justify-between mb-6">
                  <div className="flex items-center space-x-2">
                    <div className="w-3 h-3 bg-red-400 rounded-full"></div>
                    <div className="w-3 h-3 bg-yellow-400 rounded-full"></div>
                    <div className="w-3 h-3 bg-green-400 rounded-full"></div>
                  </div>
                  <div className="text-sm text-gray-500">Drug Discovery Dashboard</div>
                </div>

                {/* Molecule Grid */}
                <div className="grid grid-cols-3 gap-4 mb-6">
                  {[1, 2, 3, 4, 5, 6].map((i) => (
                    <div
                      key={i}
                      className="h-16 bg-gradient-to-br from-blue-50 to-purple-50 rounded-lg flex items-center justify-center hover:scale-105 transition-transform duration-300 animate-bounce-gentle"
                      style={{ animationDelay: `${i * 0.2}s` }}
                    >
                      <Beaker className="w-6 h-6 text-blue-500" />
                    </div>
                  ))}
                </div>

                {/* Progress Bars */}
                <div className="space-y-3">
                  <div className="flex items-center justify-between text-sm text-gray-600">
                    <span>AI Analysis Progress</span>
                    <span>87%</span>
                  </div>
                  <div className="h-2 bg-gray-100 rounded-full overflow-hidden">
                    <div
                      className="h-full bg-gradient-to-r from-blue-500 to-purple-500 rounded-full animate-pulse"
                      style={{ width: "87%" }}
                    ></div>
                  </div>

                  <div className="flex items-center justify-between text-sm text-gray-600">
                    <span>Synthesis Optimization</span>
                    <span>64%</span>
                  </div>
                  <div className="h-2 bg-gray-100 rounded-full overflow-hidden">
                    <div
                      className="h-full bg-gradient-to-r from-green-500 to-teal-500 rounded-full"
                      style={{ width: "64%" }}
                    ></div>
                  </div>

                  <div className="flex items-center justify-between text-sm text-gray-600">
                    <span>Toxicity Prediction</span>
                    <span>92%</span>
                  </div>
                  <div className="h-2 bg-gray-100 rounded-full overflow-hidden">
                    <div
                      className="h-full bg-gradient-to-r from-purple-500 to-pink-500 rounded-full"
                      style={{ width: "92%" }}
                    ></div>
                  </div>
                </div>
              </div>

              {/* Floating AI Badge */}
              <div
                className="absolute -top-4 -right-4 w-24 h-24 bg-gradient-to-br from-blue-600 to-purple-600 rounded-full flex items-center justify-center animate-pulse-glow parallax-element"
                style={{
                  transform: `translate(${(mousePosition.x - 50) * -0.03}px, ${(mousePosition.y - 50) * -0.03}px)`,
                }}
              >
                <Brain className="w-12 h-12 text-white animate-bounce-gentle" />
              </div>

              {/* Floating Stats Cards */}
              <div
                className="absolute -bottom-4 -left-4 bg-white rounded-xl p-4 shadow-lg border border-gray-100 animate-float parallax-element"
                style={{
                  transform: `translate(${(mousePosition.x - 50) * 0.01}px, ${(mousePosition.y - 50) * 0.01}px)`,
                }}
              >
                <div className="flex items-center space-x-2">
                  <TrendingUp className="w-5 h-5 text-green-500" />
                  <div>
                    <div className="text-lg font-bold text-gray-900">3.5x</div>
                    <div className="text-xs text-gray-500">Faster Discovery</div>
                  </div>
                </div>
              </div>

              <div
                className="absolute top-1/2 -right-8 bg-white rounded-xl p-4 shadow-lg border border-gray-100 animate-float parallax-element"
                style={{
                  animationDelay: "1s",
                  transform: `translate(${(mousePosition.x - 50) * -0.02}px, ${(mousePosition.y - 50) * -0.02}px)`,
                }}
              >
                <div className="flex items-center space-x-2">
                  <Target className="w-5 h-5 text-blue-500" />
                  <div>
                    <div className="text-lg font-bold text-gray-900">85%</div>
                    <div className="text-xs text-gray-500">Accuracy</div>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Stats Section */}
      <section className="py-16 bg-white border-t border-gray-100" id="stats" data-animate>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-8">
            {[
              { key: "stat1", value: counters.stat1, suffix: "%", label: "Faster Discovery", color: "text-blue-600" },
              { key: "stat2", value: counters.stat2, suffix: "%", label: "Cost Reduction", color: "text-green-600" },
              { key: "stat3", value: counters.stat3, suffix: "x", label: "More Candidates", color: "text-purple-600" },
              { key: "stat4", value: counters.stat4, suffix: "", label: "AI Assistance", color: "text-orange-600" },
            ].map((stat, index) => (
              <div
                key={index}
                className={`text-center opacity-0 ${isVisible.stats ? "animate-slide-up" : ""}`}
                style={{ animationDelay: `${index * 0.1}s` }}
              >
                <div className={`text-4xl font-bold ${stat.color} mb-2`}>
                  {stat.value}
                  {stat.suffix}
                </div>
                <div className="text-gray-600">{stat.label}</div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section id="features" className="py-20 bg-gray-50" data-animate>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className={`text-center mb-16 opacity-0 ${isVisible.features ? "animate-slide-up" : ""}`}>
            <div className="inline-flex items-center px-4 py-2 rounded-full bg-blue-100 text-blue-800 text-sm font-medium mb-6">
              Powerful Capabilities
            </div>
            <h2 className="text-4xl font-bold text-gray-900 mb-4">Accelerate Your Research</h2>
            <p className="text-xl text-gray-600 max-w-3xl mx-auto">
              Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {features.map((feature, index) => (
              <div
                key={index}
                className={`feature-card bg-white rounded-xl p-8 shadow-sm border border-gray-100 opacity-0 ${isVisible.features ? "animate-slide-up" : ""}`}
                style={{ animationDelay: `${index * 0.1}s` }}
                onMouseEnter={() => setActiveFeature(index)}
                onMouseLeave={() => setActiveFeature(null)}
              >
                <div
                  className={`inline-flex items-center justify-center w-12 h-12 rounded-lg ${feature.color} mb-6 transition-transform duration-300 ${activeFeature === index ? "scale-110 rotate-6" : ""}`}
                >
                  {feature.icon}
                </div>

                <h3 className="text-xl font-semibold text-gray-900 mb-3">{feature.title}</h3>

                <p className="text-gray-600 leading-relaxed mb-4">{feature.description}</p>

                <div
                  className={`flex items-center text-blue-600 font-medium transition-all duration-300 ${activeFeature === index ? "translate-x-2" : ""}`}
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
      <section id="testimonials" className="py-20 bg-white" data-animate>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className={`text-center mb-16 opacity-0 ${isVisible.testimonials ? "animate-slide-up" : ""}`}>
            <div className="inline-flex items-center px-4 py-2 rounded-full bg-green-100 text-green-800 text-sm font-medium mb-6">
              Success Stories
            </div>
            <h2 className="text-4xl font-bold text-gray-900 mb-4">What Researchers Are Saying</h2>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            {[
              {
                quote: "This platform reduced our initial discovery phase from 6 months to just 2 weeks. Unbelievable!",
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
                className={`hover-lift bg-gray-50 rounded-xl p-8 border border-gray-100 opacity-0 ${isVisible.testimonials ? "animate-slide-up" : ""}`}
                style={{ animationDelay: `${index * 0.2}s` }}
              >
                <div className="flex mb-4">
                  {[...Array(testimonial.rating)].map((_, i) => (
                    <Star
                      key={i}
                      className="w-5 h-5 text-yellow-400 fill-current animate-bounce-gentle"
                      style={{ animationDelay: `${i * 0.1}s` }}
                    />
                  ))}
                </div>

                <p className="text-gray-700 text-lg mb-6 leading-relaxed">"{testimonial.quote}"</p>

                <div className="flex items-center">
                  <div className="w-12 h-12 rounded-full bg-blue-600 flex items-center justify-center text-white font-semibold mr-4 hover:scale-110 transition-transform duration-300">
                    {testimonial.author.charAt(0)}
                  </div>
                  <div>
                    <h4 className="font-semibold text-gray-900">{testimonial.author}</h4>
                    <p className="text-gray-600 text-sm">{testimonial.role}</p>
                    <p className="text-gray-500 text-sm">{testimonial.affiliation}</p>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-20 bg-gradient-to-r from-blue-600 to-purple-600 relative overflow-hidden">
        <div className="absolute inset-0 bg-black opacity-10"></div>
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center relative">
          <div className="inline-flex items-center px-4 py-2 rounded-full bg-blue-500 text-blue-100 text-sm font-medium mb-6 animate-bounce-gentle">
            Limited Time Offer
          </div>

          <h2 className="text-4xl font-bold text-white mb-6 animate-slide-up">Ready to Transform Drug Discovery?</h2>

          <p className="text-xl text-blue-100 mb-10 max-w-3xl mx-auto animate-slide-up delay-200">
            Join thousands of researchers using our platform to accelerate their drug discovery process with the power
            of generative AI.
          </p>

          <div className="flex flex-col sm:flex-row justify-center gap-4 animate-slide-up delay-400">
            <button className="group bg-white text-blue-600 px-8 py-4 rounded-lg hover:bg-gray-50 transition-all duration-300 font-semibold flex items-center justify-center hover:scale-105 hover:shadow-xl">
              Get Started Now
              <ArrowRight className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-300" />
            </button>
            <button className="group border-2 border-blue-400 text-white px-8 py-4 rounded-lg hover:bg-blue-500 transition-all duration-300 font-semibold hover:scale-105">
              Contact Our Team
              <MessageCircle className="ml-2 w-5 h-5 group-hover:rotate-12 transition-transform duration-300" />
            </button>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer className="bg-gray-900 text-white py-16">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-12">
            <div className="space-y-4">
              <div className="flex items-center">
                <Atom className="w-8 h-8 text-blue-400 mr-3 animate-rotate-slow" />
                <span className="text-xl font-bold">Drug Discovery AI</span>
              </div>
              <p className="text-gray-400 leading-relaxed">
                Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more
                efficient.
              </p>
            </div>

            <div>
              <h3 className="text-lg font-semibold mb-4">Features</h3>
              <ul className="space-y-3">
                <li>
                  <a
                    href="#features"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    Molecule Generation
                  </a>
                </li>
                <li>
                  <a
                    href="#features"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    Visualization Tools
                  </a>
                </li>
                <li>
                  <a
                    href="#features"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    Research Automation
                  </a>
                </li>
              </ul>
            </div>

            <div>
              <h3 className="text-lg font-semibold mb-4">Company</h3>
              <ul className="space-y-3">
                <li>
                  <a
                    href="#"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    About Us
                  </a>
                </li>
                <li>
                  <a
                    href="#"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    Careers
                  </a>
                </li>
                <li>
                  <a
                    href="#"
                    className="text-gray-400 hover:text-white transition-all duration-300 hover:translate-x-1"
                  >
                    Blog
                  </a>
                </li>
              </ul>
            </div>

            <div id="contact">
              <h3 className="text-lg font-semibold mb-4">Contact</h3>
              <div className="space-y-3 text-gray-400">
                <div className="flex items-center hover:text-white transition-colors duration-300">
                  <Mail className="w-4 h-4 mr-2" />
                  <span>alokchaturvedi190@gmail.com</span>
                </div>
                <div className="flex items-center hover:text-white transition-colors duration-300">
                  <Phone className="w-4 h-4 mr-2" />
                  <span>+91 9975175098</span>
                </div>
                <div className="flex items-center hover:text-white transition-colors duration-300">
                  <Mail className="w-4 h-4 mr-2" />
                  <span>tbhangale9@gmail.com</span>
                </div>
                <div className="flex items-center hover:text-white transition-colors duration-300">
                  <Phone className="w-4 h-4 mr-2" />
                  <span>+91 8766816061</span>
                </div>
                <div className="flex items-center hover:text-white transition-colors duration-300">
                  <MapPin className="w-4 h-4 mr-2" />
                  <span>Pune, Maharashtra, India</span>
                </div>
              </div>
            </div>
          </div>

          <div className="mt-12 pt-8 border-t border-gray-800 flex flex-col md:flex-row justify-between items-center">
            <p className="text-gray-400">© {new Date().getFullYear()} Drug Discovery AI. All rights reserved.</p>
            <div className="mt-4 md:mt-0">
              <ul className="flex space-x-6 text-sm">
                <li>
                  <a href="#" className="text-gray-400 hover:text-white transition-colors duration-300">
                    Terms
                  </a>
                </li>
                <li>
                  <a href="#" className="text-gray-400 hover:text-white transition-colors duration-300">
                    Privacy
                  </a>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </footer>
    </div>
  )
}

export default Homepage
