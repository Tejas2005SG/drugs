"use client"

import React, { useState } from "react"
import { motion } from "framer-motion"
import { useInView } from "react-intersection-observer"
import {
  User,
UserPlus,
ArrowRight,
Atom,
} from "lucide-react";
import { Link } from "react-router-dom"
import Navbar from "../../components/Navbar.jsx"

function Homepage() {
  const [activeFeature, setActiveFeature] = useState(null)
  const [ref, inView] = useInView({
    triggerOnce: true,
    threshold: 0.1,
  })
  const features = [
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "AI-Powered Molecule Generation",
      description: "Generate novel molecules using SMILES notation with advanced generative AI algorithms.",
      color: "from-purple-500 to-indigo-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Molecule Visualization & Analysis",
      description: "Visualize and analyze molecules in real-time using RDKit for data-driven insights.",
      color: "from-blue-500 to-cyan-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Molecular Evolution",
      description: "Iteratively mutate molecules to optimize for stability, solubility, and efficacy.",
      color: "from-green-500 to-emerald-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Synthesis Cost Estimation",
      description: "Estimate real-world synthesis costs by analyzing lab materials and complexity.",
      color: "from-yellow-500 to-amber-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Research Paper Generator",
      description: "Automatically generate research papers in IEEE, APA, or Nature journal styles.",
      color: "from-red-500 to-pink-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Research Summaries",
      description: "Summarize the latest research from PubMed, arXiv, and Google Scholar using Gemini AI.",
      color: "from-indigo-500 to-blue-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Discovery Recommendations",
      description: "Get AI-driven recommendations for new molecules based on your previous work.",
      color: "from-purple-500 to-pink-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Drug Naming",
      description: "Get intelligent, systematic names for your novel drug candidates.",
      color: "from-blue-500 to-indigo-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Toxicity Prediction",
      description: "Predict potential toxicity and side effects in real-time.",
      color: "from-green-500 to-teal-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Voice-to-Text Notes",
      description: "Dictate observations, and AI will auto-transcribe and summarize them for you.",
      color: "from-yellow-500 to-orange-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "Real-Time Collaboration",
      description: "Collaborate seamlessly with your team through shared workspaces.",
      color: "from-red-500 to-purple-600",
    },
    {
      icon: <UserPlus className="w-8 h-8" />,
      title: "AI Chatbot",
      description: "A Gemini-powered assistant to explain complex molecular structures.",
      color: "from-indigo-500 to-purple-600",
    },
  ];

  const containerVariants = {
    hidden: { opacity: 0, y: 20 },
    visible: {
      opacity: 1,
      y: 0,
      transition: {
        duration: 0.6,
        staggerChildren: 0.1,
      },
    },
  }

  const itemVariants = {
    hidden: { opacity: 0, y: 20 },
    visible: {
      opacity: 1,
      y: 0,
      transition: { duration: 0.6 },
    },
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-gray-900 via-gray-800 to-gray-900 text-white">
      <Navbar />
      
      {/* Hero Section */}
      <motion.section 
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 1 }}
        className="relative overflow-hidden py-24 sm:py-32"
      >
        {/* Background Image with Overlay */}
        <div className="absolute inset-0 bg-[url('https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?ixlib=rb-1.2.1&auto=format&fit=crop&w=2970&q=80')] bg-cover bg-center opacity-10"></div>
        
        {/* Gradient Overlay */}
        <div className="absolute inset-0 bg-gradient-to-r from-blue-600/20 via-purple-600/20 to-pink-600/20"></div>
        
        {/* Floating Particles */}
        <div className="absolute inset-0 overflow-hidden pointer-events-none">
          {[...Array(30)].map((_, i) => (
            <motion.div
              key={i}
              className="absolute rounded-full"
              style={{
                width: Math.random() * 6 + 2 + 'px',
                height: Math.random() * 6 + 2 + 'px',
                background: `rgba(${Math.floor(Math.random() * 100 + 155)}, ${Math.floor(Math.random() * 100 + 155)}, ${Math.floor(Math.random() * 255)}, ${Math.random() * 0.5 + 0.2})`,
                boxShadow: `0 0 ${Math.random() * 10 + 5}px rgba(${Math.floor(Math.random() * 100 + 155)}, ${Math.floor(Math.random() * 100 + 155)}, ${Math.floor(Math.random() * 255)}, 0.5)`,
                left: `${Math.random() * 100}%`,
                top: `${Math.random() * 100}%`,
              }}
              animate={{
                y: [0, -Math.random() * 30 - 10, 0],
                x: [0, Math.random() * 20 - 10, 0],
                opacity: [0.2, 0.6, 0.2]
              }}
              transition={{
                duration: Math.random() * 5 + 3,
                repeat: Infinity,
                ease: "easeInOut",
                delay: Math.random() * 2
              }}
            />
          ))}
        </div>
        
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 relative">
          <motion.div 
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8 }}
            className="text-center md:text-left md:max-w-2xl"
          >
            <motion.div
              initial={{ opacity: 0, scale: 0.9 }}
              animate={{ opacity: 1, scale: 1 }}
              transition={{ delay: 0.1, duration: 0.5 }}
              className="inline-block mb-4 px-4 py-1 rounded-full bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 backdrop-blur-sm"
            >
              <span className="text-blue-400 font-medium">Next-Gen AI for Drug Discovery</span>
            </motion.div>
            
            <motion.h1 
              className="text-5xl md:text-7xl font-extrabold mb-6 leading-tight"
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: 0.2, duration: 0.8 }}
            >
              <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-400 via-purple-400 to-pink-400">
                Revolutionizing
              </span>
              <br />
              Drug Discovery with AI
            </motion.h1>
            
            <motion.p 
              className="text-xl md:text-2xl mb-10 text-gray-300 max-w-3xl"
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              transition={{ delay: 0.4, duration: 0.8 }}
            >
              Harness the power of generative AI to accelerate your research, reduce costs, and discover breakthrough medications faster than ever before.
            </motion.p>
            
            <motion.div 
              className="flex flex-col sm:flex-row justify-center md:justify-start gap-4"
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: 0.6, duration: 0.8 }}
            >
              <Link 
                to="/dashboard" 
                className="group relative px-8 py-4 overflow-hidden rounded-xl bg-gradient-to-r from-blue-600 to-purple-600 text-white transition-all duration-300 ease-out hover:scale-105 hover:shadow-lg hover:shadow-blue-500/25"
              >
                <span className="absolute inset-0 bg-gradient-to-r from-blue-700 to-purple-700 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></span>
                <span className="relative font-bold flex items-center justify-center">
                  Start Your Discovery Journey
                  <UserPlus className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-200" />
                </span>
              </Link>
              
              <Link 
                to="#features" 
                className="group relative px-8 py-4 overflow-hidden rounded-xl border-2 border-white/30 transition-all duration-300 ease-out hover:border-white/60 hover:scale-105 backdrop-blur-sm"
              >
                <span className="absolute inset-0 bg-white/5 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></span>
                <span className="relative font-bold flex items-center justify-center">
                  Explore Features
                  <UserPlus  className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-200" />
                </span>
              </Link>
            </motion.div>
            
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              transition={{ delay: 0.8, duration: 0.8 }}
              className="mt-12 flex items-center justify-center md:justify-start space-x-6"
            >
              <div className="flex -space-x-2">
                {[...Array(4)].map((_, i) => (
                  <div key={i} className="w-10 h-10 rounded-full border-2 border-gray-900 overflow-hidden bg-gray-700 flex items-center justify-center text-white">
                    {i+1}
                  </div>
                ))}
              </div>
              <div className="text-sm text-gray-300">
                <span className="font-semibold text-white">500+</span> researchers already using our platform
              </div>
            </motion.div>
          </motion.div>
        </div>
      </motion.section>

      {/* Stats Section */}
      <section className="py-12 bg-gradient-to-r from-gray-800/50 to-gray-900/50 border-y border-white/5">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-8">
            {[
              { value: "85%", label: "Faster Discovery" },
              { value: "60%", label: "Cost Reduction" },
              { value: "3.5x", label: "More Candidates" },
              { value: "24/7", label: "AI Assistance" }
            ].map((stat, index) => (
              <motion.div
                key={index}
                initial={{ opacity: 0, y: 20 }}
                whileInView={{ opacity: 1, y: 0 }}
                viewport={{ once: true }}
                transition={{ delay: index * 0.1, duration: 0.5 }}
                className="text-center"
              >
                <div className="text-3xl md:text-4xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-purple-400">
                  {stat.value}
                </div>
                <div className="text-gray-400 mt-1">{stat.label}</div>
              </motion.div>
            ))}
          </div>
        </div>
      </section>

      {/* Features Section */}
      <motion.section
        id="features"
        ref={ref}
        initial="hidden"
        animate={inView ? "visible" : "hidden"}
        variants={containerVariants}
        className="py-24 relative"
      >
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 relative">
          <motion.div 
            className="text-center mb-16"
            variants={itemVariants}
          >
            <motion.div
              initial={{ opacity: 0, scale: 0.9 }}
              whileInView={{ opacity: 1, scale: 1 }}
              viewport={{ once: true }}
              transition={{ duration: 0.5 }}
              className="inline-block mb-4 px-4 py-1 rounded-full bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 backdrop-blur-sm"
            >
              <span className="text-blue-400 font-medium">Powerful Capabilities</span>
            </motion.div>
            
            <h2 className="text-4xl md:text-5xl font-bold mb-6">
              <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-400 via-purple-400 to-pink-400">
                Accelerate Your Research
              </span>
            </h2>
            <p className="text-xl text-gray-300 max-w-3xl mx-auto">
              Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
            </p>
          </motion.div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {features.map((feature, index) => (
              <motion.div
                key={index}
                variants={itemVariants}
                whileHover={{ scale: 1.03, y: -5 }}
                className="group relative"
              >
                <div className={`
                  absolute inset-0 rounded-2xl bg-gradient-to-r ${feature.color} opacity-0 
                  group-hover:opacity-10 transition-all duration-500 blur-xl
                `} />
                
                <div className="relative backdrop-blur-sm p-8 rounded-2xl border border-white/10 hover:border-white/20 transition-all duration-300 h-full flex flex-col">
                  <motion.div
                    whileHover={{ scale: 1.1, rotate: 5 }}
                    className={`
                      inline-flex items-center justify-center w-14 h-14 rounded-xl 
                      bg-gradient-to-r ${feature.color} mb-6 text-white shadow-lg
                    `}
                  >
                    {feature.icon}
                  </motion.div>
                  
                  <h3 className="text-xl font-bold text-white mb-3 group-hover:text-transparent group-hover:bg-clip-text group-hover:bg-gradient-to-r group-hover:from-blue-400 group-hover:to-purple-400 transition-all duration-300">
                    {feature.title}
                  </h3>
                  
                  <p className="text-gray-300 group-hover:text-gray-200 transition-colors duration-300 flex-grow">
                    {feature.description}
                  </p>

                  <div className="mt-6 pt-4 border-t border-white/5 flex justify-between items-center opacity-0 group-hover:opacity-100 transition-opacity duration-300">
                    <span className="text-sm text-blue-400">Learn more</span>
                    <UserPlus className="w-5 h-5 text-blue-400 group-hover:translate-x-1 transition-transform duration-200" />
                  </div>
                </div>
              </motion.div>
            ))}
          </div>
        </div>
      </motion.section>

      {/* Testimonials Section */}
      <motion.section 
        className="py-24 bg-gradient-to-br from-gray-800/50 to-gray-900/50 relative"
        initial={{ opacity: 0 }}
        whileInView={{ opacity: 1 }}
        viewport={{ once: true }}
      >
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            whileInView={{ opacity: 1, scale: 1 }}
            viewport={{ once: true }}
            transition={{ duration: 0.5 }}
            className="text-center mb-16"
          >
            <motion.div
              initial={{ opacity: 0, scale: 0.9 }}
              whileInView={{ opacity: 1, scale: 1 }}
              viewport={{ once: true }}
              transition={{ duration: 0.5 }}
              className="inline-block mb-4 px-4 py-1 rounded-full bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 backdrop-blur-sm"
            >
              <span className="text-blue-400 font-medium">Success Stories</span>
            </motion.div>
            
            <h2 className="text-3xl md:text-4xl font-bold text-center mb-12 bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-purple-400">
              What Researchers Are Saying
            </h2>
          </motion.div>
          
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            {[
              {
                quote: "This platform reduced our initial discovery phase from 6 months to just 2 weeks. Unbelievable!",
                author: "Dr. Sarah Chen",
                role: "Lead Researcher",
                affiliation: "Stanford University",
              },
              {
                quote: "The AI-generated molecules had better binding affinity than our manually designed compounds.",
                author: "Prof. Raj Patel",
                role: "Principal Investigator",
                affiliation: "MIT Bioengineering",
              },
              {
                quote: "Finally, a tool that bridges the gap between computational chemistry and practical drug development.",
                author: "Dr. Elena Rodriguez",
                role: "Senior Scientist",
                affiliation: "Novartis Pharmaceuticals",
              }
            ].map((testimonial, index) => (
              <motion.div
                key={index}
                initial={{ opacity: 0, y: 20 }}
                whileInView={{ opacity: 1, y: 0 }}
                viewport={{ once: true }}
                transition={{ delay: index * 0.1, duration: 0.5 }}
                whileHover={{ y: -5 }}
                className="relative group"
              >
                <div className="absolute inset-0 bg-gradient-to-r from-blue-600/10 to-purple-600/10 rounded-2xl opacity-0 group-hover:opacity-100 transition-opacity duration-300" />
                
                <div className="relative p-8 rounded-2xl backdrop-blur-sm border border-white/10 hover:border-white/20 transition-all duration-300 h-full flex flex-col">
                  <div className="text-blue-400 text-5xl mb-4">"</div>
                  <p className="text-gray-300 text-lg mb-6 flex-grow">{testimonial.quote}</p>
                  <div className="flex items-center">
                    <div className="w-12 h-12 rounded-full bg-blue-600 flex items-center justify-center text-white font-bold mr-4">
                      {testimonial.author.charAt(0)}
                    </div>
                    <div>
                      <h4 className="font-bold text-white">{testimonial.author}</h4>
                      <p className="text-gray-400 text-sm">{testimonial.role}</p>
                      <p className="text-gray-400 text-sm">{testimonial.affiliation}</p>
                    </div>
                  </div>
                </div>
              </motion.div>
            ))}
          </div>
        </div>
      </motion.section>

      {/* CTA Section */}
      <motion.section 
        className="py-32 relative overflow-hidden"
        initial={{ opacity: 0 }}
        whileInView={{ opacity: 1 }}
        viewport={{ once: true }}
      >
        <div className="max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 relative">
          <div className="bg-gradient-to-r from-gray-900/80 to-gray-800/80 backdrop-blur-lg rounded-3xl p-12 border border-white/10 shadow-2xl">
            <motion.div
              initial={{ opacity: 0, y: 20 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true }}
              className="text-center"
            >
              <motion.div
                initial={{ opacity: 0, scale: 0.9 }}
                whileInView={{ opacity: 1, scale: 1 }}
                viewport={{ once: true }}
                transition={{ duration: 0.5 }}
                className="inline-block mb-4 px-4 py-1 rounded-full bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 backdrop-blur-sm"
              >
                <span className="text-blue-400 font-medium">Limited Time Offer</span>
              </motion.div>
              
              <h2 className="text-4xl md:text-5xl font-bold mb-6">
                Ready to Transform{' '}
                <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-purple-400">
                  Drug Discovery
                </span>?
              </h2>
              
              <p className="text-xl text-gray-300 mb-10 max-w-3xl mx-auto">
                Join thousands of researchers using our platform to accelerate their drug discovery process with the power of generative AI.
              </p>
              
              <div className="flex flex-col sm:flex-row justify-center gap-4">
                <Link 
                  to="/dashboard" 
                  className="group relative px-8 py-4 overflow-hidden rounded-xl bg-gradient-to-r from-blue-600 to-purple-600 text-white transition-all duration-300 ease-out hover:scale-105 hover:shadow-lg hover:shadow-blue-500/25"
                >
                  <span className="absolute inset-0 bg-gradient-to-r from-blue-700 to-purple-700 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></span>
                  <span className="relative font-bold flex items-center justify-center">
                    Get Started Now
                    <UserPlus className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-200" />
                  </span>
                </Link>
                
                <Link 
                  to="#contact" 
                  className="group relative px-8 py-4 overflow-hidden rounded-xl border-2 border-white/30 transition-all duration-300 ease-out hover:border-white/60 hover:scale-105 backdrop-blur-sm"
                >
                  <span className="absolute inset-0 bg-white/5 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></span>
                  <span className="relative font-bold flex items-center justify-center">
                    Contact Our Team
                    <UserPlus className="ml-2 w-5 h-5 group-hover:translate-x-1 transition-transform duration-200" />
                  </span>
                </Link>
              </div>
            </motion.div>
          </div>
        </div>
      </motion.section>

      {/* Footer */}
      <footer className="bg-gray-900 text-white py-16 border-t border-white/10">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-12">
            <div className="space-y-4">
              <h3 className="text-xl font-bold mb-4 flex items-center">
                <Atom className="w-6 h-6 mr-2 text-blue-400" />
                <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-purple-400">
                  Drug Discovery AI
                </span>
              </h3>
              <p className="text-gray-400">
                Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more efficient.
              </p>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4 text-white">Features</h3>
              <ul className="space-y-3">
                <li><a href="#features" className="text-gray-400 hover:text-blue-400 transition">Molecule Generation</a></li>
                <li><a href="#features" className="text-gray-400 hover:text-blue-400 transition">Visualization Tools</a></li>
                <li><a href="#features" className="text-gray-400 hover:text-blue-400 transition">Research Automation</a></li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4 text-white">Company</h3>
              <ul className="space-y-3">
                <li><a href="#" className="text-gray-400 hover:text-blue-400 transition">About Us</a></li>
                <li><a href="#" className="text-gray-400 hover:text-blue-400 transition">Careers</a></li>
                <li><a href="#" className="text-gray-400 hover:text-blue-400 transition">Blog</a></li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4 text-white">Contact</h3>
              <address className="not-italic text-gray-400 space-y-3">
                <p>alokchaturvedi190@gmail.com</p>
                <p>+91 9975175098</p>
                <p>tbhangale9@gmail.com</p>
                <p>+91 8766816061</p>
                <p>Pune, Maharashtra, India</p>
              </address>
            </div>
          </div>
          
          <div className="mt-12 pt-8 border-t border-gray-800 flex flex-col md:flex-row justify-between items-center">
            <p className="text-gray-400">
              © {new Date().getFullYear()} Drug Discovery AI. All rights reserved.
            </p>
            <div className="mt-4 md:mt-0">
              <ul className="flex space-x-6 text-sm">
                <li><a href="#" className="text-gray-400 hover:text-blue-400 transition">Terms</a></li>
                <li><a href="#" className="text-gray-400 hover:text-blue-400 transition">Privacy</a></li>
              </ul>
            </div>
          </div>
        </div>
      </footer>
    </div>
  )
}

export default Homepage