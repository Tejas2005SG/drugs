import React from 'react';
import {
  BeakerIcon,
  CpuChipIcon,
  CloudArrowUpIcon,
  DocumentMagnifyingGlassIcon,
  UsersIcon,
  SparklesIcon
} from '@heroicons/react/24/outline';
import { Link } from 'react-router-dom';
import Navbar from "../../components/Navbar.jsx";

function Homepage() {
  const features = [
    {
      icon: <SparklesIcon className="w-8 h-8" />,
      title: "AI-Powered Molecule Generation",
      description: "Generate novel molecules using SMILES notation with advanced generative AI algorithms.",
      color: "from-purple-500 to-indigo-600"
    },
    {
      icon: <DocumentMagnifyingGlassIcon className="w-8 h-8" />,
      title: "Molecule Visualization & Analysis",
      description: "Visualize and analyze molecules in real-time using RDKit for data-driven insights.",
      color: "from-blue-500 to-cyan-600"
    },
    {
      icon: <SparklesIcon className="w-8 h-8" />,
      title: "Molecular Evolution",
      description: "Iteratively mutate molecules to optimize for stability, solubility, and efficacy.",
      color: "from-green-500 to-emerald-600"
    },
    {
      icon: <BeakerIcon className="w-8 h-8" />,
      title: "Synthesis Cost Estimation",
      description: "Estimate real-world synthesis costs by analyzing lab materials and complexity.",
      color: "from-yellow-500 to-amber-600"
    },
    {
      icon: <DocumentMagnifyingGlassIcon className="w-8 h-8" />,
      title: "Research Paper Generator",
      description: "Automatically generate research papers in IEEE, APA, or Nature journal styles.",
      color: "from-red-500 to-pink-600"
    },
    {
      icon: <DocumentMagnifyingGlassIcon className="w-8 h-8" />,
      title: "Research Summaries",
      description: "Summarize the latest research from PubMed, arXiv, and Google Scholar using Gemini AI.",
      color: "from-indigo-500 to-blue-600"
    },
    {
      icon: <SparklesIcon className="w-8 h-8" />,
      title: "Discovery Recommendations",
      description: "Get AI-driven recommendations for new molecules based on your previous work.",
      color: "from-purple-500 to-pink-600"
    },
    {
      icon: <SparklesIcon className="w-8 h-8" />,
      title: "Drug Naming",
      description: "Get intelligent, systematic names for your novel drug candidates.",
      color: "from-blue-500 to-indigo-600"
    },
    {
      icon: <SparklesIcon className="w-8 h-8" />,
      title: "Toxicity Prediction",
      description: "Predict potential toxicity and side effects in real-time.",
      color: "from-green-500 to-teal-600"
    },
    {
      icon: <BeakerIcon className="w-8 h-8" />,
      title: "Voice-to-Text Notes",
      description: "Dictate observations, and AI will auto-transcribe and summarize them for you.",
      color: "from-yellow-500 to-orange-600"
    },
    {
      icon: <UsersIcon className="w-8 h-8" />,
      title: "Real-Time Collaboration",
      description: "Collaborate seamlessly with your team through shared workspaces.",
      color: "from-red-500 to-purple-600"
    },
    {
      icon: <CloudArrowUpIcon className="w-8 h-8" />,
      title: "AI Chatbot",
      description: "A Gemini-powered assistant to explain complex molecular structures.",
      color: "from-indigo-500 to-purple-600"
    }
  ];

  return (
    <div className="min-h-screen bg-gradient-to-br from-gray-50 to-gray-100">
      <Navbar />
      {/* Hero Section */}
      <section className="relative overflow-hidden bg-gradient-to-r from-blue-600 to-indigo-800 text-white">
        <div className="absolute inset-0 bg-[url('https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?ixlib=rb-1.2.1&auto=format&fit=crop&w=2970&q=80')] opacity-10"></div>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-32 relative">
          <div className="text-center max-w-3xl mx-auto">
            <h1 className="text-5xl md:text-6xl font-extrabold mb-6 leading-tight animate-fade-in">
              <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-200 to-white">Revolutionizing</span> Drug Discovery with AI
            </h1>
            <p className="text-xl mb-10 text-blue-100 max-w-2xl mx-auto">
              Harness the power of generative AI to accelerate your research, reduce costs, and discover breakthrough medications faster than ever before.
            </p>
            <div className="flex flex-col sm:flex-row justify-center gap-4">
              <Link 
                to="/dashboard" 
                className="bg-white text-blue-600 px-8 py-4 rounded-xl font-bold hover:bg-gray-100 transition transform hover:-translate-y-1 shadow-lg hover:shadow-xl"
              >
                Start Your Discovery Journey
              </Link>
              <Link 
                to="#features" 
                className="border-2 border-white text-white px-8 py-4 rounded-xl font-bold hover:bg-white hover:bg-opacity-10 transition transform hover:-translate-y-1"
              >
                Explore Features
              </Link>
            </div>
          </div>
        </div>
        
        {/* Animated floating molecules */}
        <div className="absolute top-0 left-0 w-full h-full overflow-hidden pointer-events-none">
          {[...Array(10)].map((_, i) => (
            <div 
              key={i}
              className="absolute rounded-full bg-white bg-opacity-10"
              style={{
                width: `${Math.random() * 10 + 5}px`,
                height: `${Math.random() * 10 + 5}px`,
                top: `${Math.random() * 100}%`,
                left: `${Math.random() * 100}%`,
                animation: `float ${Math.random() * 15 + 10}s linear infinite`,
                animationDelay: `${Math.random() * 5}s`
              }}
            />
          ))}
        </div>
      </section>

      {/* Features Section */}
      <section id="features" className="py-20 relative">
        {/* Decorative elements */}
        <div className="absolute top-0 left-0 w-full h-32 bg-gradient-to-b from-blue-600/10 to-transparent"></div>
        <div className="absolute bottom-0 left-0 w-full h-32 bg-gradient-to-t from-blue-600/10 to-transparent"></div>
        
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-4xl font-bold text-gray-800 mb-4">
              <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-indigo-600">
                Powerful Features
              </span>
            </h2>
            <p className="text-xl text-gray-600 max-w-3xl mx-auto">
              Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
            </p>
          </div>
          
          {/* Interactive feature grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {features.map((feature, index) => (
              <div 
                key={index}
                className={`bg-white p-8 rounded-2xl shadow-lg hover:shadow-xl transition-all duration-300 transform hover:-translate-y-2 border-l-8 border-${feature.color.split(' ')[1]}`}
              >
                <div className={`inline-flex items-center justify-center w-14 h-14 rounded-full bg-gradient-to-r ${feature.color} mb-6 text-white`}>
                  {feature.icon}
                </div>
                <h3 className="text-xl font-bold text-gray-800 mb-3">{feature.title}</h3>
                <p className="text-gray-600">{feature.description}</p>
              </div>
            ))}
          </div>
          
          {/* Molecular structure animation */}
          <div className="mt-20 flex justify-center">
            <div className="relative w-full max-w-2xl h-64">
              <div className="absolute inset-0 bg-gradient-to-r from-blue-100 to-indigo-100 rounded-3xl opacity-70"></div>
              <div className="absolute inset-0 flex items-center justify-center">
                <div className="relative">
                  {/* Animated molecule structure */}
                  <div className="relative w-48 h-48">
                    {[...Array(6)].map((_, i) => (
                      <div 
                        key={i}
                        className="absolute rounded-full bg-blue-600"
                        style={{
                          width: '12px',
                          height: '12px',
                          top: `${Math.sin(i * Math.PI / 3) * 40 + 40}px`,
                          left: `${Math.cos(i * Math.PI / 3) * 40 + 40}px`,
                          animation: `pulse 2s infinite ${i * 0.2}s`
                        }}
                      />
                    ))}
                    {[...Array(6)].map((_, i) => (
                      <div 
                        key={`line-${i}`}
                        className="absolute bg-blue-400"
                        style={{
                          width: '2px',
                          height: '40px',
                          top: '54px',
                          left: '54px',
                          transform: `rotate(${i * 60}deg)`,
                          transformOrigin: '0 0',
                          animation: `glow 3s infinite ${i * 0.1}s`
                        }}
                      />
                    ))}
                  </div>
                </div>
              </div>
              <div className="absolute bottom-8 left-0 right-0 text-center">
                <p className="text-lg font-medium text-gray-700">Interactive molecular visualization powered by RDKit</p>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Testimonials */}
      <section className="py-20 bg-gradient-to-br from-indigo-50 to-blue-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <h2 className="text-3xl font-bold text-center text-gray-800 mb-12">
            What Researchers Are Saying
          </h2>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            {[
              {
                quote: "This platform reduced our initial discovery phase from 6 months to just 2 weeks. Unbelievable!",
                author: "Dr. Sarah Chen, Lead Researcher",
                affiliation: "Stanford University"
              },
              {
                quote: "The AI-generated molecules had better binding affinity than our manually designed compounds.",
                author: "Prof. Raj Patel",
                affiliation: "MIT Bioengineering"
              },
              {
                quote: "Finally, a tool that bridges the gap between computational chemistry and practical drug development.",
                author: "Dr. Elena Rodriguez",
                affiliation: "Novartis Pharmaceuticals"
              }
            ].map((testimonial, index) => (
              <div key={index} className="bg-white p-8 rounded-2xl shadow-md hover:shadow-lg transition">
                <div className="text-blue-500 text-5xl mb-4">"</div>
                <p className="text-gray-700 text-lg mb-6">{testimonial.quote}</p>
                <div className="flex items-center">
                  <div className="bg-blue-100 w-12 h-12 rounded-full flex items-center justify-center text-blue-600 font-bold mr-4">
                    {testimonial.author.charAt(0)}
                  </div>
                  <div>
                    <h4 className="font-bold text-gray-800">{testimonial.author}</h4>
                    <p className="text-gray-600 text-sm">{testimonial.affiliation}</p>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Call-to-Action Section */}
      <section className="py-24 bg-gradient-to-r from-blue-600 to-indigo-700 text-white relative overflow-hidden">
        {/* Floating particles */}
        <div className="absolute inset-0 overflow-hidden">
          {[...Array(20)].map((_, i) => (
            <div 
              key={i}
              className="absolute rounded-full bg-white opacity-10"
              style={{
                width: `${Math.random() * 20 + 5}px`,
                height: `${Math.random() * 20 + 5}px`,
                top: `${Math.random() * 100}%`,
                left: `${Math.random() * 100}%`,
                animation: `float ${Math.random() * 20 + 10}s linear infinite`,
                animationDelay: `${Math.random() * 5}s`
              }}
            />
          ))}
        </div>
        
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center relative">
          <h2 className="text-4xl font-bold mb-6">
            Ready to Transform <span className="bg-clip-text text-transparent bg-gradient-to-r from-blue-200 to-cyan-200">Drug Discovery</span>?
          </h2>
          <p className="text-xl mb-10 text-blue-100 max-w-3xl mx-auto">
            Join thousands of researchers using our platform to accelerate their drug discovery process with the power of generative AI.
          </p>
          <div className="flex flex-col sm:flex-row justify-center gap-4">
            <Link 
              to="/dashboard" 
              className="bg-white text-blue-600 px-8 py-4 rounded-xl font-bold hover:bg-gray-100 transition transform hover:scale-105 shadow-lg hover:shadow-xl"
            >
              Get Started Now
            </Link>
            <Link 
              to="#contact" 
              className="border-2 border-white text-white px-8 py-4 rounded-xl font-bold hover:bg-white hover:bg-opacity-10 transition transform hover:scale-105"
            >
              Contact Our Team
            </Link>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer id="contact" className="bg-gray-900 text-white py-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-8">
            <div>
              <h3 className="text-xl font-bold mb-4 flex items-center">
                <SparklesIcon className="w-6 h-6 mr-2 text-blue-400" />
                Drug Discovery AI
              </h3>
              <p className="text-gray-400">
                Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more efficient.
              </p>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4">Features</h3>
              <ul className="space-y-2">
                <li><a href="#features" className="text-gray-400 hover:text-white transition">Molecule Generation</a></li>
                <li><a href="#features" className="text-gray-400 hover:text-white transition">Visualization Tools</a></li>
                <li><a href="#features" className="text-gray-400 hover:text-white transition">Research Automation</a></li>
                <li><a href="#features" className="text-gray-400 hover:text-white transition">Collaboration</a></li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4">Company</h3>
              <ul className="space-y-2">
                <li><a href="#" className="text-gray-400 hover:text-white transition">About Us</a></li>
                <li><a href="#" className="text-gray-400 hover:text-white transition">Careers</a></li>
                <li><a href="#" className="text-gray-400 hover:text-white transition">Blog</a></li>
                <li><a href="#" className="text-gray-400 hover:text-white transition">Press</a></li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4">Contact</h3>
              <address className="not-italic text-gray-400 space-y-2">
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
            <div className="flex space-x-6 mt-4 md:mt-0">
              <a href="#" className="text-gray-400 hover:text-white transition">
                <span className="sr-only">Twitter</span>
                <svg className="h-6 w-6" fill="currentColor" viewBox="0 0 24 24">
                  <path d="M8.29 20.251c7.547 0 11.675-6.253 11.675-11.675 0-.178 0-.355-.012-.53A8.348 8.348 0 0022 5.92a8.19 8.19 0 01-2.357.646 4.118 4.118 0 001.804-2.27 8.224 8.224 0 01-2.605.996 4.107 4.107 0 00-6.993 3.743 11.65 11.65 0 01-8.457-4.287 4.106 4.106 0 001.27 5.477A4.072 4.072 0 012.8 9.713v.052a4.105 4.105 0 003.292 4.022 4.095 4.095 0 01-1.853.07 4.108 4.108 0 003.834 2.85A8.233 8.233 0 012 18.407a11.616 11.616 0 006.29 1.84" />
                </svg>
              </a>
              <a href="#" className="text-gray-400 hover:text-white transition">
                <span className="sr-only">LinkedIn</span>
                <svg className="h-6 w-6" fill="currentColor" viewBox="0 0 24 24">
                  <path d="M19 0h-14c-2.761 0-5 2.239-5 5v14c0 2.761 2.239 5 5 5h14c2.762 0 5-2.239 5-5v-14c0-2.761-2.238-5-5-5zm-11 19h-3v-11h3v11zm-1.5-12.268c-.966 0-1.75-.79-1.75-1.764s.784-1.764 1.75-1.764 1.75.79 1.75 1.764-.783 1.764-1.75 1.764zm13.5 12.268h-3v-5.604c0-3.368-4-3.113-4 0v5.604h-3v-11h3v1.765c1.396-2.586 7-2.777 7 2.476v6.759z" />
                </svg>
              </a>
              <a href="#" className="text-gray-400 hover:text-white transition">
                <span className="sr-only">GitHub</span>
                <svg className="h-6 w-6" fill="currentColor" viewBox="0 0 24 24">
                  <path fillRule="evenodd" d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z" clipRule="evenodd" />
                </svg>
              </a>
            </div>
          </div>
        </div>
      </footer>

      {/* Global styles for animations */}
      <style jsx global>{`
        @keyframes float {
          0% {
            transform: translateY(0) rotate(0deg);
          }
          50% {
            transform: translateY(-20px) rotate(5deg);
          }
          100% {
            transform: translateY(0) rotate(0deg);
          }
        }
        @keyframes pulse {
          0% {
            transform: scale(1);
            opacity: 1;
          }
          50% {
            transform: scale(1.1);
            opacity: 0.8;
          }
          100% {
            transform: scale(1);
            opacity: 1;
          }
        }
        @keyframes glow {
          0% {
            opacity: 0.3;
          }
          50% {
            opacity: 1;
          }
          100% {
            opacity: 0.3;
          }
        }
      `}</style>
    </div>
  );
}

export default Homepage;