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
import Navbar from "../../components/Navbar.jsx"
function Homepage() {
  return (
    <div className="min-h-screen bg-gray-50">


     
      
      <section className="relative overflow-hidden bg-gradient-to-r from-blue-600 to-indigo-700 text-white">
        <div className="absolute inset-0 bg-[url('https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?ixlib=rb-1.2.1&auto=format&fit=crop&w=2970&q=80')] opacity-10"></div>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-24 relative">
          <div className="text-center max-w-3xl mx-auto">
            <h1 className="text-5xl md:text-6xl font-extrabold mb-6 leading-tight animate-fade-in">
              Revolutionizing Drug Discovery with AI
            </h1>
            <p className="text-xl mb-10 text-blue-100">
              Harness the power of generative AI to accelerate your research, reduce costs, and discover breakthrough medications faster than ever before.
            </p>
            <div className="space-x-4">
            <Link to="/dashboard" className="bg-white text-blue-600 px-6 py-3 rounded-lg font-semibold hover:bg-gray-100 transition">
            Start Your Discovery Journey
          </Link>
            </div>
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section id="features" className="py-16">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <h2 className="text-3xl font-bold text-center text-gray-800 mb-12">Powerful Features for Drug Discovery</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {/* Feature Card 1: AI-Powered Molecule Generation */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <BeakerIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">AI-Powered Molecule Generation</h3>
              <p className="text-gray-600">
                Generate novel molecules using SMILES notation with advanced generative AI algorithms.
              </p>
            </div>

            {/* Feature Card 3: Molecule Visualization */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <DocumentMagnifyingGlassIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Molecule Visualization & Analysis</h3>
              <p className="text-gray-600">
                Visualize and analyze molecules in real-time using RDKit for data-driven insights.
              </p>
            </div>

            {/* Feature Card 4: AI-Powered Molecular Evolution */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <SparklesIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">AI-Powered Molecular Evolution</h3>
              <p className="text-gray-600">
                Iteratively mutate molecules to optimize for stability, solubility, and efficacy.
              </p>
            </div>


            {/* Feature Card 15: Molecular Synthesis Cost Estimation */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <BeakerIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Synthesis Cost Estimation</h3>
              <p className="text-gray-600">
                Estimate real-world synthesis costs by analyzing lab materials and complexity.
              </p>
            </div>


            {/* Feature Card 9: AI-Powered Research Paper Generator */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <DocumentMagnifyingGlassIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">AI-Powered Research Paper Generator</h3>
              <p className="text-gray-600">
                Automatically generate research papers in IEEE, APA, or Nature journal styles.
              </p>
            </div>



            {/* Feature Card 5: Auto-Generated Research Summaries */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <DocumentMagnifyingGlassIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Auto-Generated Research Summaries</h3>
              <p className="text-gray-600">
                Summarize the latest research from PubMed, arXiv, and Google Scholar using Gemini AI.
              </p>
            </div>

            {/* Feature Card 16: Personalized Drug Discovery Recommendations */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <SparklesIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Drug Discovery Recommendations</h3>
              <p className="text-gray-600">
                Get AI-driven recommendations for new molecules based on your previous work.
              </p>
            </div>

            {/* Feature Card 13: AI-Powered Drug Naming Suggestions */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <SparklesIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">AI-Powered Drug Naming</h3>
              <p className="text-gray-600">
                Get intelligent, systematic names for your novel drug candidates.
              </p>
            </div>


            {/* Feature Card 12: Real-Time Toxicity Prediction */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <SparklesIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Real-Time Toxicity Prediction</h3>
              <p className="text-gray-600">
                Predict potential toxicity and side effects using Google ML as molecules are generated.
              </p>
            </div>

            {/* Feature Card 8: Voice-to-Text AI for Research Notes */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <BeakerIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Voice-to-Text Research Notes</h3>
              <p className="text-gray-600">
                Dictate observations, and AI will auto-transcribe and summarize them for you.
              </p>
            </div>

            {/* Feature Card 2: Real-Time Collaboration */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <UsersIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Real-Time Collaboration</h3>
              <p className="text-gray-600">
                Collaborate seamlessly with your team through group messaging and shared workspaces.
              </p>
            </div>






            {/* Feature Card 6: Conversational AI Chatbot */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <CloudArrowUpIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Conversational AI Chatbot</h3>
              <p className="text-gray-600">
                A Gemini-powered assistant to explain complex molecular structures and predictions.
              </p>
            </div>

            {/* Feature Card 7: Quantum Computing-Assisted Drug Discovery
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <CpuChipIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Quantum Computing-Assisted Discovery</h3>
              <p className="text-gray-600">
                Leverage Google’s Quantum AI for ultra-precise molecular simulations.
              </p>
            </div> */}




            {/* Feature Card 10: Auto-Suggested Research Collaborators */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <UsersIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Auto-Suggested Collaborators</h3>
              <p className="text-gray-600">
                Find potential collaborators worldwide based on research similarity.
              </p>
            </div>

            {/* Feature Card 11: AI-Powered Smart Citations */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <DocumentMagnifyingGlassIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">AI-Powered Smart Citations</h3>
              <p className="text-gray-600">
                Auto-generate citations from Google Scholar and arXiv for your findings.
              </p>
            </div>




            {/* Feature Card 14: Live AI-Based Drug News Updates */}
            <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition">
              <DocumentMagnifyingGlassIcon className="text-blue-600 w-10 h-10 mb-4" />
              <h3 className="text-xl font-semibold mb-2">Live Drug News Updates</h3>
              <p className="text-gray-600">
                Stay updated with breaking news on drug trials and FDA approvals from Google News & PubMed.
              </p>
            </div>



          </div>

        </div>
      </section>

      {/* Call-to-Action Section */}
      <section className="bg-blue-600 text-white py-16">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl font-bold mb-4">Ready to Transform Drug Discovery?</h2>
          <p className="text-lg mb-8">
            Join thousands of researchers using our platform to accelerate their drug discovery process with the power of generative AI.
          </p>
          <Link to="/dashboard" className="bg-white text-blue-600 px-6 py-3 rounded-lg font-semibold hover:bg-gray-100 transition">
            Sign Up Now
          </Link>
        </div>
      </section>

      {/* Footer */}
      <footer id="contact" className="bg-gray-800 text-white py-8">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div>
              <h3 className="text-lg font-semibold mb-4">Drug Discovery Assistant</h3>
              <p className="text-gray-400">
                Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more efficient.
              </p>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4">Quick Links</h3>
              <ul className="space-y-2">
                <li><a href="#features" className="text-gray-400 hover:text-white">Features</a></li>
                {/* <li><a href="#about" className="text-gray-400 hover:text-white">About</a></li> */}
                <li><a href="#contact" className="text-gray-400 hover:text-white">Contact</a></li>
              </ul>
            </div>
            <div>
              <h3 className="text-lg font-semibold mb-4">Contact Us</h3>

              <p className="text-gray-400">Email: alokchaturvedi190@gmail.com</p>
              <p className="text-gray-400">Phone: +91 9975175098</p>
            
              <p className="text-gray-400 mt-2">Email: tbhangale9@gmail.com</p>
              <p className="text-gray-400">Phone: +91 8766816061</p>
              
              <p className="text-gray-400 mt-2">Address: Pune, Maharasthra, India</p>
            </div>
          </div>
          <div className="mt-8 text-center text-gray-400">
            © {new Date().getFullYear()} Drug Discovery Assistant. All rights reserved.
          </div>
        </div>
      </footer>
    </div>
  );
}

export default Homepage;