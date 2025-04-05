// App.jsx
import React from 'react';
import { Routes, Route, Navigate } from 'react-router-dom';
import { Toaster } from 'react-hot-toast';
import { useAuthStore } from './Store/auth.store.js';

import Navbar from './components/Navbar.jsx';
import Homepage from './pages/Homepage/Homepage.jsx';
import Signup from './pages/SIgnup/Signuppage.jsx';
import Login from './pages/Login/Loginpage.jsx';
import VerifyPhone from './pages/VerifyPhone/VerifyPhone.jsx'
import Dashboard from './components/Dashboard.jsx';
import DashboardHome from './components/Dashboardhome.jsx';
import ProteinStructureEvolution from './pages/Protienstructureevolution/Protienstructureevolution.jsx';
import ProteinStructureApp from './pages/Proteinstructureapp/ProteinStructureApp.jsx';
import AIResearchPaperGenerator from "./pages/AIresearchgenerator/Airesearchgenerator.jsx";
import AIDrivenTargetPrediction from "./pages/AIdriventargetprediction/AIdriventargetprediction.jsx";
import Costestimation from "./pages/Costestimination/Costestimination.jsx";
import DrugDiscoveryRecommendation from "./pages/Drugdiscoveryrecommendation/Drugdiscoveryrecommendation.jsx";
import LiveNews from "./pages/Livenews/Livenews.jsx";
// import Message from "./pages/Message/Message.jsx";
import GetAlphaFoldStructure from "./pages/Alphafold/Alphafold.jsx"
import AINamingSuggestion from './pages/AINamingSuggestion/AINamingSuggestion.jsx';
import Toxicityprediction from './pages/Toxicityprediction/Toxicityprediction.jsx';

const ProtectedRoute = ({ children }) => {
  const { user } = useAuthStore();
  console.log(user)
  return user ? children : <Navigate to="/login" replace />;
};

function App() {


  return (
    <div className="min-h-screen bg-gray-100">
      {/* <Navbar /> */}
      <Routes>
        <Route path="/" element={<Homepage />} />
        <Route path="/signup" element={<Signup />} />
        <Route path="/login" element={<Login />} />
        <Route path="/verifyotp" element={<VerifyPhone />} />
        <Route path='/dashboard' element={
          <ProtectedRoute>
            <Dashboard />
          </ProtectedRoute>
        } >
          <Route index element={<DashboardHome />} /> {/* Default route for /dashboard */}
          <Route
            path="protein-structure"
            element={
              <ProtectedRoute>
                <ProteinStructureApp />
              </ProtectedRoute>
            }
          />
          <Route
            path="protein-structure-mutation"
            element={
              <ProtectedRoute>
                <ProteinStructureEvolution />
              </ProtectedRoute>
            }
          />
          <Route
            path="cost-estimation"
            element={
              <ProtectedRoute>
                <Costestimation />
              </ProtectedRoute>
            }
          />
          <Route
            path="ai-research-paper-generator"
            element={
              <ProtectedRoute>
                <AIResearchPaperGenerator />
              </ProtectedRoute>
            }
          />
          <Route
            path="ai-driven-target-prediction"
            element={
              <ProtectedRoute>
                <AIDrivenTargetPrediction />
              </ProtectedRoute>
            }
          />
          <Route
            path="ai-naming"
            element={
              <ProtectedRoute>
                <AINamingSuggestion />
              </ProtectedRoute>
            }
          />
          {/* <Route
            path="drug-discovery-recommendation"
            element={
              <ProtectedRoute>
                <DrugDiscoveryRecommendation />
              </ProtectedRoute>
            }
          /> */}
          <Route
            path="live-news"
            element={
              <ProtectedRoute>
                <LiveNews />
              </ProtectedRoute>
            }
          />
          <Route
            path="toxicity-prediction"
            element={
              <ProtectedRoute>
                <Toxicityprediction />
              </ProtectedRoute>
            }
          />
          <Route
            path="getalphafoldstrcture"
            element={
              <ProtectedRoute>
                <GetAlphaFoldStructure />
              </ProtectedRoute>
            }
          />



        </Route>

      </Routes>
      <Toaster />
    </div>
  );
}

export default App;