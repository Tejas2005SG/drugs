// App.jsx
import React, { useEffect, useState } from 'react';
import { Routes, Route, Navigate } from 'react-router-dom';
import { Toaster } from 'react-hot-toast';
import ProteinStructureApp from './pages/Proteinstructureapp/ProteinStructureApp.jsx';
import Homepage from './pages/Homepage/Homepage.jsx';
import Login from './pages/Login/Loginpage.jsx';
import Signup from './pages/SIgnup/Signuppage.jsx';
import { useAuthStore } from './Store/auth.store.js';
import Navbar from './components/Navbar.jsx';
import Dashboard from './components/Dashboard.jsx';
import DashboardHome from './components/Dashboardhome.jsx';

const ProtectedRoute = ({ children }) => {
  const { user } = useAuthStore();
  console.log(user)
  return user ? children : <Navigate to="/login" replace />;
};

function App() {


  return (
    <div className="min-h-screen bg-gray-100">
      <Navbar />
      <Routes>
        <Route path="/" element={<Homepage />} />
        <Route path="/signup" element={<Signup />} />
        <Route path="/login" element={<Login />} />
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

        </Route>

      </Routes>
      <Toaster />
    </div>
  );
}

export default App;