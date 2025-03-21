// App.js
import React from 'react';
import { Routes, Route } from 'react-router-dom';
import { Toaster } from 'react-hot-toast';

import ProteinStructureApp from './pages/Proteinstructureapp/ProteinStructureApp.jsx';
import Homepage from './pages/Homepage/Homepage.jsx';
import Login from './pages/Login/Loginpage.jsx';
import Signup from './pages/SIgnup/Signuppage.jsx';

function App() {
  return (
    <div className="min-h-screen bg-gray-100">
      <Routes>
        <Route path="/" element={<Homepage />} /> 
        <Route path="/signup" element={<Signup />} />
        <Route path="/login" element={<Login/>} />
        <Route path="/protein-structure" element={<ProteinStructureApp />} />
      </Routes>
      <Toaster />
    </div>
  );
}

export default App;