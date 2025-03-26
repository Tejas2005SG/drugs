// components/Dashboard.jsx
import React, { useEffect } from 'react';
import { useAuthStore } from '../Store/auth.store.js';
import { Outlet, Link, useNavigate } from 'react-router-dom';
import { UserPlus, LogIn, LogOut, Menu, ChevronRight, Activity, Settings, Home, Layers, Dna, DollarSign, FileText, Target, Pill, Newspaper, MessageSquare } from 'lucide-react';

function DashboardPage() {
  const { user, logout } = useAuthStore();
  const navigate = useNavigate();

  useEffect(() => {
    console.log('DashboardPage - User:', user); // Debug user object
  }, [user]);


  const navElements = [
      { 
          name: 'Dashboard Home', 
          icon: <Home size={20} className="mr-3" />, // Home fits a main dashboard
          navigation: () => navigate('/dashboard'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'Protein Structure Generation', 
          icon: <Dna size={20} className="mr-3" />, // Dna represents molecular/protein work
          navigation: () => navigate('/dashboard/protein-structure'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'Protein Structure Evolution', 
          icon: <Layers size={20} className="mr-3" />, // Layers suggests evolution or stacking changes
          navigation: () => navigate('/dashboard/protein-structure-mutation'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'Cost Estimation', 
          icon: <DollarSign size={20} className="mr-3" />, // DollarSign for anything cost-related
          navigation: () => navigate('/dashboard/cost-estimation'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'AI Research Paper Generator', 
          icon: <FileText size={20} className="mr-3" />, // FileText for generating papers/documents
          navigation: () => navigate('/dashboard/ai-research-paper-generator'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'AI Driven Target Prediction', 
          icon: <Target size={20} className="mr-3" />, // Target for predicting targets
          navigation: () => navigate('/dashboard/ai-driven-target-prediction'), 
          roles: ['admin', 'citizen', 'guest']
      },
      // { 
      //     name: 'Drug Discovery Recommendation', 
      //     icon: <Pill size={20} className="mr-3" />, // Pill for drug-related features
      //     navigation: () => navigate('/dashboard/drug-discovery-recommendation'), 
      //     roles: ['admin', 'citizen', 'guest']
      // },
      { 
          name: 'Live News', 
          icon: <Newspaper size={20} className="mr-3" />, // Newspaper for news updates
          navigation: () => navigate('/dashboard/live-news'), 
          roles: ['admin', 'citizen', 'guest']
      },
      { 
          name: 'Message Board', 
          icon: <MessageSquare size={20} className="mr-3" />, // MessageSquare for communication
          navigation: () => navigate('/dashboard/message'), 
          roles: ['admin', 'citizen', 'guest']
      },
  ];;

  const listNav = navElements
    .filter(navElement => {
      const userRole = user?.role || 'guest';
      console.log('Filtering navElement:', navElement.name, 'User Role:', userRole, 'Allowed Roles:', navElement.roles);
      return navElement.roles.includes(userRole);
    })
    .map((navElement, index) => (
      <li
        key={index}
        onClick={() => {
          navElement.navigation();
        }}
        className="w-full bg-white hover:bg-blue-50 text-gray-700 hover:text-blue-600 rounded-lg p-3 flex items-center my-1 cursor-pointer transition-all duration-200 shadow-sm"
      >
        {navElement.icon}
        <span className="text-base font-medium">{navElement.name}</span>
        <ChevronRight size={16} className="ml-auto text-gray-400" />
      </li>
    ));

  const handleLogout = () => {
    logout();
    navigate('/');
  };

  return (
    <div className="flex flex-row min-h-screen bg-gray-100">
      {/* Sidebar */}
      <div className="md:block w-64 bg-white shadow-lg z-10 md:relative h-[calc(100vh-48px)] mt-[48px]">
        <div className="flex flex-col h-full">
          {/* User Welcome Section */}
          <div className="p-6 bg-gradient-to-r from-blue-500 to-blue-600 text-white">
            <div className="flex items-center space-x-3">
              <div className="w-12 h-12 rounded-full bg-white flex items-center justify-center text-blue-600 font-bold text-xl">
                {user?.firstName?.charAt(0) || 'U'}
              </div>
              <div>
                <h2 className="text-lg font-semibold">Welcome,</h2>
                <h1 className="text-xl font-bold">{user ? user.firstName : 'Guest'}</h1>
              </div>
            </div>
          </div>

          {/* Navigation */}
          <div className="flex-1 p-4 overflow-y-auto">
            {/* <h3 className="text-xs uppercase text-gray-500 font-semibold tracking-wider mb-4 ml-3">Main Navigation</h3> */}
            <ul className="space-y-2">
              {listNav.length > 0 ? (
                listNav
              ) : (
                <li className="text-gray-500 text-sm">No navigation items available for your role.</li>
              )}
            </ul>
          </div>

          {/* Logout Section */}
          <div className="p-4 border-t border-gray-200">
            {user ? (
              <button
                className="w-full bg-red-300 border border-red-800 hover:bg-red-400 text-gray-800 font-bold py-3 px-4 rounded-lg flex items-center justify-center transition duration-300 ease-in-out"
                onClick={handleLogout}
              >
                <LogOut size={18} className="mr-2" />
                <span>Log Out</span>
              </button>
            ) : (
              <div className="flex flex-col space-y-2">
                <Link
                  to="/signup"
                  className="w-full bg-blue-500 hover:bg-blue-600 text-white py-2 px-4 rounded-lg flex items-center justify-center transition duration-300"
                >
                  <UserPlus className="mr-2" size={18} />
                  Sign Up
                </Link>
                <Link
                  to="/login"
                  className="w-full bg-gray-600 hover:bg-gray-700 text-white py-2 px-4 rounded-lg flex items-center justify-center transition duration-300"
                >
                  <LogIn className="mr-2" size={18} />
                  Login
                </Link>
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 p-6 mt-[48px] h-[calc(100vh-48px)] overflow-y-auto">
        <Outlet />
      </div>
    </div>
  );
}

export default DashboardPage;