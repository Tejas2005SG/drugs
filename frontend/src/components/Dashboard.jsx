import React, { useEffect, useRef, useState } from 'react';
import { useAuthStore } from '../Store/auth.store.js';
import { Outlet, Link, useNavigate } from 'react-router-dom';
import { UserPlus, LogIn, LogOut, Menu, ChevronRight, Activity, Settings, Home, Layers, Dna, DollarSign, FileText, Target, Pill, Newspaper, MessageSquare, X } from 'lucide-react';
import { toast, Toaster } from 'react-hot-toast';

function DashboardPage() {
  const { user, logout } = useAuthStore();
  const navigate = useNavigate();
  const toastShown = useRef(false);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);

  // useEffect(() => {
  //   if (!toastShown.current) {
  //     toast(
  //       "Our system is under development, integrating Gemini, MolMIM, and Grok models to train the AI model with the right mentorship. Stay tuned!",
  //       {
  //         position: "top-right",
  //         duration: 8000,
  //         style: {
  //           background: "#fefcbf",
  //           color: "#92400e",
  //           border: "1px solid #f59e0b",
  //         },
  //         icon: "⚠️",
  //       }
  //     );
  //     toastShown.current = true;
  //   }
  // }, []);

  useEffect(() => {
    console.log('DashboardPage - User:', user);
  }, [user]);

  // Close sidebar when route changes
  useEffect(() => {
    setIsSidebarOpen(false);
  }, [location.pathname]);

  const navElements = [
    { 
      name: 'Dashboard Home', 
      icon: <Home size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'Alphafold Structure Generation', 
      icon: <MessageSquare size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/getalphafoldstrcture'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'Protein Structure Generation', 
      icon: <Dna size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/protein-structure'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'New Drug Discovery', 
      icon: <Layers size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/protein-structure-mutation'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'Cost Estimation', 
      icon: <DollarSign size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/cost-estimation'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'AI Research Paper Generator', 
      icon: <FileText size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/ai-research-paper-generator'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'AI Driven Target Prediction', 
      icon: <Target size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/ai-driven-target-prediction'), 
      roles: ['admin', 'citizen', 'guest']
    },
    { 
      name: 'AI Naming suggestions', 
      icon: <Target size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/ai-naming'), 
      roles: ['admin', 'citizen', 'guest']
    },
   
    { 
      name: 'Live News', 
      icon: <Newspaper size={20} className="mr-3" />,
      navigation: () => navigate('/dashboard/live-news'), 
      roles: ['admin', 'citizen', 'guest']
    },
   
  ];

  const listNav = navElements
    .filter(navElement => {
      const userRole = user?.role || 'guest';
      return navElement.roles.includes(userRole);
    })
    .map((navElement, index) => (
      <li
        key={index}
        onClick={() => {
          navElement.navigation();
          setIsSidebarOpen(false);
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
    <div className="flex flex-col md:flex-row h-screen overflow-hidden bg-gray-100">
      {/* Mobile Header */}
      <div className="fixed top-0 left-0 right-0 h-[48px] bg-white shadow-md flex items-center justify-between px-4 md:hidden z-50">
        <button
          onClick={() => setIsSidebarOpen(!isSidebarOpen)}
          className="p-2 hover:bg-gray-100 rounded-lg"
        >
          {isSidebarOpen ? <X size={24} /> : <Menu size={24} />}
        </button>
        <span className="font-semibold text-gray-800">Dashboard</span>
        <div className="w-8" /> {/* Spacer for balance */}
      </div>

      {/* Overlay for mobile */}
      {isSidebarOpen && (
        <div
          className="fixed inset-0 bg-black bg-opacity-50 z-40 md:hidden"
          onClick={() => setIsSidebarOpen(false)}
        />
      )}

      {/* Sidebar */}
      <div className={`
        fixed md:static w-[280px] md:w-64 bg-white shadow-lg z-40 h-[calc(100vh-48px)] md:h-screen
        transition-transform duration-300 ease-in-out
        ${isSidebarOpen ? 'translate-x-0' : '-translate-x-full md:translate-x-0'}
        mt-[48px] md:mt-0 overflow-hidden
      `}>
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
      <div className="flex-1 flex flex-col h-screen md:h-screen">
        <div className="flex-1 overflow-y-auto pt-[48px] md:pt-0">
          <div className="p-4 md:p-6 h-full">
            <Outlet />
          </div>
        </div>
      </div>

      {/* Toaster */}
      <Toaster />
    </div>
  );
}

export default DashboardPage;