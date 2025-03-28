import React, { useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../Store/auth.store.js';
import { LogOut, UserPlus, LogIn, Menu, X } from 'lucide-react';

function Navbar() {
  const { user, logout } = useAuthStore();
  const navigate = useNavigate();
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  const handleLogout = async () => {
    await logout();
    navigate('/');
    setMobileMenuOpen(false); // Close mobile menu on logout
  };

  return (
    <nav className="bg-white shadow-lg fixed w-full z-50">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex justify-between h-16">
          {/* Logo */}
          <div className="flex items-center">
            <Link to="/" className="text-xl md:text-2xl font-bold text-blue-600">
              Drug Discovery Assistant
            </Link>
          </div>

          {/* Desktop Navigation */}
          <div className="hidden md:flex items-center space-x-4">
            <a href="#features" className="text-gray-600 hover:text-blue-600 font-medium px-3 py-2">
              Features
            </a>
            <a href="#contact" className="text-gray-600 hover:text-blue-600 font-medium px-3 py-2">
              Contact
            </a>
            
            {user ? (
              <button
                className="bg-gray-700 hover:bg-gray-600 text-white py-2 px-4 rounded-md flex items-center transition duration-300 ease-in-out"
                onClick={handleLogout}
              >
                <LogOut size={18} />
                <span className="ml-2">Log Out</span>
              </button>
            ) : (
              <>
                <Link
                  to="/signup"
                  className="bg-blue-600 hover:bg-blue-700 text-white py-2 px-4 rounded-md flex items-center transition duration-300 ease-in-out"
                >
                  <UserPlus className="mr-2" size={18} />
                  Sign Up
                </Link>
                <Link
                  to="/login"
                  className="bg-gray-700 hover:bg-gray-600 text-white py-2 px-4 rounded-md flex items-center transition duration-300 ease-in-out"
                >
                  <LogIn className="mr-2" size={18} />
                  Login
                </Link>
              </>
            )}
          </div>

          {/* Mobile menu button */}
          <div className="md:hidden flex items-center">
            <button
              onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
              className="inline-flex items-center justify-center p-2 rounded-md text-gray-600 hover:text-blue-600 focus:outline-none"
              aria-expanded="false"
            >
              <span className="sr-only">Open main menu</span>
              {mobileMenuOpen ? (
                <X size={24} />
              ) : (
                <Menu size={24} />
              )}
            </button>
          </div>
        </div>
      </div>

      {/* Mobile Menu */}
      <div className={`md:hidden ${mobileMenuOpen ? 'block' : 'hidden'}`}>
        <div className="px-2 pt-2 pb-3 space-y-1 sm:px-3 bg-white shadow-lg">
          <a
            href="#features"
            className="block px-3 py-2 rounded-md text-base font-medium text-gray-600 hover:text-blue-600 hover:bg-gray-50"
            onClick={() => setMobileMenuOpen(false)}
          >
            Features
          </a>
          <a
            href="#contact"
            className="block px-3 py-2 rounded-md text-base font-medium text-gray-600 hover:text-blue-600 hover:bg-gray-50"
            onClick={() => setMobileMenuOpen(false)}
          >
            Contact
          </a>
          
          {user ? (
            <button
              className="w-full flex items-center justify-start px-3 py-2 rounded-md text-base font-medium text-white bg-gray-700 hover:bg-gray-600"
              onClick={handleLogout}
            >
              <LogOut className="mr-2" size={18} />
              Log Out
            </button>
          ) : (
            <div className="space-y-2 pt-2">
              <Link
                to="/signup"
                className="w-full flex items-center justify-start px-3 py-2 rounded-md text-base font-medium text-white bg-blue-600 hover:bg-blue-700"
                onClick={() => setMobileMenuOpen(false)}
              >
                <UserPlus className="mr-2" size={18} />
                Sign Up
              </Link>
              <Link
                to="/login"
                className="w-full flex items-center justify-start px-3 py-2 rounded-md text-base font-medium text-white bg-gray-700 hover:bg-gray-600"
                onClick={() => setMobileMenuOpen(false)}
              >
                <LogIn className="mr-2" size={18} />
                Login
              </Link>
            </div>
          )}
        </div>
      </div>
    </nav>
  );
}

export default Navbar;