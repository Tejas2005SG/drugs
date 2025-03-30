import React, { useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../Store/auth.store.js';
import { LogOut, UserPlus, LogIn, Menu, X, Atom, User, ChevronRight  } from 'lucide-react';

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
    <nav className="bg-zinc-900 shadow-lg fixed w-full z-50">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex justify-between h-16">
          <div className="flex items-center">
            <Link to="/" className="flex items-center gap-2">
              <div className="flex items-center justify-center w-10 h-10 rounded-full bg-gradient-to-r from-blue-500 to-purple-600 text-white">
                <Atom size={22} />
              </div>
              <span className="text-xl md:text-2xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-purple-700">
                Drug Discovery Assistant
              </span>
            </Link>
          </div>

          {/* Desktop Navigation */}
          <div className="hidden md:flex items-center space-x-6">
            <a
              href="#features"
              className="text-white hover:text-blue-600 transition-colors duration-300 font-medium relative group"
            >
              Features
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-blue-500 group-hover:w-full transition-all duration-300"></span>
            </a>
            <a
              href="#contact"
              className="text-white hover:text-blue-600 transition-colors duration-300 font-medium relative group"
            >
              Contact
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-blue-600 group-hover:w-full transition-all duration-300"></span>
            </a>

            {user ? (
              <div className="flex items-center ml-4">
                <div className="mr-4 flex items-center">
                  <div className="w-8 h-8 rounded-full bg-gradient-to-r from-blue-500 to-purple-600 flex items-center justify-center text-white mr-2">
                    <User size={16} />
                  </div>
                  <div>
                    <div className="text-xs text-zinc-400">Welcome,</div>
                    <div className="font-medium text-white">{user.firstName || "User"}</div>
                  </div>
                </div>
                <button
                  className="group relative px-4 py-2 overflow-hidden rounded-md bg-gradient-to-r from-gray-700 to-gray-800 text-white transition-all duration-300 ease-out hover:shadow-md hover:from-gray-800 hover:to-gray-900"
                  onClick={handleLogout}
                >
                  <span className="relative font-medium flex items-center">
                    Log Out
                    <LogOut size={16} className="ml-2 group-hover:translate-x-1 transition-transform duration-200" />
                  </span>
                </button>
              </div>
            ) : (
              <div className="flex items-center space-x-4 ml-4">
                <Link
                  to="/login"
                  className="group relative px-4 py-2 overflow-hidden rounded-md border border-gray-300 text-white transition-all duration-300 ease-out hover:border-gray-400 hover:bg-gray-50 hover:text-black"
                >
                  <span className="relative font-medium flex items-center">
                    Login
                    <LogIn size={16} className="ml-2 group-hover:translate-x-1 transition-transform duration-200" />
                  </span>
                </Link>
                <Link
                  to="/signup"
                  className="group relative px-4 py-2 overflow-hidden rounded-md bg-gradient-to-r from-blue-600 to-purple-600 text-white transition-all duration-300 ease-out hover:shadow-md hover:from-blue-700 hover:to-purple-700"
                >
                  <span className="relative font-medium flex items-center">
                    Sign Up
                    <ChevronRight
                      size={16}
                      className="ml-2 group-hover:translate-x-1 transition-transform duration-200"
                    />
                  </span>
                </Link>
              </div>
            )}
          </div>

          {/* Mobile menu button */}
          {/* Mobile menu button */}
          <div className="md:hidden flex items-center">
            <button
              onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
              className="inline-flex items-center justify-center p-2 rounded-md text-gray-600 hover:text-blue-600 focus:outline-none transition-colors duration-300"
              aria-expanded="false"
            >
              <span className="sr-only">Open main menu</span>
              {mobileMenuOpen ? <X size={24} /> : <Menu size={24} />}
            </button>
          </div>
        </div>
      </div>

      {/* Mobile Menu */}
      {/* <div className={`md:hidden ${mobileMenuOpen ? 'block' : 'hidden'}`}>
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
      </div> */}

      {/* Mobile Menu Dropdown */}
      <div
        className={`md:hidden fixed inset-x-0 top-16 z-50 bg-zinc-900 text-white shadow-lg transition-all duration-300 ease-in-out transform ${mobileMenuOpen
            ? "opacity-100 translate-y-0"
            : "opacity-0 -translate-y-4 pointer-events-none"
          }`}
      >
        <div className="px-2 pt-2 pb-4 space-y-1 border-t border-gray-200">
          <a
            href="#features"
            className="block px-3 py-2 rounded-md text-base font-medium text-white  hover:bg-zinc-800 transition-colors"
            onClick={() => setMobileMenuOpen(false)}
          >
            Features
          </a>
          <a
            href="#contact"
            className="block px-3 py-2 rounded-md text-base font-medium text-white  hover:bg-zinc-800 transition-colors"
            onClick={() => setMobileMenuOpen(false)}
          >
            Contact
          </a>

          {user ? (
            <div className="pt-2 border-t border-gray-200">
              <div className="px-3 py-2 flex items-center">
                <div className="w-8 h-8 rounded-full bg-gradient-to-r from-blue-500 to-purple-600 flex items-center justify-center text-white mr-3">
                  <User className="w-4 h-4" />
                </div>
                <div>
                  <div className="text-xs text-zinc-300">Logged in as</div>
                  <div className="font-medium text-white">{user.firstName || "User"}</div>
                </div>
              </div>
              <button
                className="w-full mt-2 flex items-center justify-between px-3 py-2 rounded-md text-base font-medium text-white bg-zinc-800 hover:from-zinc-700 transition-colors"
                onClick={handleLogout}
              >
                <span>Log Out</span>
                <LogOut className="w-4 h-4" />
              </button>
            </div>
          ) : (
            <div className="space-y-2 pt-2 border-t border-zinc-900">
              <Link
                to="/login"
                className="w-full flex items-center justify-between px-3 py-2 rounded-md text-base font-medium text-gray-700 border border-gray-300 hover:bg-gray-50 transition-colors"
                onClick={() => setMobileMenuOpen(false)}
              >
                <span>Login</span>
                <LogIn className="w-4 h-4" />
              </Link>
              <Link
                to="/signup"
                className="w-full flex items-center justify-between px-3 py-2 rounded-md text-base font-medium text-white bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 transition-colors"
                onClick={() => setMobileMenuOpen(false)}
              >
                <span>Sign Up</span>
                <UserPlus className="w-4 h-4" />
              </Link>
            </div>
          )}
        </div>
      </div>
    </nav>
  );
}

export default Navbar;