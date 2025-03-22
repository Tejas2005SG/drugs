import React from 'react';
import { Link, useNavigate } from 'react-router-dom'; // Add useNavigate for logout redirect
import { useAuthStore } from '../Store/auth.store.js';
import { LogOut, UserPlus, LogIn } from 'lucide-react'; // Import icons from lucide-react

function Navbar() {
  const { user, logout } = useAuthStore();
  const navigate = useNavigate(); // For redirecting after logout

  const handleLogout = async () => {
    await logout(); // Call the logout function from the store
    navigate('/'); // Redirect to homepage after logout
  };

  return (
    <div>
      <nav className="bg-white shadow-lg fixed w-full z-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between h-16">
            <div className="flex items-center">
              <Link to="/" className="text-2xl font-bold text-blue-600">Drug Discovery Assistant</Link>
            </div>
            <div className="flex items-center space-x-4">
              <a href="#features" className="text-gray-600 hover:text-blue-600 font-medium">
                Features
              </a>
              <a href="#about" className="text-gray-600 hover:text-blue-600 font-medium">
                About
              </a>
              <a href="#contact" className="text-gray-600 hover:text-blue-600 font-medium">
                Contact
              </a>
              {user ? (
                <button
                  className="bg-gray-700 hover:bg-gray-600 text-white py-2 px-4 rounded-md flex items-center transition duration-300 ease-in-out"
                  onClick={handleLogout}
                >
                  <LogOut size={18} />
                  <span className="hidden sm:inline ml-2">Log Out</span>
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
          </div>
        </div>
      </nav>
    </div>
  );
}

export default Navbar;