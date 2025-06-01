"use client"

import { useState, useEffect } from "react"
import { Link, useNavigate } from "react-router-dom"
import { useAuthStore } from "../Store/auth.store.js"
import { LogOut, UserPlus, LogIn, Menu, X, Atom, User, ChevronRight, Bell, Settings } from "lucide-react"
import { motion, AnimatePresence, useScroll, useMotionValueEvent } from "framer-motion"

function Navbar() {
  const { user, logout } = useAuthStore()
  const navigate = useNavigate()
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)
  const [isScrolled, setIsScrolled] = useState(false)
  const [activeSection, setActiveSection] = useState("")
  const { scrollY } = useScroll()

  // Handle scroll effects
  useMotionValueEvent(scrollY, "change", (latest) => {
    setIsScrolled(latest > 50)
  })

  // Track active section for navigation highlighting
  useEffect(() => {
    const handleScroll = () => {
      const sections = ["features", "testimonials", "contact"]
      const scrollPosition = window.scrollY + 100

      for (const section of sections) {
        const element = document.getElementById(section)
        if (element) {
          const offsetTop = element.offsetTop
          const offsetHeight = element.offsetHeight

          if (scrollPosition >= offsetTop && scrollPosition < offsetTop + offsetHeight) {
            setActiveSection(section)
            break
          }
        }
      }
    }

    window.addEventListener("scroll", handleScroll)
    return () => window.removeEventListener("scroll", handleScroll)
  }, [])

  const handleLogout = async () => {
    await logout()
    navigate("/")
    setMobileMenuOpen(false)
  }

  const navItems = [
    { href: "#features", label: "Features" },
    { href: "#testimonials", label: "Testimonials" },
    { href: "#contact", label: "Contact" },
  ]

  return (
    <>
      <style jsx>{`
        .navbar-blur {
          backdrop-filter: blur(20px);
          -webkit-backdrop-filter: blur(20px);
        }
        
        .gradient-border {
          background: linear-gradient(90deg, #3B82F6, #8B5CF6, #3B82F6);
          background-size: 200% 200%;
          animation: gradient-shift 3s ease infinite;
        }
        
        @keyframes gradient-shift {
          0%, 100% { background-position: 0% 50%; }
          50% { background-position: 100% 50%; }
        }
        
        .nav-link-active {
          background: linear-gradient(90deg, #3B82F6, #8B5CF6);
          -webkit-background-clip: text;
          background-clip: text;
          -webkit-text-fill-color: transparent;
        }
      `}</style>

      <motion.nav
        className={`fixed w-full z-50 transition-all duration-300 ${
          isScrolled ? "bg-white/80 navbar-blur shadow-lg border-b border-gray-200/50" : "bg-white/95 shadow-sm"
        }`}
        initial={{ y: -100 }}
        animate={{ y: 0 }}
        transition={{ duration: 0.6, ease: "easeOut" }}
      >
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between h-16">
            {/* Logo Section */}
            <motion.div
              className="flex items-center"
              whileHover={{ scale: 1.02 }}
              transition={{ type: "spring", stiffness: 400, damping: 10 }}
            >
              <Link to="/" className="flex items-center gap-3 group">
                <motion.div
                  className="relative flex items-center justify-center w-10 h-10 rounded-xl bg-gradient-to-br from-blue-500 to-purple-600 text-white shadow-lg"
                  whileHover={{ rotate: 360, scale: 1.1 }}
                  transition={{ duration: 0.6, ease: "easeInOut" }}
                >
                  <Atom size={22} />
                  <motion.div
                    className="absolute inset-0 rounded-xl bg-gradient-to-br from-blue-400 to-purple-500 opacity-0 group-hover:opacity-100"
                    initial={{ scale: 0.8, opacity: 0 }}
                    whileHover={{ scale: 1.2, opacity: 0.3 }}
                    transition={{ duration: 0.3 }}
                  />
                </motion.div>
                <motion.span
                  className="text-xl md:text-2xl font-bold bg-gradient-to-r from-blue-600 to-purple-700 bg-clip-text text-transparent"
                  initial={{ opacity: 0, x: -20 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ delay: 0.2, duration: 0.5 }}
                >
                  Drug Discovery AI
                </motion.span>
              </Link>
            </motion.div>

            {/* Desktop Navigation */}
            <div className="hidden md:flex items-center space-x-8">
              {navItems.map((item, index) => (
                <motion.a
                  key={item.href}
                  href={item.href}
                  className={`relative font-medium transition-all duration-300 group ${
                    activeSection === item.href.slice(1) ? "nav-link-active" : "text-gray-700 hover:text-blue-600"
                  }`}
                  initial={{ opacity: 0, y: -20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: 0.1 * index, duration: 0.5 }}
                  whileHover={{ y: -2 }}
                >
                  {item.label}
                  <motion.span
                    className="absolute bottom-0 left-0 h-0.5 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full"
                    initial={{ width: 0 }}
                    whileHover={{ width: "100%" }}
                    transition={{ duration: 0.3, ease: "easeInOut" }}
                  />
                  {activeSection === item.href.slice(1) && (
                    <motion.span
                      className="absolute bottom-0 left-0 h-0.5 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full"
                      initial={{ width: 0 }}
                      animate={{ width: "100%" }}
                      transition={{ duration: 0.3 }}
                    />
                  )}
                </motion.a>
              ))}

              {/* User Section */}
              {user ? (
                <motion.div
                  className="flex items-center ml-6 space-x-4"
                  initial={{ opacity: 0, scale: 0.8 }}
                  animate={{ opacity: 1, scale: 1 }}
                  transition={{ delay: 0.3, duration: 0.5 }}
                >
                  {/* Notifications */}
                  <motion.button
                    className="relative p-2 text-gray-600 hover:text-blue-600 transition-colors duration-200"
                    whileHover={{ scale: 1.1 }}
                    whileTap={{ scale: 0.95 }}
                  >
                    <Bell size={20} />
                    <motion.span
                      className="absolute -top-1 -right-1 w-3 h-3 bg-red-500 rounded-full"
                      animate={{ scale: [1, 1.2, 1] }}
                      transition={{ duration: 2, repeat: Number.POSITIVE_INFINITY }}
                    />
                  </motion.button>

                  {/* User Profile */}
                  <div className="flex items-center space-x-3">
                    <motion.div
                      className="w-8 h-8 rounded-full bg-gradient-to-br from-blue-500 to-purple-600 flex items-center justify-center text-white shadow-md"
                      whileHover={{ scale: 1.1, rotate: 5 }}
                    >
                      <User size={16} />
                    </motion.div>
                    <div className="hidden lg:block">
                      <div className="text-xs text-gray-500">Welcome back,</div>
                      <div className="font-semibold text-gray-800">{user.firstName || "User"}</div>
                    </div>
                  </div>

                  {/* Settings */}
                  <motion.button
                    className="p-2 text-gray-600 hover:text-blue-600 transition-colors duration-200"
                    whileHover={{ scale: 1.1, rotate: 90 }}
                    whileTap={{ scale: 0.95 }}
                  >
                    <Settings size={20} />
                  </motion.button>

                  {/* Logout Button */}
                  <motion.button
                    className="group relative px-4 py-2 overflow-hidden rounded-lg bg-gradient-to-r from-gray-100 to-gray-200 text-gray-700 transition-all duration-300 hover:shadow-md hover:from-red-50 hover:to-red-100 hover:text-red-600"
                    onClick={handleLogout}
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.95 }}
                  >
                    <span className="relative font-medium flex items-center">
                      Log Out
                      <motion.div
                        className="ml-2"
                        whileHover={{ x: 3 }}
                        transition={{ type: "spring", stiffness: 400, damping: 10 }}
                      >
                        <LogOut size={16} />
                      </motion.div>
                    </span>
                  </motion.button>
                </motion.div>
              ) : (
                <motion.div
                  className="flex items-center space-x-4 ml-6"
                  initial={{ opacity: 0, x: 20 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ delay: 0.3, duration: 0.5 }}
                >
                  <Link to="/login">
                    <motion.button
                      className="group relative px-6 py-2 overflow-hidden rounded-lg border-2 border-gray-300 text-gray-700 transition-all duration-300 hover:border-blue-500 hover:text-blue-600"
                      whileHover={{ scale: 1.05 }}
                      whileTap={{ scale: 0.95 }}
                    >
                      <span className="relative font-medium flex items-center">
                        Login
                        <motion.div
                          className="ml-2"
                          whileHover={{ x: 3 }}
                          transition={{ type: "spring", stiffness: 400, damping: 10 }}
                        >
                          <LogIn size={16} />
                        </motion.div>
                      </span>
                    </motion.button>
                  </Link>

                  <Link to="/signup">
                    <motion.button
                      className="group relative px-6 py-2 overflow-hidden rounded-lg bg-gradient-to-r from-blue-600 to-purple-600 text-white transition-all duration-300 hover:shadow-lg hover:from-blue-700 hover:to-purple-700"
                      whileHover={{ scale: 1.05, boxShadow: "0 10px 25px -5px rgba(59, 130, 246, 0.4)" }}
                      whileTap={{ scale: 0.95 }}
                    >
                      <motion.div
                        className="absolute inset-0 bg-gradient-to-r from-blue-400 to-purple-500 opacity-0 group-hover:opacity-100 transition-opacity duration-300"
                        initial={{ scale: 0.8 }}
                        whileHover={{ scale: 1 }}
                      />
                      <span className="relative font-medium flex items-center">
                        Sign Up
                        <motion.div
                          className="ml-2"
                          whileHover={{ x: 3 }}
                          transition={{ type: "spring", stiffness: 400, damping: 10 }}
                        >
                          <ChevronRight size={16} />
                        </motion.div>
                      </span>
                    </motion.button>
                  </Link>
                </motion.div>
              )}
            </div>

            {/* Mobile menu button */}
            <div className="md:hidden flex items-center">
              <motion.button
                onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
                className="inline-flex items-center justify-center p-2 rounded-lg text-gray-600 hover:text-blue-600 hover:bg-gray-100 focus:outline-none transition-all duration-200"
                whileHover={{ scale: 1.1 }}
                whileTap={{ scale: 0.95 }}
                aria-expanded="false"
              >
                <span className="sr-only">Open main menu</span>
                <AnimatePresence mode="wait">
                  {mobileMenuOpen ? (
                    <motion.div
                      key="close"
                      initial={{ rotate: -90, opacity: 0 }}
                      animate={{ rotate: 0, opacity: 1 }}
                      exit={{ rotate: 90, opacity: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <X size={24} />
                    </motion.div>
                  ) : (
                    <motion.div
                      key="menu"
                      initial={{ rotate: 90, opacity: 0 }}
                      animate={{ rotate: 0, opacity: 1 }}
                      exit={{ rotate: -90, opacity: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <Menu size={24} />
                    </motion.div>
                  )}
                </AnimatePresence>
              </motion.button>
            </div>
          </div>
        </div>

        {/* Mobile Menu */}
        <AnimatePresence>
          {mobileMenuOpen && (
            <motion.div
              className="md:hidden fixed inset-x-0 top-16 z-50 bg-white/95 navbar-blur shadow-xl border-b border-gray-200"
              initial={{ opacity: 0, y: -20, scale: 0.95 }}
              animate={{ opacity: 1, y: 0, scale: 1 }}
              exit={{ opacity: 0, y: -20, scale: 0.95 }}
              transition={{ duration: 0.2, ease: "easeOut" }}
            >
              <div className="px-4 pt-4 pb-6 space-y-3">
                {/* Navigation Links */}
                {navItems.map((item, index) => (
                  <motion.a
                    key={item.href}
                    href={item.href}
                    className="block px-4 py-3 rounded-lg text-base font-medium text-gray-700 hover:text-blue-600 hover:bg-blue-50 transition-all duration-200"
                    onClick={() => setMobileMenuOpen(false)}
                    initial={{ opacity: 0, x: -20 }}
                    animate={{ opacity: 1, x: 0 }}
                    transition={{ delay: index * 0.1, duration: 0.3 }}
                    whileHover={{ x: 5 }}
                  >
                    {item.label}
                  </motion.a>
                ))}

                {/* User Section */}
                {user ? (
                  <motion.div
                    className="pt-4 border-t border-gray-200"
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.3, duration: 0.3 }}
                  >
                    <div className="px-4 py-3 flex items-center space-x-3">
                      <motion.div
                        className="w-10 h-10 rounded-full bg-gradient-to-br from-blue-500 to-purple-600 flex items-center justify-center text-white"
                        whileHover={{ scale: 1.1, rotate: 5 }}
                      >
                        <User size={18} />
                      </motion.div>
                      <div>
                        <div className="text-xs text-gray-500">Logged in as</div>
                        <div className="font-semibold text-gray-800">{user.firstName || "User"}</div>
                      </div>
                    </div>

                    <motion.button
                      className="w-full mt-3 flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium text-red-600 bg-red-50 hover:bg-red-100 transition-colors duration-200"
                      onClick={handleLogout}
                      whileHover={{ scale: 1.02 }}
                      whileTap={{ scale: 0.98 }}
                    >
                      <span>Log Out</span>
                      <LogOut size={18} />
                    </motion.button>
                  </motion.div>
                ) : (
                  <motion.div
                    className="space-y-3 pt-4 border-t border-gray-200"
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.3, duration: 0.3 }}
                  >
                    <Link to="/login" onClick={() => setMobileMenuOpen(false)}>
                      <motion.button
                        className="w-full flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium text-gray-700 border-2 border-gray-300 hover:border-blue-500 hover:text-blue-600 hover:bg-blue-50 transition-all duration-200"
                        whileHover={{ scale: 1.02 }}
                        whileTap={{ scale: 0.98 }}
                      >
                        <span>Login</span>
                        <LogIn size={18} />
                      </motion.button>
                    </Link>

                    <Link to="/signup" onClick={() => setMobileMenuOpen(false)}>
                      <motion.button
                        className="w-full flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium text-white bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 transition-all duration-200 shadow-md"
                        whileHover={{ scale: 1.02, boxShadow: "0 8px 20px -5px rgba(59, 130, 246, 0.4)" }}
                        whileTap={{ scale: 0.98 }}
                      >
                        <span>Sign Up</span>
                        <UserPlus size={18} />
                      </motion.button>
                    </Link>
                  </motion.div>
                )}
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </motion.nav>
    </>
  )
}

export default Navbar
