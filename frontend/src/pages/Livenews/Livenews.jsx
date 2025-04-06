import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { useAuthStore } from '../../Store/auth.store.js';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode==="development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

function LiveNews() {
  const [liveArticles, setLiveArticles] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  const token = user ? user.token || user._id || JSON.stringify(user) : null;

  useEffect(() => {
    const fetchLiveNews = async () => {
      if (!token) return;
      try {
        setLoading(true);
        const response = await axiosInstance.get(
          '/news/top-headlines?q=drugs&category=health',
          {
            headers: {
              Authorization: `Bearer ${token}`,
            },
          }
        );
        setLiveArticles(response.data.articles || []);
        setLoading(false);
      } catch (err) {
        setError(err.response?.data?.message || 'Failed to fetch live news');
        setLoading(false);
      }
    };

    const initializeLiveNews = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError('Please log in to view live news');
        setLoading(false);
        return;
      }
      fetchLiveNews();
      const interval = setInterval(fetchLiveNews, 300000);
      return () => clearInterval(interval);
    };

    initializeLiveNews();
  }, [checkAuth, token]);

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-50">
        <div className="text-center space-y-4">
          <div className="inline-block animate-spin rounded-full h-12 w-12 border-4 border-t-transparent border-blue-600"></div>
          <p className="text-gray-600 text-lg">Verifying authentication...</p>
        </div>
      </div>
    );
  }

  if (!token) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-50">
        <div className="text-center p-8 bg-white rounded-2xl shadow-lg max-w-md space-y-6">
          <div className="space-y-4">
            <div className="text-blue-600 mx-auto animate-pulse">
              <svg className="w-20 h-20" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="1.5" d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z"/>
              </svg>
            </div>
            <h2 className="text-3xl font-bold text-gray-800">Access Restricted</h2>
            <p className="text-gray-600 text-lg">Sign in to view real-time health updates</p>
          </div>
          <button
            className="w-full px-6 py-3 bg-blue-600 text-white rounded-xl 
                      hover:bg-blue-700 transition-all duration-300 font-semibold 
                      shadow-md hover:shadow-lg"
            onClick={() => window.location.href = '/login'}
          >
            Continue to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-6xl mx-auto my-8 bg-white rounded-2xl shadow-sm overflow-hidden 
                   transition-all duration-300 hover:shadow-md">
      <header className="p-6 bg-gradient-to-r from-blue-600 to-purple-600 text-white">
        <div className="flex items-center gap-4">
          <div className="p-3 bg-white/10 rounded-xl backdrop-blur-sm">
            <svg className="w-8 h-8" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" 
                    d="M19 14l-7 7m0 0l-7-7m7 7V3"/>
            </svg>
          </div>
          <div>
            <h1 className="text-2xl font-bold tracking-tight">Drugs News Live</h1>
            <p className="mt-1 text-sm opacity-90 font-medium">Real-time medical breakthroughs and updates</p>
          </div>
        </div>
      </header>

      <main className="relative">
        {loading ? (
          <div className="grid gap-6 p-6">
            {[...Array(5)].map((_, i) => (
              <div key={i} className="animate-pulse space-y-4">
                <div className="flex gap-4">
                  <div className="h-6 w-24 bg-gray-200 rounded-full"></div>
                  <div className="flex-1 space-y-2">
                    <div className="h-5 bg-gray-300 rounded w-3/4"></div>
                    <div className="h-4 bg-gray-200 rounded w-full"></div>
                    <div className="h-4 bg-gray-200 rounded w-1/2"></div>
                  </div>
                </div>
                <div className="border-b border-gray-100"></div>
              </div>
            ))}
          </div>
        ) : error ? (
          <div className="p-6 m-6 bg-red-50 rounded-xl flex items-center gap-4">
            <svg className="w-8 h-8 text-red-600 flex-shrink-0" 
                 fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" 
                    d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"/>
            </svg>
            <div>
              <h3 className="font-medium text-red-700">Error loading news</h3>
              <p className="text-sm text-red-600 mt-1">{error}</p>
            </div>
          </div>
        ) : (
          <div className="divide-y divide-gray-100">
            {liveArticles.slice(0, 10).map((article, index) => (
              <article 
                key={index}
                className="group p-6 hover:bg-gray-50 transition-colors duration-200"
              >
                <div className="flex gap-4">
                  <div className="w-28 flex-shrink-0">
                    <time className="text-sm font-medium text-gray-500">
                      {new Date(article.publishedAt).toLocaleTimeString([], {
                        hour: 'numeric',
                        minute: '2-digit',
                        hour12: true,
                      })}
                    </time>
                  </div>
                  <div className="flex-1 space-y-3">
                    <h2 className="text-lg font-semibold text-gray-800">
                      <a 
                        href={article.url} 
                        target="_blank" 
                        rel="noopener noreferrer"
                        className="hover:text-blue-600 transition-colors"
                      >
                        {article.title}
                      </a>
                    </h2>
                    {article.urlToImage && (
                      <img 
                        src={article.urlToImage} 
                        alt={article.title}
                        className="w-full h-48 object-cover rounded-lg mb-3 shadow-sm"
                        loading="lazy"
                        onError={(e) => e.target.style.display = 'none'}
                      />
                    )}
                    <p className="text-gray-600 leading-relaxed">
                      {article.description}
                    </p>
                    <div className="flex items-center justify-between text-sm">
                      <span className="font-medium text-blue-600">
                        {article.source?.name}
                      </span>
                      <a
                        href={article.url}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="flex items-center gap-2 text-gray-500 hover:text-gray-700 transition-colors"
                      >
                        Read Full Article
                        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M17 8l4 4m0 0l-4 4m4-4H3"/>
                        </svg>
                      </a>
                    </div>
                  </div>
                </div>
              </article>
            ))}
          </div>
        )}

        <footer className="p-4 text-center text-sm text-gray-500 bg-gray-50">
          <div className="flex items-center justify-center gap-2">
            <span className="w-2 h-2 bg-blue-500 rounded-full animate-pulse"></span>
            <span className="w-2 h-2 bg-blue-500 rounded-full animate-pulse delay-100"></span>
            <span className="w-2 h-2 bg-blue-500 rounded-full animate-pulse delay-200"></span>
            <span className="ml-2">As this news API offers limited access in the free version, only one or two news articles may be available, and they may not include timestamps(Will change this api soon). </span>
          </div>
        </footer>
      </main>
    </div>
  );
}

export default LiveNews;