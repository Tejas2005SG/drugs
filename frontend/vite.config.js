import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig(({ mode }) => {
  // Load env variables based on mode (development/production)
  const env = loadEnv(mode, process.cwd(), '');

  // Determine API base URL based on environment
  const apiBaseUrl = mode === 'development' 
    ? 'http://localhost:5000' 
    : env.VITE_API_BASE_URL || 'https://your-render-backend-url.onrender.com';

  return {
    plugins: [react()],
    server: {
      proxy: {
        '/api': {
          target: apiBaseUrl,
          changeOrigin: true,
          secure: false,
          rewrite: (path) => path.replace(/^\/api/, '')
        },
        '/socket.io': {
          target: apiBaseUrl,
          ws: true,
          changeOrigin: true
        }
      }
    },
    // Ensure environment variables are available in client code
    define: {
      'process.env': {
        VITE_API_BASE_URL: JSON.stringify(apiBaseUrl)
      }
    }
  };
});