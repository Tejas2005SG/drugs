import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/vite'; // Correct import for Tailwind CSS plugin

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    react(),
    tailwindcss(), // Ensure Tailwind CSS is correctly integrated
  ],
  server: {
    port: 5173, // Default Vite port, adjust if needed
    proxy: {
      '/api': {
        target: 'http://localhost:5000', // Backend server URL
        changeOrigin: true, // Needed for localhost to work properly
        rewrite: (path) => path.replace(/^\/api/, '/api'), // Optional: preserves '/api' prefix
      },
    },
  },
});