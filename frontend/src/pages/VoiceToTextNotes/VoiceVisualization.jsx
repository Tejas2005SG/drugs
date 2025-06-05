import React, { useRef, useEffect } from 'react';

const VoiceVisualization = ({ audioLevel = 0, theme }) => {
  const canvasRef = useRef(null);
  
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    const isDarkMode = theme === 'dark';
    
    // Set canvas dimensions
    canvas.width = canvas.offsetWidth;
    canvas.height = canvas.offsetHeight;
    
    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    // Set colors
    const gradient = ctx.createLinearGradient(0, 0, canvas.width, 0);
    gradient.addColorStop(0, isDarkMode ? '#a78bfa' : '#8b5cf6');
    gradient.addColorStop(1, isDarkMode ? '#60a5fa' : '#3b82f6');
    
    // Number of bars
    const barCount = 40;
    const barWidth = canvas.width / barCount - 2;
    
    // Draw visualization
    for (let i = 0; i < barCount; i++) {
      // Calculate bar height based on audio level and position
      let intensity = Math.sin((i / barCount) * Math.PI) * audioLevel;
      // Add some randomness
      intensity *= 0.7 + Math.random() * 0.3;
      
      // Ensure minimum height
      const barHeight = Math.max(4, canvas.height * intensity);
      
      // Position in the middle vertically
      const y = (canvas.height - barHeight) / 2;
      
      // Draw bar
      ctx.fillStyle = gradient;
      ctx.fillRect(i * (barWidth + 2), y, barWidth, barHeight);
    }
  }, [audioLevel, theme]);
  
  return (
    <div className={`w-full h-12 rounded-lg overflow-hidden ${
      theme === 'light' ? 'bg-gray-100' : 'bg-gray-700/50'
    }`}>
      <canvas ref={canvasRef} className="w-full h-full" />
    </div>
  );
};

export default VoiceVisualization;