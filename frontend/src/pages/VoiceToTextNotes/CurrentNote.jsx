import React, { useState } from 'react';
import { MessageSquare, AlertTriangle, CheckCircle, Info } from 'lucide-react';

const CurrentNote = ({ currentNote, setCurrentNote, sentiment, summary, isEditing, theme }) => {
  const [showMarkdownHelp, setShowMarkdownHelp] = useState(false);
  
  const renderSentimentIcon = () => {
    if (!sentiment) return null;
    
    switch(sentiment) {
      case 'positive':
        return <CheckCircle size={16} className="text-green-500" />;
      case 'negative':
        return <AlertTriangle size={16} className="text-red-500" />;
      default:
        return <Info size={16} className="text-yellow-500" />;
    }
  };
  
  const handleNoteChange = (e) => {
    setCurrentNote(e.target.value);
  };
  
  return (
    <div className={`p-4 md:p-6 border-b ${
      theme === 'light' ? 'border-gray-200' : 'border-gray-700'
    }`}>
      <div className="flex justify-between items-center mb-4">
        <h2 className={`flex items-center gap-2 text-lg font-semibold ${
          theme === 'light' ? 'text-gray-700' : 'text-gray-200'
        }`}>
          <MessageSquare size={18} />
          {isEditing ? 'Edit Note' : 'Current Note'}
        </h2>
        {/* <button
          className={`text-sm ${
            theme === 'light' 
              ? 'text-indigo-600 hover:bg-indigo-50' 
              : 'text-indigo-400 hover:bg-indigo-900/30'
          } rounded px-2 py-1 transition-colors`}
          onClick={() => setShowMarkdownHelp(!showMarkdownHelp)}
        >
          Formatting Help
        </button> */}
      </div>
      
      {/* {showMarkdownHelp && (
        <div className={`mb-4 p-3 rounded-md text-sm ${
          theme === 'light' ? 'bg-gray-50' : 'bg-gray-700/50'
        }`}>
          <h3 className="font-medium mb-1">Markdown Formatting</h3>
          <ul className="pl-5 list-disc">
            <li><code className={`px-1 py-0.5 rounded ${
              theme === 'light' ? 'bg-gray-200' : 'bg-gray-600'
            }`}># Heading</code> - Creates a heading</li>
            <li><code className={`px-1 py-0.5 rounded ${
              theme === 'light' ? 'bg-gray-200' : 'bg-gray-600'
            }`}>**Bold**</code> - Makes text bold</li>
            <li><code className={`px-1 py-0.5 rounded ${
              theme === 'light' ? 'bg-gray-200' : 'bg-gray-600'
            }`}>*Italic*</code> - Makes text italic</li>
            <li><code className={`px-1 py-0.5 rounded ${
              theme === 'light' ? 'bg-gray-200' : 'bg-gray-600'
            }`}>- Item</code> - Creates a list item</li>
          </ul>
        </div>
      )} */}
      
      <div className="mb-4">
        <textarea
          value={currentNote}
          onChange={handleNoteChange}
          placeholder="Start speaking or type your note here..."
          className={`w-full min-h-[150px] p-4 border rounded-lg resize-y focus:ring-2 focus:outline-none ${
            theme === 'light'
              ? 'bg-white border-gray-300 text-gray-800 focus:ring-indigo-200 focus:border-indigo-500'
              : 'bg-gray-700 border-gray-600 text-gray-100 focus:ring-indigo-500/30 focus:border-indigo-500'
          }`}
        ></textarea>
      </div>
      
      <div className={`flex flex-wrap gap-4 text-sm ${
        theme === 'light' ? 'text-gray-600' : 'text-gray-400'
      }`}>
        <div className="flex items-center">
          {currentNote.trim() ? currentNote.trim().split(/\s+/).length : 0} words
        </div>
        
        {sentiment && (
          <div className={`flex items-center gap-1 ${
            sentiment === 'positive' 
              ? 'text-green-500' 
              : sentiment === 'negative' 
                ? 'text-red-500' 
                : 'text-yellow-500'
          }`}>
            {renderSentimentIcon()}
            <span>Sentiment: {sentiment}</span>
          </div>
        )}
        
        {summary && (
          <div className="flex items-center">
            <span>Summary: {summary}</span>
          </div>
        )}
      </div>
    </div>
  );
};

export default CurrentNote;