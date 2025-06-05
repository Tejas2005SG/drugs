import React from 'react';
import { Mic, MicOff, Save, Play, Loader } from 'lucide-react';
import VoiceVisualization from './VoiceVisualization';

const RecordingControls = ({
  isRecording,
  startRecording,
  stopRecording,
  saveNote,
  showContinue,
  error,
  isProcessing,
  hasUnsavedChanges,
  currentNote,
  audioLevel,
  theme,
  isEditing,
}) => {
  return (
    <div className={`p-4 md:p-6 border-b ${theme === 'light' ? 'border-gray-200' : 'border-gray-700'}`}>
      <div className="flex flex-wrap gap-3 justify-center">
        {!showContinue && !isRecording && !isEditing && (
          <button
            onClick={startRecording}
            className={`flex items-center gap-2 px-4 py-2 rounded-full font-medium text-white bg-green-500 hover:bg-green-600 transition-colors shadow-md disabled:opacity-50 disabled:cursor-not-allowed`}
            aria-label="Start recording"
            title="Start recording"
            disabled={!!error}
          >
            <Mic size={18} />
            <span>Start Recording</span>
          </button>
        )}

        {isRecording && (
          <button
            onClick={stopRecording}
            className="flex items-center gap-2 px-4 py-2 rounded-full font-medium text-white bg-red-500 hover:bg-red-600 transition-colors shadow-md disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label="Stop recording"
            title="Stop recording"
          >
            <MicOff size={18} />
            <span>Stop Recording</span>
          </button>
        )}

        {showContinue && currentNote.trim() && !isRecording && !isEditing && (
          <button
            onClick={startRecording}
            className="flex items-center gap-2 px-4 py-2 rounded-full font-medium text-white bg-blue-500 hover:bg-blue-600 transition-colors shadow-md disabled:opacity-50 disabled:cursor-not-allowed"
            aria-label="Continue recording"
            title="Continue recording"
            disabled={!!error}
          >
            <Play size={18} />
            <span>Continue Recording</span>
          </button>
        )}

        <button
          onClick={saveNote}
          className={`flex items-center gap-2 px-4 py-2 rounded-full font-medium text-white bg-purple-500 hover:bg-purple-600 transition-colors shadow-md disabled:opacity-50 disabled:cursor-not-allowed ${
            hasUnsavedChanges ? 'animate-pulse' : ''
          }`}
          aria-label="Save note"
          title="Save note"
          disabled={!currentNote.trim() || !!error}
        >
          {isProcessing ? <Loader className="animate-spin" size={18} /> : <Save size={18} />}
          <span>Save Note</span>
        </button>
      </div>

      {isRecording && (
        <div className="mt-4 flex flex-col items-center">
          <VoiceVisualization audioLevel={audioLevel} theme={theme} />
          <div className="mt-2 flex items-center text-sm gap-2 text-red-500">
            <span className="h-2 w-2 rounded-full bg-red-500 animate-pulse"></span>
            <span>Recording...</span>
          </div>
        </div>
      )}

      {error && (
        <div className={`mt-4 p-3 rounded-md text-sm ${
          theme === 'light' ? 'bg-red-100 text-red-700' : 'bg-red-900/30 text-red-400'
        }`}>
          {error}
        </div>
      )}

      <div className="mt-4">
        <details className="text-sm">
          <summary className={`cursor-pointer ${
            theme === 'light' ? 'text-indigo-600 hover:text-indigo-800' : 'text-indigo-400 hover:text-indigo-300'
          }`}>
            Voice Commands
          </summary>
          <ul className={`mt-2 pl-5 list-disc ${
            theme === 'light' ? 'text-gray-700' : 'text-gray-300'
          }`}>
            <li><strong>"Save note"</strong> - Saves current note</li>
            <li><strong>"Cancel"</strong> - Cancels editing</li>
            <li><strong>"Clear"</strong> - Clears current text</li>
            <li><strong>"Stop recording"</strong> - Stops the recording</li>
          </ul>
        </details>
      </div>
    </div>
  );
};

export default RecordingControls;