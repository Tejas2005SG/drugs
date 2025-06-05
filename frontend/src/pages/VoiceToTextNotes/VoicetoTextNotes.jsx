import React, { useState, useEffect } from 'react';
import { Mic, MicOff, AlertCircle, Sun, Moon } from 'lucide-react';
import RecordingControls from './RecordingControls';
import CurrentNote from './CurrentNote';
import NotesList from './NotesList';
import { useVoiceRecognition } from './hooks/useVoiceRecognition';
import { useNotes } from './hooks/useNotes';
import { useSpeechCommands } from './hooks/useSpeechCommands';

function VoiceToTextNotes() {
  const [theme, setTheme] = useState(() => {
    const savedTheme = localStorage.getItem('voice_notes_theme');
    return savedTheme || 'light';
  });
  const [showAlert, setShowAlert] = useState(false);
  const [alertMessage, setAlertMessage] = useState('');

  const {
    isRecording,
    currentNote,
    setCurrentNote,
    startRecording,
    stopRecording,
    showContinue,
    error,
    audioLevel,
    isEditing,
    setShowContinue,
    setWasLastActionSave, // Add setWasLastActionSave from useVoiceRecognition
  } = useVoiceRecognition();

  const {
    notes,
    sentiment,
    summary,
    editingNoteId,
    saveNote,
    editNote,
    updateNote,
    cancelEditing,
    isProcessing,
    hasUnsavedChanges,
    addTag,
    exportNotes,
    tags,
  } = useNotes(currentNote, setCurrentNote);

  // Setup speech commands
  useSpeechCommands({
    isRecording,
    saveNote,
    cancelEditing,
    currentNote,
    setCurrentNote
  });

  // Toggle theme
  const toggleTheme = () => {
    const newTheme = theme === 'light' ? 'dark' : 'light';
    setTheme(newTheme);
    localStorage.setItem('voice_notes_theme', newTheme);
    document.documentElement.classList.toggle('dark', newTheme === 'dark');
  };

  // Set theme class on document
  useEffect(() => {
    document.documentElement.classList.toggle('dark', theme === 'dark');
  }, [theme]);

  // Show alert
  const displayAlert = (message) => {
    setAlertMessage(message);
    setShowAlert(true);
    setTimeout(() => setShowAlert(false), 3000);
  };

  // Handle save note with feedback
  const handleSaveNote = async () => {
    if (isRecording) {
      stopRecording();
    }
    const result = await saveNote();
    if (result) {
      displayAlert('Note saved successfully!');
      localStorage.removeItem('voice_notes_draft'); // Clear draft
      setCurrentNote('');
      setShowContinue(false); // Disable continue recording
      setWasLastActionSave(true); // Mark last action as save
    } else {
      displayAlert('Failed to save note. Please try again.');
    }
  };

  // Unsaved changes warning
  useEffect(() => {
    const handleBeforeUnload = (e) => {
      if (hasUnsavedChanges && currentNote.trim()) {
        e.preventDefault();
        e.returnValue = '';
        return '';
      }
    };

    window.addEventListener('beforeunload', handleBeforeUnload);
    return () => window.removeEventListener('beforeunload', handleBeforeUnload);
  }, [hasUnsavedChanges, currentNote]);

  return (
    <div className={`min-h-screen transition-colors duration-300 ${
      theme === 'light' 
        ? 'bg-gradient-to-br from-blue-50 to-purple-50 text-gray-800' 
        : 'bg-gradient-to-br from-gray-900 to-indigo-900 text-gray-100'
    }`}>
      <div className="max-w-4xl mx-auto p-4 md:p-6">
        <div className={`rounded-xl overflow-hidden shadow-xl ${
          theme === 'light' ? 'bg-white' : 'bg-gray-800'
        }`}>
          <header className={`px-6 py-4 border-b ${
            theme === 'light' ? 'border-gray-200' : 'border-gray-700'
          } flex justify-between items-center`}>
            <div className="flex items-center gap-2">
              {isRecording ? (
                <Mic className="h-5 w-5 text-red-500 animate-pulse" />
              ) : (
                <MicOff className="h-5 w-5" />
              )}
              <h1 className={`text-xl font-bold bg-gradient-to-r ${
                theme === 'light'
                  ? 'from-purple-600 to-blue-600'
                  : 'from-purple-400 to-blue-400'
              } bg-clip-text text-transparent`}>
                Voice-to-Text Notes
              </h1>
            </div>
            
            {/* <button 
              onClick={toggleTheme}
              className={`p-2 rounded-full transition-colors ${
                theme === 'light'
                  ? 'hover:bg-gray-100 text-gray-700'
                  : 'hover:bg-gray-700 text-gray-300'
              }`}
              aria-label="Toggle theme"
            >
              {theme === 'light' ? <Moon size={18} /> : <Sun size={18} />}
            </button> */}
          </header>

          {showAlert && (
            <div className={`mx-4 mt-4 px-4 py-3 rounded-md flex items-center gap-2 ${
              alertMessage.includes('Failed') 
                ? 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400' 
                : 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400'
            }`}>
              <AlertCircle size={16} />
              <span>{alertMessage}</span>
            </div>
          )}

          <RecordingControls 
            isRecording={isRecording}
            startRecording={() => startRecording(currentNote, false)}
            stopRecording={stopRecording}
            saveNote={handleSaveNote}
            showContinue={showContinue}
            error={error}
            isProcessing={isProcessing}
            hasUnsavedChanges={hasUnsavedChanges}
            currentNote={currentNote}
            audioLevel={audioLevel}
            theme={theme}
            isEditing={isEditing}
          />

          <CurrentNote 
            currentNote={currentNote}
            setCurrentNote={setCurrentNote}
            sentiment={sentiment}
            summary={summary}
            isEditing={!!editingNoteId}
            theme={theme}
          />

          <NotesList 
            notes={notes}
            editNote={(noteId) => {
              if (isRecording) stopRecording();
              editNote(noteId);
            }}
            editingNoteId={editingNoteId}
            updateNote={updateNote}
            cancelEditing={cancelEditing}
            startEditingWithRecording={(noteId) => {
              if (isRecording) stopRecording();
              const note = notes.find(n => n._id === noteId);
              if (note) {
                editNote(noteId);
                startRecording(note.content, true);
              }
            }}
            addTag={addTag}
            exportNotes={exportNotes}
            tags={tags}
            theme={theme}
            currentNote={currentNote}
            setCurrentNote={setCurrentNote}
            stopRecording={stopRecording}
            setShowContinue={setShowContinue}
            setWasLastActionSave={setWasLastActionSave} // Pass setWasLastActionSave to NotesList
          />
        </div>
      </div>
    </div>
  );
}

export default VoiceToTextNotes;