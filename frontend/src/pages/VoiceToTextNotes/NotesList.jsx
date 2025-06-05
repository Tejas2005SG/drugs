import React, { useState } from 'react';
import { Edit, Mic, Save, X, Tag, Download, Search } from 'lucide-react';

const NotesList = ({
  notes,
  editNote,
  editingNoteId,
  updateNote,
  cancelEditing,
  startEditingWithRecording,
  addTag,
  exportNotes,
  tags,
  theme,
  currentNote,
  setCurrentNote,
  stopRecording,
  setShowContinue, // Add setShowContinue to props
}) => {
  const [searchTerm, setSearchTerm] = useState('');
  const [selectedTag, setSelectedTag] = useState('');
  const [newTag, setNewTag] = useState('');

  // Filter notes based on search term and selected tag
  const filteredNotes = notes.filter(note => {
    const matchesSearch = note.content.toLowerCase().includes(searchTerm.toLowerCase());
    const matchesTag = selectedTag ? (note.tags && note.tags.includes(selectedTag)) : true;
    return matchesSearch && matchesTag;
  });

  const handleAddTag = (noteId) => {
    if (newTag.trim()) {
      addTag(noteId, newTag.trim());
      setNewTag('');
    }
  };

  const formatDate = (dateString) => {
    if (!dateString || isNaN(new Date(dateString))) {
      return 'Unknown Date';
    }
    const date = new Date(dateString);
    return new Intl.DateTimeFormat('en-US', {
      month: 'short',
      day: 'numeric',
      year: 'numeric',
      hour: '2-digit',
      minute: '2-digit'
    }).format(date);
  };

  return (
    <div className="p-4 md:p-6">
      <div className="flex flex-wrap justify-between items-center mb-6 gap-4">
        <h2 className={`text-lg font-semibold ${
          theme === 'light' ? 'text-gray-700' : 'text-gray-200'
        }`}>
          Your Notes
        </h2>

        <div className="flex flex-wrap items-center gap-3">
          <div className={`relative flex items-center rounded-md ${
            theme === 'light' ? 'bg-gray-100' : 'bg-gray-700'
          }`}>
            <Search size={16} className="absolute left-3 text-gray-400" />
            <input
              type="text"
              placeholder="Search notes..."
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className={`py-2 pl-9 pr-3 rounded-md w-full sm:w-auto ${
                theme === 'light' 
                  ? 'bg-gray-100 text-gray-800 placeholder-gray-400 focus:ring-indigo-300 focus:border-indigo-500' 
                  : 'bg-gray-700 text-gray-100 placeholder-gray-400 focus:ring-indigo-500'
              } border-none focus:outline-none focus:ring-2`}
            />
          </div>

          <select
            value={selectedTag}
            onChange={(e) => setSelectedTag(e.target.value)}
            className={`py-2 px-3 rounded-md border ${
              theme === 'light'
                ? 'bg-white border-gray-300 text-gray-700'
                : 'bg-gray-700 border-gray-600 text-gray-200'
            }`}
          >
            <option value="">All Tags</option>
            {tags.map(tag => (
              <option key={tag} value={tag}>{tag}</option>
            ))}
          </select>

          <button
            onClick={() => exportNotes(filteredNotes)}
            className={`flex items-center gap-1 py-2 px-3 rounded-md ${
              theme === 'light'
                ? 'bg-gray-100 hover:bg-gray-200 text-gray-700'
                : 'bg-gray-700 hover:bg-gray-600 text-gray-200'
            } transition-colors`}
            title="Export notes"
          >
            <Download size={16} />
            <span className="hidden sm:inline">Export</span>
          </button>
        </div>
      </div>

      {filteredNotes.length === 0 && (
        <div className={`text-center py-8 ${
          theme === 'light' ? 'text-gray-500' : 'text-gray-400'
        }`}>
          {searchTerm || selectedTag ? 'No matching notes found.' : 'No notes saved yet.'}
        </div>
      )}

      <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
        {filteredNotes.map((note) => (
          <div 
            key={note._id} 
            className={`rounded-lg overflow-hidden border ${
              editingNoteId === note._id 
                ? theme === 'light' ? 'border-indigo-300 ring-2 ring-indigo-200' : 'border-indigo-500 ring-2 ring-indigo-500/30'
                : theme === 'light' ? 'border-gray-200 hover:border-gray-300' : 'border-gray-700 hover:border-gray-600'
            } transition-all`}
          >
            {editingNoteId === note._id ? (
              <div className={`p-4 ${
                theme === 'light' ? 'bg-white' : 'bg-gray-800'
              }`}>
                <textarea
                  value={currentNote}
                  onChange={(e) => setCurrentNote(e.target.value)}
                  className={`w-full min-h-[120px] p-3 mb-3 border rounded resize-y ${
                    theme === 'light'
                      ? 'bg-white border-gray-300 text-gray-800 focus:ring-indigo-300 focus:border-indigo-500'
                      : 'bg-gray-700 border-gray-600 text-gray-100 focus:ring-indigo-500/30 focus:border-indigo-500'
                  } focus:outline-none focus:ring-2`}
                ></textarea>

                <div className="flex items-center gap-2 mb-3">
                  <input
                    type="text"
                    placeholder="Add tag..."
                    value={newTag}
                    onChange={(e) => setNewTag(e.target.value)}
                    className={`flex-1 p-2 text-sm border rounded ${
                      theme === 'light'
                        ? 'bg-white border-gray-300 text-gray-800 focus:ring-indigo-300 focus:border-indigo-500'
                        : 'bg-gray-700 border-gray-600 text-gray-100 focus:ring-indigo-500/30 focus:border-indigo-500'
                    } focus:outline-none focus:ring-2`}
                  />
                  <button
                    onClick={() => handleAddTag(note._id)}
                    className={`p-2 rounded ${
                      !newTag.trim()
                        ? 'opacity-50 cursor-not-allowed'
                        : theme === 'light'
                          ? 'bg-indigo-100 text-indigo-700 hover:bg-indigo-200'
                          : 'bg-indigo-900/50 text-indigo-400 hover:bg-indigo-900/70'
                    } transition-colors`}
                    disabled={!newTag.trim()}
                  >
                    <Tag size={16} />
                  </button>
                </div>

                <div className="flex gap-2">
                  <button
                    onClick={async () => {
                      stopRecording(); // Stop any active recording
                      const success = await updateNote(note._id, currentNote);
                      if (success) {
                        cancelEditing();
                        setShowContinue(false); // Prevent "Continue Recording" button
                        localStorage.removeItem('voice_notes_draft'); // Clear draft to avoid continue option
                      }
                    }}
                    className={`flex items-center gap-1 py-1 px-3 rounded font-medium ${
                      theme === 'light'
                        ? 'bg-green-100 text-green-700 hover:bg-green-200'
                        : 'bg-green-900/30 text-green-400 hover:bg-green-900/50'
                    } transition-colors`}
                  >
                    <Save size={14} />
                    <span>Save</span>
                  </button>

                  <button
                    onClick={() => startEditingWithRecording(note._id)}
                    className={`flex items-center gap-1 py-1 px-3 rounded font-medium ${
                      theme === 'light'
                        ? 'bg-blue-100 text-blue-700 hover:bg-blue-200'
                        : 'bg-blue-900/30 text-blue-400 hover:bg-blue-900/50'
                    } transition-colors`}
                  >
                    <Mic size={14} />
                    <span>Record</span>
                  </button>

                  <button
                    onClick={cancelEditing}
                    className={`flex items-center gap-1 py-1 px-3 rounded font-medium ${
                      theme === 'light'
                        ? 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                        : 'bg-gray-700 text-gray-300 hover:bg-gray-600'
                    } transition-colors`}
                  >
                    <X size={14} />
                    <span>Cancel</span>
                  </button>
                </div>
              </div>
            ) : (
              <div className={`p-4 h-full flex flex-col ${
                theme === 'light' ? 'bg-white' : 'bg-gray-800'
              }`}>
                <div className="mb-2 flex-grow">
                  <p className={theme === 'light' ? 'text-gray-800' : 'text-gray-200'}>
                    {note.content.length > 200 
                      ? `${note.content.substring(0, 200)}...` 
                      : note.content}
                  </p>
                </div>

                <div className={`text-xs mb-2 ${
                  theme === 'light' ? 'text-gray-500' : 'text-gray-400'
                }`}>
                  {formatDate(note.createdAt)}
                </div>

                {/* {note.sentiment && (
                  <div className={`text-xs mb-1 ${
                    note.sentiment === 'positive' 
                      ? 'text-green-500' 
                      : note.sentiment === 'negative' 
                        ? 'text-red-500' 
                        : 'text-yellow-500'
                  }`}>
                    Sentiment: {note.sentiment}
                  </div>
                )} */}

                {note.tags && note.tags.length > 0 && (
                  <div className="flex flex-wrap gap-1 mb-2">
                    {note.tags.map(tag => (
                      <span 
                        key={tag} 
                        className={`text-xs px-2 py-0.5 rounded-full ${
                          theme === 'light'
                            ? 'bg-indigo-100 text-indigo-700'
                            : 'bg-indigo-900/50 text-indigo-300'
                        }`}
                      >
                        {tag}
                      </span>
                    ))}
                  </div>
                )}

                {/* {note.summary && (
                  <div className={`text-xs mb-3 ${
                    theme === 'light' ? 'text-gray-600' : 'text-gray-300'
                  }`}>
                    <span className="font-medium">Summary:</span> {note.summary}
                  </div>
                )} */}

                <div className="flex gap-2 mt-auto pt-2">
                  <button
                    onClick={() => editNote(note._id)}
                    className={`p-2 rounded ${
                      theme === 'light'
                        ? 'bg-blue-50 text-blue-600 hover:bg-blue-100'
                        : 'bg-blue-900/20 text-blue-400 hover:bg-blue-900/30'
                    } transition-colors`}
                    title="Edit note"
                  >
                    <Edit size={16} />
                  </button>

                  <button
                    onClick={() => startEditingWithRecording(note._id)}
                    className={`p-2 rounded ${
                      theme === 'light'
                        ? 'bg-indigo-50 text-indigo-600 hover:bg-indigo-100'
                        : 'bg-indigo-900/20 text-indigo-400 hover:bg-indigo-900/30'
                    } transition-colors`}
                    title="Continue recording"
                  >
                    <Mic size={16} />
                  </button>
                </div>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
};

export default NotesList;