import mongoose from 'mongoose';

const NoteSchema = new mongoose.Schema({    
    title: {
        
        type: String,
        required: true,
        trim: true
    },    
    description: {
        type: String,
        required: true,
        trim: true
    },
    userId: {
        type: mongoose.Schema.Types.ObjectId,
        ref: 'User',
        required: true
    },
    createdAt: {
        type: Date,
        default: Date.now
    }
});

const Note = mongoose.model('Note', NoteSchema);      

export default Note;