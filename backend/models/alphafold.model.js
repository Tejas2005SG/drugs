// models/alphafold.model.js
import mongoose from 'mongoose';

const JobSchema = new mongoose.Schema({
  uniprot_id: {  // Changed from sequence to uniprot_id
    type: String,
    required: [true, 'UniProt ID is required'],
    trim: true,
    validate: {
      validator: (v) => /^[A-Z0-9]{6,10}$/i.test(v),
      message: 'Invalid UniProt ID format'
    }
  },
  status: {
    type: String,
    enum: ['pending', 'processing', 'completed', 'failed'],
    default: 'pending'
  },
  pythonJobId: {  // Renamed to internalJobId for clarity
    type: String,
    index: true,
    sparse: true
  },
  pdbUrl: {
    type: String,
    default: null
  },
  error: {
    type: String,
    default: null
  },
  createdAt: {
    type: Date,
    default: Date.now
  },
  updatedAt: {
    type: Date,
    default: Date.now
  },
  completedAt: {
    type: Date,
    default: null
  }
}, {
  timestamps: false,
  versionKey: false
});

JobSchema.pre('save', function(next) {
  this.updatedAt = new Date();
  if (this.isModified('status') && this.status === 'completed') {
    this.completedAt = new Date();
  }
  next();
});

JobSchema.index({ status: 1 });
JobSchema.index({ createdAt: -1 });

export const Job = mongoose.model('Job', JobSchema);