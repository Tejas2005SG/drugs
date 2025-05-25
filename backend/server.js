// Add this line with the other route imports
import symptomsRoutes from './routes/symptoms.routes.js';

// Add this line with the other app.use statements
app.use('/api/symptoms', symptomsRoutes);