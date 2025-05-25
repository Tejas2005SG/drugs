// Add this import with other route imports
import SymptomChecker from "./pages/SymptomChecker/SymptomChecker.jsx";

// Add this route inside your Routes component
<Route
  path="dashboard/symptom-checker"
  element={
    <ProtectedRoute>
      <SymptomChecker />
    </ProtectedRoute>
  }
/>