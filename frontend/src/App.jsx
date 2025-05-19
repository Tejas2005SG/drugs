// Update the routes in App.jsx to include the new Summary component
import Summary from "./pages/Summary/Summary.jsx";

// Add this to your routes array inside the Dashboard route
<Route
  path="summary"
  element={
    <ProtectedRoute>
      <Summary />
    </ProtectedRoute>
  }
/>

export default Summary