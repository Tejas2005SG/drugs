// Add Summary to the navElements array in Dashboard.jsx
{
  name: 'Summary',
  icon: <FileDown size={20} className="mr-3" />,
  navigation: () => navigate('/dashboard/summary'),
  roles: ['admin', 'citizen', 'guest']
},