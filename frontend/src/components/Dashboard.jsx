// Add this to your navElements array
{
  name: 'Symptom Checker',
  icon: <Activity size={20} className="mr-3" />,
  navigation: () => navigate('/dashboard/symptom-checker'),
  roles: ['admin', 'citizen', 'guest']
},