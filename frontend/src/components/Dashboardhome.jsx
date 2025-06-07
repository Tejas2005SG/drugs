import React from 'react';
import { useAuthStore } from '../store/auth.store';
import { Link } from 'react-router-dom';
import {
  Box,
  Typography,
  Card,
  CardContent,
  Grid,
  Button,
  Stack,
  IconButton,
  useTheme,
  Chip,
  CardActions,
} from '@mui/material';
import {
  Menu as MenuIcon,
  Home,
  Message as MessageSquare,
  Biotech,
  DirectionsRun as Activity,
  Layers,
  AttachMoney,
  Psychology,
  TrendingUp,
  ArrowForward,
  Description,
  CenterFocusStrong as Target,
  Science,
  FileDownload,
  KeyboardVoice as RecordVoiceOver,
  Newspaper,
  ChevronLeft,
  Logout,
  Dashboard,
  MonetizationOn,
  Assignment,
  Warning,
} from '@mui/icons-material';
import { Toaster } from 'react-hot-toast';

function DashboardHome() {
  const { user } = useAuthStore();
  const theme = useTheme();

  const quickActions = [
    {
      title: 'Alphafold Structure Prediction',
      description: 'Predict and visualize protein structures with atom and bond',
      icon: <Science sx={{ fontSize: 28 }} />,
      link: '/dashboard/getalphafoldstrcture',
      gradient: 'linear-gradient(135deg, #00C9B7, #00F5D4)',
    },
    {
      title: 'Protein Structure Generation',
      description: 'Generate and analyze protein structures',
      icon: <Biotech sx={{ fontSize: 28 }} />,
      link: '/dashboard/protein-structure',
      gradient: 'linear-gradient(135deg, #00A8C5, #00D4F5)',
    },
    {
      title: 'New Drug Discovery',
      description: 'Explore New Drugs Discovery',
      icon: <Psychology sx={{ fontSize: 28 }} />,
      link: '/dashboard/protein-structure-mutation',
      gradient: 'linear-gradient(135deg, #009B95, #00C9B7)',
    },
    {
      title: 'Cost Estimation',
      description: 'Estimate costs for drug development',
      icon: <MonetizationOn sx={{ fontSize: 28 }} />,
      link: '/dashboard/cost-estimation',
      gradient: 'linear-gradient(135deg, #008B8F, #00A8C5)',
    },
    {
      title: 'AI Research Paper Generator',
      description: 'Generate research papers using AI',
      icon: <Assignment sx={{ fontSize: 28 }} />,
      link: '/dashboard/ai-research-paper-generator',
      gradient: 'linear-gradient(135deg, #007A7A, #009B95)',
    },
    {
      title: 'AI Driven Target Prediction',
      description: 'Predict drug targets using AI',
      icon: <Target sx={{ fontSize: 28 }} />,
      link: '/dashboard/ai-driven-target-prediction',
      gradient: 'linear-gradient(135deg, #006666, #008B8F)',
    },
    {
      title: 'AI Naming Suggestions',
      description: 'Let AI suggest the name for your New Drug',
      icon: <Psychology sx={{ fontSize: 28 }} />,
      link: '/dashboard/ai-naming',
      gradient: 'linear-gradient(135deg, #005151, #007A7A)',
    },
    {
      title: 'Toxicity Prediction',
      description: 'Get to know about New Drug',
      icon: <Warning sx={{ fontSize: 28 }} />,
      link: '/dashboard/toxicityPrediction',
      gradient: 'linear-gradient(135deg, #003D3D, #006666)',
    },
    {
      title: 'Live News',
      description: 'Stay updated with the latest news',
      icon: <Newspaper sx={{ fontSize: 28 }} />,
      link: '/dashboard/live-news',
      gradient: 'linear-gradient(135deg, #002A2A, #005151)',
    },
  ];

  return (
    <Box sx={{ 
      maxWidth: 1400, 
      mx: 'auto', 
      p: 3,
      background: 'radial-gradient(circle at 10% 20%, #0A192F 0%, #0D1B36 100%)',
      minHeight: '100vh'
    }}>
      {/* Header Section */}
      <Card
        sx={{
          mb: 4,
          background: 'linear-gradient(135deg, #0A192F, #172A45)',
          border: '1px solid rgba(0, 245, 212, 0.1)',
          borderRadius: 4,
          boxShadow: '0 8px 32px rgba(0, 245, 212, 0.1)',
          backdropFilter: 'blur(8px)',
          overflow: 'hidden',
          position: 'relative',
          '&:before': {
            content: '""',
            position: 'absolute',
            top: 0,
            left: 0,
            right: 0,
            height: 4,
            background: 'linear-gradient(90deg, #00F5D4, #5E81F4)',
          }
        }}
      >
        <CardContent sx={{ p: 4 }}>
          <Stack
            direction={{ xs: 'column', md: 'row' }}
            justifyContent="space-between"
            alignItems={{ xs: 'flex-start', md: 'center' }}
            spacing={3}
          >
            <Box>
              <Typography
                variant="h3"
                sx={{
                  fontWeight: 800,
                  mb: 1,
                  fontFamily: '"Barlow", sans-serif',
                  background: 'linear-gradient(135deg, #00F5D4, #5E81F4)',
                  backgroundClip: 'text',
                  WebkitBackgroundClip: 'text',
                  WebkitTextFillColor: 'transparent',
                  textShadow: '0 2px 10px rgba(0, 245, 212, 0.3)',
                  letterSpacing: '-0.5px'
                }}
              >
                Welcome back, {user?.firstName || 'Researcher'}!
              </Typography>
              <Typography 
                variant="h6" 
                sx={{
                  color: '#A0A0A0',
                  fontFamily: '"Roboto", sans-serif',
                  fontWeight: 400,
                }}
              >
                Your drug discovery dashboard is ready for today's research.
              </Typography>
            </Box>
            <Button
              component={Link}
              to="/"
              variant="outlined"
              startIcon={<Home />}
              sx={{
                borderRadius: 3,
                px: 4,
                py: 1.5,
                borderColor: 'rgba(0, 245, 212, 0.3)',
                color: '#00F5D4',
                fontFamily: '"Barlow", sans-serif',
                fontWeight: 500,
                letterSpacing: '0.5px',
                transition: 'all 0.3s ease',
                '&:hover': {
                  borderColor: '#00F5D4',
                  backgroundColor: 'rgba(0, 245, 212, 0.1)',
                  transform: 'translateY(-2px)',
                  boxShadow: '0 4px 15px rgba(0, 245, 212, 0.2)'
                },
              }}
            >
              Back to Home
            </Button>
          </Stack>
        </CardContent>
      </Card>

      {/* Quick Actions Section */}
      <Box sx={{ mb: 4 }}>
        <Typography
          variant="h4"
          sx={{
            fontWeight: 700,
            color: '#E0E0E0',
            mb: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            fontFamily: '"Barlow", sans-serif',
            letterSpacing: '-0.25px'
          }}
        >
          <TrendingUp sx={{ 
            color: '#00F5D4',
            fontSize: 36
          }} />
          Quick Actions
        </Typography>
        
        <Grid container spacing={3}>
          {quickActions.map((action, index) => (
            <Grid item xs={12} sm={6} lg={4} key={index}>
              <Card
                sx={{
                  height: '100%',
                  background: 'linear-gradient(135deg, #172A45, #0A192F)',
                  border: '1px solid rgba(0, 245, 212, 0.1)',
                  borderRadius: 3,
                  transition: 'all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275)',
                  cursor: 'pointer',
                  boxShadow: '0 4px 20px rgba(0, 0, 0, 0.2)',
                  '&:hover': {
                    transform: 'translateY(-8px)',
                    boxShadow: '0 12px 35px rgba(0, 245, 212, 0.25)',
                    borderColor: 'rgba(0, 245, 212, 0.3)',
                  },
                }}
                component={Link}
                to={action.link}
                style={{ textDecoration: 'none' }}
              >
                <CardContent sx={{ 
                  p: 3, 
                  height: '100%', 
                  display: 'flex', 
                  flexDirection: 'column',
                  position: 'relative',
                  overflow: 'hidden',
                  '&:before': {
                    content: '""',
                    position: 'absolute',
                    top: 0,
                    left: 0,
                    right: 0,
                    height: 4,
                    background: action.gradient,
                  }
                }}>
                  <Box
                    sx={{
                      width: 56,
                      height: 56,
                      borderRadius: 2,
                      background: action.gradient,
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      mb: 2,
                      color: '#0A192F',
                      mt: 2,
                      boxShadow: '0 4px 15px rgba(0, 245, 212, 0.3)'
                    }}
                  >
                    {action.icon}
                  </Box>
                  
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: '#E0E0E0',
                      mb: 1,
                      lineHeight: 1.3,
                      fontFamily: '"Barlow", sans-serif',
                      fontSize: '1.1rem'
                    }}
                  >
                    {action.title}
                  </Typography>
                  
                  <Typography
                    variant="body2"
                    sx={{ 
                      mb: 2, 
                      flex: 1,
                      color: '#A0A0A0',
                      fontFamily: '"Roboto", sans-serif',
                      lineHeight: 1.5
                    }}
                  >
                    {action.description}
                  </Typography>
                  
                  <Stack direction="row" justifyContent="space-between" alignItems="center">
                    <Chip
                      label="Available"
                      size="small"
                      sx={{
                        backgroundColor: 'rgba(112, 224, 0, 0.1)',
                        color: '#70E000',
                        fontWeight: 500,
                        fontFamily: '"Roboto", sans-serif',
                        fontSize: '0.75rem'
                      }}
                    />
                    <IconButton
                      size="small"
                      sx={{
                        color: '#00F5D4',
                        backgroundColor: 'rgba(0, 245, 212, 0.1)',
                        transition: 'all 0.3s ease',
                        '&:hover': {
                          backgroundColor: 'rgba(0, 245, 212, 0.2)',
                          transform: 'translateX(2px)'
                        },
                      }}
                    >
                      <ArrowForward fontSize="small" />
                    </IconButton>
                  </Stack>
                </CardContent>
              </Card>
            </Grid>
          ))}
        </Grid>
      </Box>

      <Toaster 
        position="bottom-right"
        toastOptions={{
          style: {
            background: '#172A45',
            color: '#E0E0E0',
            border: '1px solid rgba(0, 245, 212, 0.3)',
            borderRadius: '12px',
            fontFamily: '"Roboto", sans-serif'
          },
        }}
      />
    </Box>
  );
}

export default DashboardHome;