import React, { useEffect, useState } from 'react';
import { useAuthStore } from '../Store/auth.store.js';
import { Outlet, Link, useNavigate, useLocation } from 'react-router-dom';
import {
  UserPlus, LogIn, LogOut, Menu, X, ChevronRight, Activity, BrainCog,
  Settings, Home, Layers, Dna, DollarSign, FileText, Target, Pill,
  Newspaper, MessageSquare, Torus, FileDown, ChevronLeft, 
 
  
      FlaskConical,

  
  AlertTriangle,
  Mic,
  FileBox,

  Atom
  
} from 'lucide-react';
import { Toaster } from 'react-hot-toast';
import {
  AppBar,
  Toolbar,
  Drawer,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  IconButton,
  Typography,
  Avatar,
  Box,
  Divider,
  Button,
  Tooltip,
  useMediaQuery,
  useTheme,
  Chip,
  Fade,
  GlobalStyles
} from '@mui/material';
import { createTheme, ThemeProvider, styled } from '@mui/material/styles';
import ToxicityPrediction from '../pages/ToxicityPrediction/ToxicityPrediction.jsx';
import JarvisBot from '../pages/Chatbot/Chatbot.jsx';

// Global styles for custom scrollbar
const globalStyles = (
  <GlobalStyles
    styles={{
      '*': {
        scrollbarWidth: 'thin',
        scrollbarColor: '#00F5D4 rgba(0, 245, 212, 0.1)',
      },
      '*::-webkit-scrollbar': {
        width: '8px',
        height: '8px',
      },
      '*::-webkit-scrollbar-track': {
        background: 'rgba(0, 245, 212, 0.05)',
        borderRadius: '10px',
      },
      '*::-webkit-scrollbar-thumb': {
        background: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
        borderRadius: '10px',
        border: '2px solid transparent',
        backgroundClip: 'padding-box',
      },
      '*::-webkit-scrollbar-thumb:hover': {
        background: 'linear-gradient(135deg, #5E81F4 0%, #00F5D4 100%)',
      },
      '*::-webkit-scrollbar-corner': {
        background: 'transparent',
      },
      // Remove horizontal scrollbar
      'html, body': {
        overflowX: 'hidden',
        maxWidth: '100vw',
      },
      // Ensure no horizontal overflow
      '*': {
        boxSizing: 'border-box',
      },
    }}
  />
);

// Custom theme based on your specifications
const customTheme = createTheme({
  palette: {
    mode: 'dark',
    primary: {
      main: '#00F5D4', // AI Teal
      light: '#5E81F4', // Soft Blue
      dark: '#0A192F', // Navy
    },
    secondary: {
      main: '#5E81F4', // Soft Blue
      light: '#00F5D4',
      dark: '#172A45',
    },
    background: {
      default: '#0A192F', // Primary BG
      paper: '#172A45', // Secondary BG
    },
    text: {
      primary: '#E0E0E0', // Primary text
      secondary: '#A0A0A0', // Secondary text
    },
    error: {
      main: '#FF4D6D', // Coral Pink
    },
    success: {
      main: '#70E000', // Vibrant Green
    },
  },
  typography: {
    fontFamily: '"Inter", "Roboto", "Helvetica", "Arial", sans-serif',
    h1: {
      fontFamily: '"Inter", sans-serif',
      fontWeight: 700,
    },
    h2: {
      fontFamily: '"Inter", sans-serif',
      fontWeight: 600,
    },
    h3: {
      fontFamily: '"Inter", sans-serif',
      fontWeight: 600,
    },
    body1: {
      fontFamily: '"Roboto", sans-serif',
    },
    body2: {
      fontFamily: '"Roboto", sans-serif',
    },
  },
  components: {
    MuiDrawer: {
      styleOverrides: {
        paper: {
          background: 'linear-gradient(135deg, #0A192F 0%, #172A45 100%)',
          borderRight: '1px solid rgba(0, 245, 212, 0.1)',
        },
      },
    },
    MuiAppBar: {
      styleOverrides: {
        root: {
          background: 'linear-gradient(90deg, #0A192F 0%, #172A45 100%)',
          backdropFilter: 'blur(10px)',
          borderBottom: '1px solid rgba(0, 245, 212, 0.1)',
        },
      },
    },
    MuiListItem: {
      styleOverrides: {
        root: {
          borderRadius: '12px',
          margin: '4px 8px',
          transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
          '&:hover': {
            backgroundColor: 'rgba(0, 245, 212, 0.08)',
            transform: 'translateX(4px)',
            boxShadow: '0 4px 20px rgba(0, 245, 212, 0.15)',
          },
          '&.active': {
            backgroundColor: 'rgba(0, 245, 212, 0.12)',
            borderLeft: '3px solid #00F5D4',
            '& .MuiListItemIcon-root': {
              color: '#00F5D4',
            },
            '& .MuiListItemText-primary': {
              color: '#00F5D4',
              fontWeight: 600,
            },
          },
        },
      },
    },
    MuiButton: {
      styleOverrides: {
        root: {
          borderRadius: '12px',
          textTransform: 'none',
          fontWeight: 600,
          padding: '10px 20px',
        },
      },
    },
  },
});

// Styled components
const StyledAppBar = styled(AppBar)(({ theme }) => ({
  zIndex: theme.zIndex.drawer + 1,
  background: 'linear-gradient(90deg, #0A192F 0%, #172A45 100%)',
  backdropFilter: 'blur(10px)',
  borderBottom: '1px solid rgba(0, 245, 212, 0.1)',
  boxShadow: '0 4px 20px rgba(0, 0, 0, 0.3)',
  width: '100%',
  maxWidth: '100vw',
}));

const StyledDrawer = styled(Drawer)(({ theme }) => ({
  '& .MuiDrawer-paper': {
    background: 'linear-gradient(135deg, #0A192F 0%, #172A45 100%)',
    borderRight: '1px solid rgba(0, 245, 212, 0.1)',
    backdropFilter: 'blur(10px)',
    overflowX: 'hidden',
    maxWidth: '100vw',
    // Custom scrollbar for drawer
    '&::-webkit-scrollbar': {
      width: '6px',
    },
    '&::-webkit-scrollbar-track': {
      background: 'rgba(0, 245, 212, 0.05)',
      borderRadius: '10px',
    },
    '&::-webkit-scrollbar-thumb': {
      background: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
      borderRadius: '10px',
    },
    '&::-webkit-scrollbar-thumb:hover': {
      background: 'linear-gradient(135deg, #5E81F4 0%, #00F5D4 100%)',
    },
  },
}));

const UserSection = styled(Box)(({ theme }) => ({
  background: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
  padding: theme.spacing(3),
  margin: theme.spacing(2),
  borderRadius: '16px',
  boxShadow: '0 8px 32px rgba(0, 245, 212, 0.2)',
  color: '#0A192F',
  position: 'relative',
  overflow: 'hidden',
  '&::before': {
    content: '""',
    position: 'absolute',
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    background: 'linear-gradient(135deg, rgba(255,255,255,0.1) 0%, rgba(255,255,255,0.05) 100%)',
    borderRadius: '16px',
  },
}));

const CollapseButton = styled(IconButton)(({ theme }) => ({
  position: 'absolute',
  right: -12,
  top: 20,
  backgroundColor: '#00F5D4',
  color: '#0A192F',
  width: 24,
  height: 24,
  boxShadow: '0 4px 12px rgba(0, 245, 212, 0.3)',
  '&:hover': {
    backgroundColor: '#5E81F4',
    transform: 'scale(1.1)',
  },
  transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
}));

// Custom scrollable container
const ScrollableContainer = styled(Box)(({ theme }) => ({
  overflowY: 'auto',
  overflowX: 'hidden',
  maxHeight: '100%',
  paddingRight: '4px',
  '&::-webkit-scrollbar': {
    width: '6px',
  },
  '&::-webkit-scrollbar-track': {
    background: 'rgba(0, 245, 212, 0.05)',
    borderRadius: '10px',
    margin: '4px 0',
  },
  '&::-webkit-scrollbar-thumb': {
    background: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
    borderRadius: '10px',
    '&:hover': {
      background: 'linear-gradient(135deg, #5E81F4 0%, #00F5D4 100%)',
    },
  },
}));

const DRAWER_WIDTH = 280;
const COLLAPSED_WIDTH = 80;

function DashboardPage() {
  const { user, logout } = useAuthStore();
  const navigate = useNavigate();
  const location = useLocation();
  // const toastShown = useRef(false);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isCollapsed, setIsCollapsed] = useState(true);
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('md'));

  useEffect(() => {
    console.log('DashboardPage - User:', user);
  }, [user]);

  useEffect(() => {
    if (isMobile) {
      setIsSidebarOpen(false);
    }
  }, [location.pathname, isMobile]);

  
const navElements = [
  {
    name: "Dashboard Home",
    icon: <Home size={20} />,
    path: "/dashboard",
    navigation: () => navigate("/dashboard"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Discover New Drugs",
    icon: <FlaskConical size={20} />,
    path: "/dashboard/newdrug-discovery",
    navigation: () => navigate("/dashboard/newdrug-discovery"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Cost Estimation",
    icon: <DollarSign size={20} />,
    path: "/dashboard/cost-estimation",
    navigation: () => navigate("/dashboard/cost-estimation"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Protein Structure Generation",
    icon: <Dna size={20} />,
    path: "/dashboard/protein-structure",
    navigation: () => navigate("/dashboard/protein-structure"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Alphafold Structure",
    icon: <Atom size={20} />,
    path: "/dashboard/getalphafoldstrcture",
    navigation: () => navigate("/dashboard/getalphafoldstrcture"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "AI Naming Suggestion",
    icon: <BrainCog size={20} />,
    path: "/dashboard/ai-naming",
    navigation: () => navigate("/dashboard/ai-naming"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Research Paper Generation",
    icon: <FileText size={20} />,
    path: "/dashboard/ai-research-paper-generator",
    navigation: () => navigate("/dashboard/ai-research-paper-generator"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Side Effects Prediction",
    icon: <AlertTriangle size={20} />,
    path: "/dashboard/sideeffect-prediction",
    navigation: () => navigate("/dashboard/sideeffect-prediction"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Voice Notes",
    icon: <Mic size={20} />,
    path: "/dashboard/voice-text-notes",
    navigation: () => navigate("/dashboard/voice-text-notes"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Summary",
    icon: <FileBox size={20} />,
    path: "/dashboard/summary",
    navigation: () => navigate("/dashboard/summary"),
    roles: ["admin", "citizen", "guest"],
  },
  {
    name: "Live News",
    icon: <Newspaper size={20} />,
    path: "/dashboard/live-news",
    navigation: () => navigate("/dashboard/live-news"),
    roles: ["admin", "citizen", "guest"],
  },
];

  const filteredNavElements = navElements.filter(navElement => {
    const userRole = user?.role || 'guest';
    return navElement.roles.includes(userRole);
  });

  const handleLogout = () => {
    logout();
    navigate('/');
  };

  const toggleSidebar = () => {
    if (isMobile) {
      setIsSidebarOpen(!isSidebarOpen);
    } else {
      setIsCollapsed(!isCollapsed);
    }
  };

  const isActiveRoute = (path) => {
    return location.pathname === path;
  };

  const drawerWidth = isMobile ? DRAWER_WIDTH : (isCollapsed ? COLLAPSED_WIDTH : DRAWER_WIDTH);

  const DrawerContent = () => (
    <Box sx={{
      height: '100%',
      display: 'flex',
      flexDirection: 'column',
      position: 'relative',
      overflowX: 'hidden',
      maxWidth: '100%',
    }}>


      {/* User Section */}
      {(!isCollapsed || isMobile) && (
        <Fade in={!isCollapsed || isMobile} timeout={300}>
          <UserSection>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, position: 'relative', zIndex: 1 }}>
              <Avatar
                sx={{
                  bgcolor: '#0A192F',
                  color: '#00F5D4',
                  fontWeight: 'bold',
                  width: 48,
                  height: 48,
                  border: '2px solid rgba(255, 255, 255, 0.2)',
                }}
              >
                {user?.firstName?.charAt(0) || 'U'}
              </Avatar>
              <Box>
                <Typography variant="body2" sx={{ opacity: 0.8, fontSize: '0.875rem' }}>
                  Welcome,
                </Typography>
                <Typography variant="h6" sx={{ fontWeight: 700, fontSize: '1.1rem' }}>
                  {user ? user.firstName : 'Guest'}
                </Typography>
              </Box>
            </Box>
            {user?.role && (
              <Chip
                label={user.role.charAt(0).toUpperCase() + user.role.slice(1)}
                size="small"
                sx={{
                  mt: 1,
                  bgcolor: 'rgba(10, 25, 47, 0.3)',
                  color: '#0A192F',
                  fontWeight: 600,
                  position: 'relative',
                  zIndex: 1,
                }}
              />
            )}
          </UserSection>
        </Fade>
      )}

      {/* Collapsed User Avatar */}
      {isCollapsed && !isMobile && (
        <Box sx={{ display: 'flex', justifyContent: 'center', py: 2 }}>
          <Tooltip title={`Welcome, ${user?.firstName || 'Guest'}`} placement="right">
            <Avatar
              sx={{
                bgcolor: '#00F5D4',
                color: '#0A192F',
                fontWeight: 'bold',
                width: 40,
                height: 40,
                border: '2px solid rgba(0, 245, 212, 0.2)',
                boxShadow: '0 4px 12px rgba(0, 245, 212, 0.3)',
              }}
            >
              {user?.firstName?.charAt(0) || 'U'}
            </Avatar>
          </Tooltip>
        </Box>
      )}

      {/* Navigation List */}
      <ScrollableContainer sx={{ flex: 1 }}>
        <List sx={{ py: 1 }}>
          {filteredNavElements.map((navElement, index) => (
            <Tooltip
              key={index}
              title={isCollapsed && !isMobile ? navElement.name : ''}
              placement="right"
            >
              <ListItem
                button
                onClick={() => {
                  navElement.navigation();
                  if (isMobile) setIsSidebarOpen(false);
                }}
                className={isActiveRoute(navElement.path) ? 'active' : ''}
                sx={{
                  minHeight: 48,
                  justifyContent: isCollapsed ? 'center' : 'initial',
                  px: isCollapsed ? 0 : 2,
                  mx: 1,
                  mb: 0.5,
                  borderRadius: '12px',
                  transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
                  '&:hover': {
                    backgroundColor: 'rgba(0, 245, 212, 0.08)',
                    transform: 'translateX(4px)',
                    boxShadow: '0 4px 20px rgba(0, 245, 212, 0.15)',
                  },
                  '&.active': {
                    backgroundColor: 'rgba(0, 245, 212, 0.12)',
                    borderLeft: '3px solid #00F5D4',
                    '& .MuiListItemIcon-root': {
                      color: '#00F5D4',
                    },
                    '& .MuiListItemText-primary': {
                      color: '#00F5D4',
                      fontWeight: 600,
                    },
                  },
                }}
              >
                <ListItemIcon
                  sx={{
                    minWidth: 0,
                    mr: isCollapsed ? 0 : 2,
                    justifyContent: 'center',
                    color: isActiveRoute(navElement.path) ? '#00F5D4' : '#A0A0A0',
                  }}
                >
                  {navElement.icon}
                </ListItemIcon>
                {(!isCollapsed || isMobile) && (
                  <ListItemText
                    primary={navElement.name}
                    primaryTypographyProps={{
                      fontSize: '0.95rem',
                      fontWeight: isActiveRoute(navElement.path) ? 600 : 400,
                      color: isActiveRoute(navElement.path) ? '#00F5D4' : '#E0E0E0',
                    }}
                  />
                )}
                {(!isCollapsed || isMobile) && !isActiveRoute(navElement.path) && (
                  <ChevronRight size={16} style={{ color: '#A0A0A0', marginLeft: 'auto' }} />
                )}
              </ListItem>
            </Tooltip>
          ))}
        </List>
      </ScrollableContainer>

      <Divider sx={{ bgcolor: 'rgba(0, 245, 212, 0.1)' }} />

      {/* Logout Section */}
      <Box sx={{ p: 2 }}>
        {user ? (
          <Button
            fullWidth
            variant="outlined"
            color="error"
            onClick={handleLogout}
            startIcon={!isCollapsed || isMobile ? <LogOut size={18} /> : null}
            sx={{
              borderColor: '#FF4D6D',
              color: '#FF4D6D',
              '&:hover': {
                borderColor: '#FF4D6D',
                backgroundColor: 'rgba(255, 77, 109, 0.1)',
              },
              minWidth: isCollapsed && !isMobile ? 40 : 'auto',
              px: isCollapsed && !isMobile ? 1 : 2,
            }}
          >
            {isCollapsed && !isMobile ? <LogOut size={18} /> : 'Log Out'}
          </Button>
        ) : (
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
            <Button
              component={Link}
              to="/signup"
              variant="contained"
              color="primary"
              startIcon={!isCollapsed || isMobile ? <UserPlus size={18} /> : null}
              fullWidth
              sx={{ minWidth: isCollapsed && !isMobile ? 40 : 'auto' }}
            >
              {(!isCollapsed || isMobile) ? 'Sign Up' : <UserPlus size={18} />}
            </Button>
            <Button
              component={Link}
              to="/login"
              variant="outlined"
              color="secondary"
              startIcon={!isCollapsed || isMobile ? <LogIn size={18} /> : null}
              fullWidth
              sx={{ minWidth: isCollapsed && !isMobile ? 40 : 'auto' }}
            >
              {(!isCollapsed || isMobile) ? 'Login' : <LogIn size={18} />}
            </Button>
          </Box>
        )}
      </Box>
    </Box>
  );

  return (
    <ThemeProvider theme={customTheme}>
      {globalStyles}
      <Box sx={{
        display: 'flex',
        bgcolor: 'background.default',
        minHeight: '100vh',
        overflowX: 'hidden',
        maxWidth: '100vw',
      }}>
        {/* App Bar */}
        <StyledAppBar position="fixed">
          <Toolbar>
            <IconButton
              color="inherit"
              aria-label="toggle drawer"
              onClick={toggleSidebar}
              edge="start"
              sx={{
                marginRight: 2,
                color: '#00F5D4',
                '&:hover': { backgroundColor: 'rgba(0, 245, 212, 0.1)' },
              }}
            >
              {isMobile && isSidebarOpen ? (
                <X size={24} /> // Lucide's X icon
              ) : (
                <Menu size={24} /> // Lucide's Menu icon
              )}
            </IconButton>
            <Typography
              variant="h6"
              noWrap
              component="div"
              sx={{
                fontWeight: 700,
                background: 'linear-gradient(45deg, #00F5D4 30%, #5E81F4 90%)',
                backgroundClip: 'text',
                WebkitBackgroundClip: 'text',
                WebkitTextFillColor: 'transparent',
              }}
            >
              BioGemini
            </Typography>
          </Toolbar>
        </StyledAppBar>

        {/* Drawer */}
        <StyledDrawer
          variant={isMobile ? 'temporary' : 'permanent'}
          open={isMobile ? isSidebarOpen : true}
          onClose={() => setIsSidebarOpen(false)}
          ModalProps={{ keepMounted: true }}
          sx={{
            width: drawerWidth,
            flexShrink: 0,
            '& .MuiDrawer-paper': {
              width: drawerWidth,
              boxSizing: 'border-box',
              transition: 'width 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
              overflowX: 'hidden',
            },
          }}
        >
          <Toolbar />
          <DrawerContent />
        </StyledDrawer>

        {/* Main Content */}
        <Box
          component="main"
          sx={{
            flexGrow: 1,
            transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
            bgcolor: 'background.default',
            minHeight: '100vh',
            overflowX: 'hidden',
            maxWidth: '100%',
          }}
        >
          <Toolbar />
          <ScrollableContainer sx={{
            p: { xs: 2, md: 3 },
            height: 'calc(100vh - 64px)',
            overflowY: 'auto',
            overflowX: 'hidden',
          }}>
            <Outlet />
          </ScrollableContainer>
        </Box>

        {/* Chatbot */}
        <JarvisBot />

        {/* Toast notifications */}
        <Toaster position="top-right" />
      </Box>
    </ThemeProvider>
  );
}

export default DashboardPage;