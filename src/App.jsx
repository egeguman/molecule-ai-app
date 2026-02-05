import React, { useState, useEffect } from 'react';
import './index.css';
import Login from './components/Login';
import Register from './components/Register';

const App = () => {
  const [theme, setTheme] = useState(() => localStorage.getItem('theme') || 'light');

  useEffect(() => {
    localStorage.setItem('theme', theme);
  }, [theme]);

  const styles = getStyles(theme);
  const [showIntro, setShowIntro] = useState(true);
  const [authScreen, setAuthScreen] = useState(null); // 'register', 'login', or null (dashboard)
  const [isLoggedIn, setIsLoggedIn] = useState(false);
  const [smilesInput, setSmilesInput] = useState('');
  const [promptInput, setPromptInput] = useState('');
  const [currentMolecules, setCurrentMolecules] = useState([]);
  const [analysisRequests, setAnalysisRequests] = useState([]);
  const [activeTab, setActiveTab] = useState('input');
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [identifierType, setIdentifierType] = useState('smiles'); // 'smiles', 'chembl', or 'name'

  const handleAddMolecule = () => {
    if (smilesInput.trim()) {
      const newMolecule = {
        id: Date.now() + Math.random(),
        smiles: smilesInput.trim(),
        name: `Molecule ${currentMolecules.length + 1}`,
        identifierType: identifierType,
        identifierValue: smilesInput.trim(),
      };
      setCurrentMolecules(prev => [...prev, newMolecule]);
      setSmilesInput('');
    }
  };

  const handleRemoveCurrentMolecule = (id) => {
    setCurrentMolecules(prev => prev.filter(m => m.id !== id));
  };

  const handleSubmitRequest = () => {
    if (currentMolecules.length > 0 && promptInput.trim()) {
      const newRequest = {
        id: Date.now() + Math.random(),
        molecules: [...currentMolecules],
        prompt: promptInput.trim(),
        createdAt: new Date().toLocaleString(),
        status: 'pending'
      };
      setAnalysisRequests(prev => [...prev, newRequest]);
      setCurrentMolecules([]);
      setPromptInput('');
    }
  };

  const handleDeleteRequest = (id) => {
    setAnalysisRequests(prev => prev.filter(r => r.id !== id));
  };

  const handleAnalyzeRequest = (request) => {
    setIsAnalyzing(true);
    setTimeout(() => {
      setAnalysisRequests(prev => prev.map(r =>
        r.id === request.id ? { ...r, status: 'completed' } : r
      ));
      setIsAnalyzing(false);
      setActiveTab('results');
    }, 2000);
  };

  const handleAnalyzeAll = () => {
    setIsAnalyzing(true);
    setTimeout(() => {
      setAnalysisRequests(prev => prev.map(r => ({ ...r, status: 'completed' })));
      setIsAnalyzing(false);
      setActiveTab('results');
    }, 2000);
  };

  const exampleSmiles = [
    { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
    { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' }
  ];

  const examplePrompts = [
    { label: 'Compare properties', prompt: 'Compare the drug properties of these molecules including solubility, bioavailability, and potential toxicity' },
    { label: 'Binding analysis', prompt: 'Analyze and compare the binding affinity of these molecules to common drug targets' },
    { label: 'ADMET comparison', prompt: 'Compare ADMET properties and drug-likeness scores between these molecules' },
    { label: 'Similarity analysis', prompt: 'Identify structural similarities and differences between these molecules and predict similar biological activity' }
  ];

  const canSubmit = currentMolecules.length > 0 && promptInput.trim();
  const pendingRequests = analysisRequests.filter(r => r.status === 'pending');
  // Handle getting started from intro page
  const handleGetStarted = () => {
    setShowIntro(false);
    setAuthScreen('register'); // Show register screen after intro
  };

  // Handle successful login
  const handleLogin = (data) => {
    console.log('Logged in', data);
    setIsLoggedIn(true);
    setAuthScreen(null); // Go to dashboard
  };

  // Handle successful registration
  const handleRegister = (data) => {
    console.log('Registered', data);
    setAuthScreen('login'); // Go to login after registration
  };

  // Handle switching between auth screens
  const handleSwitchToLogin = () => {
    setAuthScreen('login');
  };

  const handleSwitchToRegister = () => {
    setAuthScreen('register');
  };

  // Handle exit/logout from dashboard
  const handleExit = () => {
    setIsLoggedIn(false);
    setShowIntro(true);
    setAuthScreen(null);
  };

  // If showing intro page, render the landing page
  if (showIntro) {
    return (
      <div style={introStyles.container}>
        {/* Animated Background */}
        <div style={introStyles.backgroundGradient}></div>
        <div style={introStyles.backgroundOrbs}>
          <div style={introStyles.orb1}></div>
          <div style={introStyles.orb2}></div>
          <div style={introStyles.orb3}></div>
        </div>

        {/* Navigation */}
        <nav style={introStyles.nav}>
          <div style={introStyles.navLogo}>
            <div style={introStyles.navLogoIcon}>
              <svg width="28" height="28" viewBox="0 0 32 32" fill="none">
                <circle cx="16" cy="16" r="6" fill="#16a34a" />
                <circle cx="8" cy="10" r="4" fill="#16a34a" opacity="0.8" />
                <circle cx="24" cy="10" r="4" fill="#16a34a" opacity="0.8" />
                <circle cx="10" cy="24" r="3" fill="#16a34a" opacity="0.6" />
                <circle cx="22" cy="24" r="3" fill="#16a34a" opacity="0.6" />
                <line x1="12" y1="12" x2="14" y2="14" stroke="#16a34a" strokeWidth="2" />
                <line x1="20" y1="12" x2="18" y2="14" stroke="#16a34a" strokeWidth="2" />
              </svg>
            </div>
            <span style={introStyles.navLogoText}>MoleculeAI</span>
          </div>
          <button style={introStyles.navCta} onClick={handleGetStarted}>
            Launch App →
          </button>
        </nav>

        {/* Hero Section */}
        <section style={introStyles.hero}>
          <div style={introStyles.heroContent}>
            <div style={introStyles.heroBadge}>
              <span style={introStyles.heroBadgeIcon}>✨</span>
              <span>AI-Powered Drug Discovery</span>
            </div>
            <h1 style={introStyles.heroTitle}>
              Accelerate Your
              <span style={introStyles.heroTitleGradient}> Drug Discovery</span>
              <br />Pipeline with AI
            </h1>
            <p style={introStyles.heroSubtitle}>
              Analyze molecular structures, predict ADMET properties, and discover
              potential drug candidates using state-of-the-art machine learning models.
            </p>
            <div style={introStyles.heroButtons}>
              <button style={introStyles.heroPrimaryBtn} onClick={handleGetStarted}>
                Start Analyzing
                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" style={{ marginLeft: '8px' }}>
                  <path d="M5 12H19M19 12L12 5M19 12L12 19" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
                </svg>
              </button>
              <button style={introStyles.heroSecondaryBtn}>
                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" style={{ marginRight: '8px' }}>
                  <circle cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="2" />
                  <polygon points="10,8 16,12 10,16" fill="currentColor" />
                </svg>
                Watch Demo
              </button>
            </div>
          </div>

          {/* Molecule Visualization */}
          <div style={introStyles.heroVisual}>
            <div style={introStyles.moleculeContainer}>
              <div style={introStyles.moleculeCard}>
                <div style={introStyles.moleculeGlow}></div>
                <svg width="280" height="280" viewBox="0 0 280 280" style={introStyles.moleculeSvg}>
                  {/* Central complex molecule structure */}
                  <defs>
                    <linearGradient id="atomGradient" x1="0%" y1="0%" x2="100%" y2="100%">
                      <stop offset="0%" stopColor="#22c55e" />
                      <stop offset="100%" stopColor="#16a34a" />
                    </linearGradient>
                    <linearGradient id="bondGradient" x1="0%" y1="0%" x2="100%" y2="100%">
                      <stop offset="0%" stopColor="#16a34a" />
                      <stop offset="100%" stopColor="#15803d" />
                    </linearGradient>
                    <filter id="glow">
                      <feGaussianBlur stdDeviation="3" result="coloredBlur" />
                      <feMerge>
                        <feMergeNode in="coloredBlur" />
                        <feMergeNode in="SourceGraphic" />
                      </feMerge>
                    </filter>
                  </defs>

                  {/* Bonds */}
                  <line x1="140" y1="140" x2="90" y2="80" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="140" y1="140" x2="190" y2="80" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="140" y1="140" x2="70" y2="160" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="140" y1="140" x2="210" y2="160" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="140" y1="140" x2="100" y2="220" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="140" y1="140" x2="180" y2="220" stroke="url(#bondGradient)" strokeWidth="3" opacity="0.8" />
                  <line x1="90" y1="80" x2="50" y2="50" stroke="url(#bondGradient)" strokeWidth="2.5" opacity="0.6" />
                  <line x1="190" y1="80" x2="230" y2="50" stroke="url(#bondGradient)" strokeWidth="2.5" opacity="0.6" />
                  <line x1="70" y1="160" x2="30" y2="180" stroke="url(#bondGradient)" strokeWidth="2.5" opacity="0.6" />
                  <line x1="210" y1="160" x2="250" y2="180" stroke="url(#bondGradient)" strokeWidth="2.5" opacity="0.6" />

                  {/* Atoms */}
                  <circle cx="140" cy="140" r="20" fill="url(#atomGradient)" filter="url(#glow)" />
                  <circle cx="90" cy="80" r="14" fill="url(#atomGradient)" filter="url(#glow)" />
                  <circle cx="190" cy="80" r="14" fill="url(#atomGradient)" filter="url(#glow)" />
                  <circle cx="70" cy="160" r="12" fill="#3b82f6" filter="url(#glow)" />
                  <circle cx="210" cy="160" r="12" fill="#3b82f6" filter="url(#glow)" />
                  <circle cx="100" cy="220" r="10" fill="#ef4444" filter="url(#glow)" />
                  <circle cx="180" cy="220" r="10" fill="#ef4444" filter="url(#glow)" />
                  <circle cx="50" cy="50" r="8" fill="#16a34a" opacity="0.7" />
                  <circle cx="230" cy="50" r="8" fill="#16a34a" opacity="0.7" />
                  <circle cx="30" cy="180" r="7" fill="#3b82f6" opacity="0.6" />
                  <circle cx="250" cy="180" r="7" fill="#3b82f6" opacity="0.6" />

                  {/* Atom labels */}
                  <text x="140" y="145" textAnchor="middle" fill="white" fontSize="14" fontWeight="bold">C</text>
                  <text x="90" y="84" textAnchor="middle" fill="white" fontSize="11" fontWeight="bold">C</text>
                  <text x="190" y="84" textAnchor="middle" fill="white" fontSize="11" fontWeight="bold">C</text>
                  <text x="70" y="164" textAnchor="middle" fill="white" fontSize="10" fontWeight="bold">N</text>
                  <text x="210" y="164" textAnchor="middle" fill="white" fontSize="10" fontWeight="bold">N</text>
                  <text x="100" y="224" textAnchor="middle" fill="white" fontSize="9" fontWeight="bold">O</text>
                  <text x="180" y="224" textAnchor="middle" fill="white" fontSize="9" fontWeight="bold">O</text>
                </svg>
                <div style={introStyles.moleculeLabel}>
                  <span style={introStyles.moleculeName}>Caffeine</span>
                  <code style={introStyles.moleculeSmiles}>C8H10N4O2</code>
                </div>
              </div>
            </div>

            {/* Floating stats cards */}
            <div style={introStyles.floatingCard1}>
              <div style={introStyles.floatingCardIcon}>📊</div>
              <div>
                <div style={introStyles.floatingCardValue}>98.5%</div>
                <div style={introStyles.floatingCardLabel}>Accuracy</div>
              </div>
            </div>
            <div style={introStyles.floatingCard2}>
              <div style={introStyles.floatingCardIcon}>⚡</div>
              <div>
                <div style={introStyles.floatingCardValue}>&lt;2s</div>
                <div style={introStyles.floatingCardLabel}>Analysis Time</div>
              </div>
            </div>
            <div style={introStyles.floatingCard3}>
              <div style={introStyles.floatingCardIcon}>🧬</div>
              <div>
                <div style={introStyles.floatingCardValue}>1M+</div>
                <div style={introStyles.floatingCardLabel}>Molecules</div>
              </div>
            </div>
          </div>
        </section>

        {/* Features Section */}
        <section style={introStyles.features}>
          <div style={introStyles.featuresHeader}>
            <h2 style={introStyles.featuresTitle}>Powerful Features for Drug Discovery</h2>
            <p style={introStyles.featuresSubtitle}>
              Everything you need to accelerate your research workflow
            </p>
          </div>

          <div style={introStyles.featuresGrid}>
            <div style={introStyles.featureCard}>
              <div style={introStyles.featureIcon}>
                <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                  <rect x="4" y="4" width="24" height="24" rx="6" stroke="#16a34a" strokeWidth="2" />
                  <path d="M10 16L14 20L22 12" stroke="#16a34a" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round" />
                </svg>
              </div>
              <h3 style={introStyles.featureTitle}>ADMET Prediction</h3>
              <p style={introStyles.featureDesc}>
                Predict absorption, distribution, metabolism, excretion, and toxicity properties with high accuracy.
              </p>
            </div>

            <div style={introStyles.featureCard}>
              <div style={introStyles.featureIcon}>
                <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                  <circle cx="16" cy="16" r="10" stroke="#3b82f6" strokeWidth="2" />
                  <circle cx="16" cy="16" r="4" fill="#3b82f6" />
                  <path d="M16 6V10M16 22V26M6 16H10M22 16H26" stroke="#3b82f6" strokeWidth="2" strokeLinecap="round" />
                </svg>
              </div>
              <h3 style={introStyles.featureTitle}>Binding Affinity</h3>
              <p style={introStyles.featureDesc}>
                Analyze molecular interactions and predict binding affinity to target proteins.
              </p>
            </div>

            <div style={introStyles.featureCard}>
              <div style={introStyles.featureIcon}>
                <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                  <path d="M4 24L10 18L16 22L22 14L28 8" stroke="#f59e0b" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round" />
                  <circle cx="4" cy="24" r="2.5" fill="#f59e0b" />
                  <circle cx="16" cy="22" r="2.5" fill="#f59e0b" />
                  <circle cx="28" cy="8" r="2.5" fill="#f59e0b" />
                </svg>
              </div>
              <h3 style={introStyles.featureTitle}>Property Analysis</h3>
              <p style={introStyles.featureDesc}>
                Compare molecular properties including solubility, bioavailability, and drug-likeness scores.
              </p>
            </div>

            <div style={introStyles.featureCard}>
              <div style={introStyles.featureIcon}>
                <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                  <path d="M16 4L20 12H28L22 18L24 26L16 22L8 26L10 18L4 12H12L16 4Z" stroke="#ec4899" strokeWidth="2" strokeLinejoin="round" />
                </svg>
              </div>
              <h3 style={introStyles.featureTitle}>AI Insights</h3>
              <p style={introStyles.featureDesc}>
                Get intelligent recommendations and insights powered by advanced machine learning models.
              </p>
            </div>
          </div>
        </section>

        {/* CTA Section */}
        <section style={introStyles.cta}>
          <div style={introStyles.ctaCard}>
            <h2 style={introStyles.ctaTitle}>Ready to Transform Your Research?</h2>
            <p style={introStyles.ctaText}>
              Start analyzing molecules and discovering potential drug candidates today.
            </p>
            <button style={introStyles.ctaButton} onClick={handleGetStarted}>
              Get Started Free
              <svg width="20" height="20" viewBox="0 0 24 24" fill="none" style={{ marginLeft: '8px' }}>
                <path d="M5 12H19M19 12L12 5M19 12L12 19" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
              </svg>
            </button>
          </div>
        </section>

        {/* Footer */}
        <footer style={introStyles.footer}>
          <div style={introStyles.footerContent}>
            <div style={introStyles.footerLogo}>
              <svg width="24" height="24" viewBox="0 0 32 32" fill="none">
                <circle cx="16" cy="16" r="6" fill="#16a34a" />
                <circle cx="8" cy="10" r="4" fill="#16a34a" opacity="0.6" />
                <circle cx="24" cy="10" r="4" fill="#16a34a" opacity="0.6" />
              </svg>
              <span style={introStyles.footerLogoText}>MoleculeAI</span>
            </div>
            <p style={introStyles.footerCopyright}>© 2025 MoleculeAI. Advancing drug discovery with AI.</p>
          </div>
        </footer>
      </div>
    );
  }

  // If showing auth screen (register or login)
  if (authScreen) {
    return (
      <div style={styles.authContainer}>
        {/* Background with gradient */}
        <div style={introStyles.backgroundGradient}></div>
        <div style={introStyles.backgroundOrbs}>
          <div style={introStyles.orb1}></div>
          <div style={introStyles.orb2}></div>
          <div style={introStyles.orb3}></div>
        </div>

        {/* Back to intro button */}
        <button
          style={styles.backToIntroButton}
          onClick={() => { setShowIntro(true); setAuthScreen(null); }}
        >
          ← Back to Home
        </button>

        {authScreen === 'register' ? (
          <Register
            onRegister={handleRegister}
            onSwitchToLogin={handleSwitchToLogin}
          />
        ) : (
          <Login
            onLogin={handleLogin}
            onSwitchToRegister={handleSwitchToRegister}
          />
        )}
      </div>
    );
  }

  return (
    <div style={styles.container}>
      {/* Header */}
      <header style={styles.header}>
        <div style={styles.logoSection}>
          <div style={styles.logoIcon}>
            <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
              <circle cx="16" cy="16" r="6" fill="#16a34a" />
              <circle cx="8" cy="10" r="4" fill="#16a34a" />
              <circle cx="24" cy="10" r="4" fill="#16a34a" />
              <circle cx="10" cy="24" r="3" fill="#16a34a" />
              <circle cx="22" cy="24" r="3" fill="#16a34a" />
              <line x1="12" y1="12" x2="14" y2="14" stroke="#16a34a" strokeWidth="2" />
              <line x1="20" y1="12" x2="18" y2="14" stroke="#16a34a" strokeWidth="2" />
              <line x1="12" y1="22" x2="14" y2="19" stroke="#16a34a" strokeWidth="2" />
              <line x1="20" y1="22" x2="18" y2="19" stroke="#16a34a" strokeWidth="2" />
            </svg>
          </div>
          <div>
            <h1 style={styles.logoTitle}>MoleculeAI</h1>
            <p style={styles.logoSubtitle}>Drug Discovery Platform</p>
          </div>
        </div>
        <nav style={styles.nav}>
          <button
            style={{ ...styles.navButton, ...(activeTab === 'input' ? styles.navButtonActive : {}) }}
            onClick={() => setActiveTab('input')}
          >
            New Analysis
          </button>
          <button
            style={{ ...styles.navButton, ...(activeTab === 'library' ? styles.navButtonActive : {}) }}
            onClick={() => setActiveTab('library')}
          >
            Requests ({analysisRequests.length})
          </button>
          <button
            style={{ ...styles.navButton, ...(activeTab === 'results' ? styles.navButtonActive : {}) }}
            onClick={() => setActiveTab('results')}
          >
            Results
          </button>
          <div style={styles.navDivider}></div>
          <button
            style={styles.exitButton}
            onClick={handleExit}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" style={{ marginRight: '6px' }}>
              <path d="M9 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h4"></path>
              <polyline points="16 17 21 12 16 7"></polyline>
              <line x1="21" y1="12" x2="9" y2="12"></line>
            </svg>
            Exit
          </button>
          <button
            style={styles.themeToggle}
            onClick={() => setTheme(t => t === 'light' ? 'dark' : 'light')}
            title={`Switch to ${theme === 'light' ? 'dark' : 'light'} mode`}
          >
            {theme === 'light' ? (
              <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>
              </svg>
            ) : (
              <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                <circle cx="12" cy="12" r="5"></circle>
                <line x1="12" y1="1" x2="12" y2="3"></line>
                <line x1="12" y1="21" x2="12" y2="23"></line>
                <line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line>
                <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line>
                <line x1="1" y1="12" x2="3" y2="12"></line>
                <line x1="21" y1="12" x2="23" y2="12"></line>
                <line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line>
                <line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line>
              </svg>
            )}
          </button>
        </nav>
      </header>

      <main style={styles.main}>
        {/* Input Tab */}
        {activeTab === 'input' && (
          <div style={styles.inputLayoutVertical}>
            {/* Main Input Card */}
            <section style={styles.mainCard}>
              <div style={styles.cardHeader}>
                <h2 style={styles.cardTitle}>Create Analysis Request</h2>
              </div>

              {/* Molecule Identifier Input Area */}
              <div style={styles.smilesSection}>
                <div style={styles.sectionHeader}>
                  <span style={styles.sectionTitle}>1. Add Molecules</span>
                </div>

                {/* Identifier Type Dropdown */}
                <div style={styles.identifierTypeSection}>
                  <label style={styles.dropdownLabel}>Identifier Type:</label>
                  <select
                    value={identifierType}
                    onChange={(e) => setIdentifierType(e.target.value)}
                    style={styles.dropdown}
                  >
                    <option value="smiles">SMILES</option>
                    <option value="chembl">ChEMBL ID</option>
                    <option value="name">Molecule Name</option>
                  </select>
                </div>

                {/* Dynamic Input Field */}
                <input
                  type="text"
                  value={smilesInput}
                  onChange={(e) => setSmilesInput(e.target.value)}
                  placeholder={
                    identifierType === 'smiles'
                      ? 'Enter SMILES code and press Enter (e.g., CC(=O)OC1=CC=CC=C1C(=O)O)'
                      : identifierType === 'chembl'
                        ? 'Enter ChEMBL ID and press Enter (e.g., CHEMBL25)'
                        : 'Enter Molecule Name and press Enter (e.g., Aspirin)'
                  }
                  style={styles.smilesInputFieldFull}
                  onKeyPress={(e) => e.key === 'Enter' && handleAddMolecule()}
                />

                <div style={styles.quickAdd}>
                  <span style={styles.quickAddLabel}>Quick add:</span>
                  {exampleSmiles.map((ex) => (
                    <button
                      key={ex.name}
                      style={styles.chip}
                      onClick={() => {
                        const newMolecule = {
                          id: Date.now() + Math.random(),
                          smiles: ex.smiles,
                          name: ex.name,
                          identifierType: 'smiles',
                          identifierValue: ex.smiles,
                        };
                        setCurrentMolecules([...currentMolecules, newMolecule]);
                      }}
                    >
                      + {ex.name}
                    </button>
                  ))}
                </div>

                {/* Current Molecules List */}
                {currentMolecules.length > 0 && (
                  <div style={styles.moleculesList}>
                    {currentMolecules.map((mol, index) => (
                      <div key={mol.id} style={styles.moleculeTag}>
                        <span style={styles.moleculeTagIndex}>{index + 1}</span>
                        <div style={styles.moleculeTagContent}>
                          <span style={styles.moleculeTagName}>{mol.name}</span>
                          <div style={styles.moleculeTagIdentifier}>
                            <span style={styles.moleculeTagType}>
                              {mol.identifierType === 'smiles' ? 'SMILES' :
                                mol.identifierType === 'chembl' ? 'ChEMBL ID' :
                                  'Name'}:
                            </span>
                            <code style={styles.moleculeTagSmiles}>
                              {mol.identifierValue || mol.smiles}
                            </code>
                          </div>
                        </div>
                        <button
                          style={styles.moleculeTagRemove}
                          onClick={() => handleRemoveCurrentMolecule(mol.id)}
                        >
                          ×
                        </button>
                      </div>
                    ))}
                  </div>
                )}
              </div>

              {/* Divider */}
              <div style={styles.divider}></div>

              {/* Prompt Input Area */}
              <div style={styles.promptSection}>
                <div style={styles.sectionHeader}>
                  <span style={styles.sectionTitle}>2. Analysis Prompt</span>
                  {promptInput.trim() && <span style={styles.promptReady}>Ready</span>}
                </div>

                <textarea
                  value={promptInput}
                  onChange={(e) => setPromptInput(e.target.value)}
                  placeholder="Describe what you want to analyze. For example: 'Compare the drug properties of these molecules including solubility, bioavailability, and potential toxicity...'"
                  style={styles.promptTextarea}
                  rows={4}
                />

                <div style={styles.promptSuggestions}>
                  <span style={styles.quickAddLabel}>Suggested:</span>
                  <div style={styles.suggestionChips}>
                    {examplePrompts.map((ex) => (
                      <button
                        key={ex.label}
                        style={styles.suggestionChip}
                        onClick={() => setPromptInput(ex.prompt)}
                      >
                        {ex.label}
                      </button>
                    ))}
                  </div>
                </div>
              </div>

              {/* Submit Button */}
              <div style={styles.submitSection}>
                <button
                  style={{ ...styles.submitButton, ...(!canSubmit ? styles.buttonDisabled : {}) }}
                  onClick={handleSubmitRequest}
                  disabled={!canSubmit}
                >
                  Add to Analysis Queue
                </button>
                {!canSubmit && (
                  <p style={styles.submitHint}>
                    {currentMolecules.length === 0 && !promptInput.trim()
                      ? 'Add at least one molecule and enter a prompt'
                      : currentMolecules.length === 0
                        ? 'Add at least one molecule'
                        : 'Enter an analysis prompt'}
                  </p>
                )}
              </div>
            </section>

            {/* Queue Section Below */}
            <section style={styles.queueSection}>
              <div style={styles.queueSectionHeader}>
                <div style={styles.queueTitleRow}>
                  <h3 style={styles.queueSectionTitle}>Analysis Queue</h3>
                  <span style={styles.queueCount}>{pendingRequests.length}</span>
                </div>
                {pendingRequests.length > 0 && (
                  <div style={styles.queueActions}>
                    <button
                      style={styles.runAllButtonInline}
                      onClick={handleAnalyzeAll}
                      disabled={isAnalyzing}
                    >
                      {isAnalyzing ? 'Analyzing...' : 'Run All Analyses'}
                    </button>
                    <button
                      style={styles.clearAllButton}
                      onClick={() => setAnalysisRequests(prev => prev.filter(r => r.status === 'completed'))}
                    >
                      Clear Queue
                    </button>
                  </div>
                )}
              </div>

              {pendingRequests.length === 0 ? (
                <div style={styles.queueEmptyHorizontal}>
                  <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
                    <rect x="4" y="8" width="24" height="16" rx="3" stroke="#d1d5db" strokeWidth="1.5" strokeDasharray="3 3" />
                    <path d="M10 14H22M10 18H18" stroke="#d1d5db" strokeWidth="1.5" strokeLinecap="round" />
                  </svg>
                  <div>
                    <p style={styles.queueEmptyTextInline}>No requests in queue</p>
                    <p style={styles.queueEmptyHintInline}>Create an analysis request above to get started</p>
                  </div>
                </div>
              ) : (
                <div style={styles.queueGrid}>
                  {pendingRequests.map((req, index) => (
                    <div key={req.id} style={styles.queueCard}>
                      <div style={styles.queueCardHeader}>
                        <div style={styles.queueCardTitleRow}>
                          <span style={styles.queueCardNumber}>#{index + 1}</span>
                          <span style={styles.queueCardMolCount}>{req.molecules.length} molecules</span>
                        </div>
                        <button
                          style={styles.queueCardRemove}
                          onClick={() => handleDeleteRequest(req.id)}
                        >
                          ×
                        </button>
                      </div>
                      <div style={styles.queueCardMolecules}>
                        {req.molecules.map(mol => (
                          <span key={mol.id} style={styles.queueCardMolChip}>{mol.name}</span>
                        ))}
                      </div>
                      <p style={styles.queueCardPrompt}>
                        {req.prompt.length > 80 ? req.prompt.substring(0, 80) + '...' : req.prompt}
                      </p>
                      <button
                        style={styles.queueCardRunButton}
                        onClick={() => handleAnalyzeRequest(req)}
                        disabled={isAnalyzing}
                      >
                        Run Analysis
                      </button>
                    </div>
                  ))}
                </div>
              )}
            </section>
          </div>
        )}

        {/* Library Tab */}
        {activeTab === 'library' && (
          <div style={styles.libraryContainer}>
            <section style={styles.card}>
              <div style={styles.cardHeader}>
                <h2 style={styles.cardTitle}>Analysis Requests</h2>
                <span style={styles.countBadge}>{analysisRequests.length} requests</span>
              </div>

              {analysisRequests.length === 0 ? (
                <div style={styles.emptyState}>
                  <div style={styles.emptyIcon}>
                    <svg width="48" height="48" viewBox="0 0 48 48" fill="none">
                      <rect x="8" y="12" width="32" height="24" rx="4" stroke="#d1d5db" strokeWidth="2" strokeDasharray="4 4" />
                      <path d="M16 22H32M16 28H28" stroke="#d1d5db" strokeWidth="2" strokeLinecap="round" />
                    </svg>
                  </div>
                  <p style={styles.emptyText}>No analysis requests</p>
                  <button
                    style={styles.primaryButton}
                    onClick={() => setActiveTab('input')}
                  >
                    Create Request
                  </button>
                </div>
              ) : (
                <div style={styles.requestsList}>
                  {analysisRequests.map((request, index) => (
                    <div key={request.id} style={styles.requestCard}>
                      <div style={styles.requestHeader}>
                        <div style={styles.requestTitleRow}>
                          <span style={styles.requestNumber}>#{index + 1}</span>
                          <span style={styles.requestTime}>{request.createdAt}</span>
                        </div>
                        <span style={{
                          ...styles.statusBadge,
                          ...(request.status === 'completed' ? styles.statusCompleted : styles.statusPending)
                        }}>
                          {request.status}
                        </span>
                      </div>

                      <div style={styles.requestContent}>
                        <div style={styles.requestMolecules}>
                          <span style={styles.requestLabel}>Molecules ({request.molecules.length})</span>
                          <div style={styles.requestMoleculesList}>
                            {request.molecules.map((mol) => (
                              <div key={mol.id} style={styles.requestMoleculeItem}>
                                <span style={styles.requestMolName}>{mol.name}</span>
                                <code style={styles.requestMolSmiles}>{mol.smiles}</code>
                              </div>
                            ))}
                          </div>
                        </div>

                        <div style={styles.requestPrompt}>
                          <span style={styles.requestLabel}>Prompt</span>
                          <p style={styles.requestPromptText}>{request.prompt}</p>
                        </div>
                      </div>

                      <div style={styles.requestActions}>
                        <button
                          style={styles.actionButton}
                          onClick={() => request.status === 'completed' ? setActiveTab('results') : handleAnalyzeRequest(request)}
                          disabled={isAnalyzing && request.status !== 'completed'}
                        >
                          {request.status === 'completed' ? 'View Results' : 'Run Analysis'}
                        </button>
                        <button
                          style={{ ...styles.actionButton, ...styles.dangerButton }}
                          onClick={() => handleDeleteRequest(request.id)}
                        >
                          Delete
                        </button>
                      </div>
                    </div>
                  ))}
                </div>
              )}

              {analysisRequests.length > 0 && (
                <div style={styles.bulkActions}>
                  <button
                    style={styles.primaryButton}
                    onClick={handleAnalyzeAll}
                    disabled={isAnalyzing}
                  >
                    {isAnalyzing ? 'Analyzing...' : 'Run All Analyses'}
                  </button>
                  <button
                    style={styles.secondaryButton}
                    onClick={() => setAnalysisRequests([])}
                  >
                    Clear All
                  </button>
                </div>
              )}
            </section>
          </div>
        )}

        {/* Results Tab */}
        {activeTab === 'results' && (
          <div style={styles.resultsContainer}>
            <section style={styles.card}>
              <div style={styles.cardHeader}>
                <h2 style={styles.cardTitle}>Analysis Results</h2>
                <span style={styles.badge}>Preview</span>
              </div>

              {analysisRequests.filter(r => r.status === 'completed').length === 0 ? (
                <div style={styles.resultsPlaceholder}>
                  <div style={styles.resultsIcon}>
                    <svg width="64" height="64" viewBox="0 0 64 64" fill="none">
                      <rect x="8" y="8" width="48" height="48" rx="8" stroke="#16a34a" strokeWidth="2" />
                      <path d="M20 40L28 32L36 38L44 24" stroke="#16a34a" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
                      <circle cx="20" cy="40" r="3" fill="#16a34a" />
                      <circle cx="28" cy="32" r="3" fill="#16a34a" />
                      <circle cx="36" cy="38" r="3" fill="#16a34a" />
                      <circle cx="44" cy="24" r="3" fill="#16a34a" />
                    </svg>
                  </div>
                  <h3 style={styles.resultsTitle}>No completed analyses</h3>
                  <p style={styles.resultsText}>
                    Create analysis requests and run them to see results here.
                  </p>
                  <button
                    style={styles.primaryButton}
                    onClick={() => setActiveTab('input')}
                  >
                    Create Request
                  </button>
                </div>
              ) : (
                <div style={styles.completedResults}>
                  {analysisRequests.filter(r => r.status === 'completed').map((request, index) => (
                    <div key={request.id} style={styles.resultCard}>
                      <div style={styles.resultHeader}>
                        <span style={styles.resultTitle}>Analysis #{index + 1}</span>
                        <span style={styles.statusCompletedBadge}>Completed</span>
                      </div>

                      {/* Request Summary - Shows molecules + prompt bundled together */}
                      <div style={styles.requestSummary}>
                        <div style={styles.summaryHeader}>
                          <svg width="16" height="16" viewBox="0 0 16 16" fill="none">
                            <path d="M4 6L8 10L12 6" stroke="#6b7280" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
                          </svg>
                          <span style={styles.summaryTitle}>Request Details</span>
                        </div>

                        <div style={styles.summaryContent}>
                          <div style={styles.summaryMolecules}>
                            <span style={styles.summaryLabel}>Molecules:</span>
                            <div style={styles.summaryMolList}>
                              {request.molecules.map(mol => (
                                <span key={mol.id} style={styles.summaryMolChip}>{mol.name}</span>
                              ))}
                            </div>
                          </div>
                          <div style={styles.summaryPromptSection}>
                            <span style={styles.summaryLabel}>Prompt:</span>
                            <p style={styles.summaryPromptText}>{request.prompt}</p>
                          </div>
                        </div>
                      </div>

                      {/* Results Section */}
                      <div style={styles.resultMetrics}>
                        <div style={styles.metricCard}>
                          <span style={styles.metricLabel}>Drug-likeness</span>
                          <span style={styles.metricValue}>0.85</span>
                          <div style={styles.progressBar}>
                            <div style={{ ...styles.progressFill, width: '85%' }}></div>
                          </div>
                        </div>
                        <div style={styles.metricCard}>
                          <span style={styles.metricLabel}>Toxicity Risk</span>
                          <span style={styles.metricValue}>Low</span>
                          <div style={styles.progressBar}>
                            <div style={{ ...styles.progressFill, width: '20%', backgroundColor: '#16a34a' }}></div>
                          </div>
                        </div>
                        <div style={styles.metricCard}>
                          <span style={styles.metricLabel}>Solubility</span>
                          <span style={styles.metricValue}>High</span>
                          <div style={styles.progressBar}>
                            <div style={{ ...styles.progressFill, width: '78%' }}></div>
                          </div>
                        </div>
                      </div>

                      <div style={styles.resultInsight}>
                        <span style={styles.insightLabel}>AI Insight (Preview)</span>
                        <p style={styles.insightText}>
                          Based on the molecular structures provided, the analysis indicates favorable drug-like properties.
                          Connect to backend for detailed comparison results.
                        </p>
                      </div>
                    </div>
                  ))}
                </div>
              )}
            </section>
          </div>
        )}

      </main>

      {/* Footer */}
      <footer style={styles.footer}>
        <p style={styles.footerText}>MoleculeAI Drug Discovery Platform • Frontend Preview</p>
      </footer>

    </div>
  );
};

const themeColors = {
  light: {
    bg: '#f8fafc',
    card: '#ffffff',
    text: '#1f2937',
    textSec: '#6b7280',
    border: '#e5e7eb',
    accent: '#16a34a',
    accentLight: '#f0fdf4',
    inputBg: '#fafafa',
    divider: '#f3f4f6',
    navText: '#4b5563',
    navBorder: '#e5e7eb',
    logoBg: '#f0fdf4',
    logoTitle: '#111827',
    chipBg: '#f3f4f6',
    chipBorder: '#bbf7d0',
    buttonDisabled: '#d1d5db',
    progressBg: '#e5e7eb',
    secondaryBg: '#f9fafb',
    footerText: '#9ca3af',
  },
  dark: {
    bg: '#0a0a0f',
    card: 'rgba(255, 255, 255, 0.05)',
    text: '#ffffff',
    textSec: 'rgba(255, 255, 255, 0.7)',
    border: 'rgba(255, 255, 255, 0.1)',
    accent: '#22c55e',
    accentLight: 'rgba(22, 163, 74, 0.15)',
    inputBg: 'rgba(255, 255, 255, 0.05)',
    divider: 'rgba(255, 255, 255, 0.1)',
    navText: '#e5e7eb',
    navBorder: 'rgba(255, 255, 255, 0.15)',
    logoBg: 'rgba(22, 163, 74, 0.15)',
    logoTitle: '#ffffff',
    chipBg: 'rgba(255, 255, 255, 0.1)',
    chipBorder: 'rgba(22, 163, 74, 0.5)',
    buttonDisabled: 'rgba(255, 255, 255, 0.1)',
    progressBg: 'rgba(255, 255, 255, 0.1)',
    secondaryBg: 'rgba(255, 255, 255, 0.05)',
    footerText: 'rgba(255, 255, 255, 0.4)',
  }
};

const getStyles = (theme) => {
  const colors = themeColors[theme];
  return {
    container: {
      minHeight: '100vh',
      backgroundColor: colors.bg,
      fontFamily: "'DM Sans', sans-serif",
      color: colors.text,
      display: 'flex',
      flexDirection: 'column',
    },
    header: {
      backgroundColor: colors.card,
      borderBottom: `1px solid ${colors.border}`,
      padding: '16px 32px',
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      position: 'sticky',
      top: 0,
      zIndex: 100,
    },
    logoSection: {
      display: 'flex',
      alignItems: 'center',
      gap: '12px',
    },
    logoIcon: {
      width: '48px',
      height: '48px',
      backgroundColor: colors.accentLight,
      borderRadius: '12px',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
    },
    logoTitle: {
      fontSize: '20px',
      fontWeight: '700',
      color: colors.text,
      margin: 0,
      letterSpacing: '-0.02em',
    },
    logoSubtitle: {
      fontSize: '13px',
      color: colors.textSec,
      margin: 0,
    },
    nav: {
      display: 'flex',
      gap: '8px',
    },
    navButton: {
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      backgroundColor: 'transparent',
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      color: colors.navText,
      fontFamily: "'DM Sans', sans-serif",
    },
    navButtonActive: {
      backgroundColor: colors.accent,
      borderColor: colors.accent,
      color: '#ffffff',
    },
    navDivider: {
      width: '1px',
      backgroundColor: colors.progressBg,
      height: '24px',
      margin: '0 8px',
      alignSelf: 'center',
    },
    authContainer: {
      minHeight: '100vh',
      display: 'flex',
      flexDirection: 'column',
      alignItems: 'center',
      justifyContent: 'center',
      position: 'relative',
      overflow: 'hidden',
      backgroundColor: '#0a0a0a',
    },
    backToIntroButton: {
      position: 'absolute',
      top: '24px',
      left: '24px',
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      backgroundColor: 'rgba(255, 255, 255, 0.1)',
      border: '1px solid rgba(255, 255, 255, 0.2)',
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.2s ease',
      color: '#ffffff',
      fontFamily: "'DM Sans', sans-serif",
      zIndex: 10,
    },
    exitButton: {
      display: 'flex',
      alignItems: 'center',
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      backgroundColor: '#fee2e2',
      border: '1px solid #fecaca',
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      color: '#dc2626',
      fontFamily: "'DM Sans', sans-serif",
    },
    main: {
      flex: 1,
      padding: '32px',
      maxWidth: '1400px',
      margin: '0 auto',
      width: '100%',
    },
    inputLayoutVertical: {
      display: 'flex',
      flexDirection: 'column',
      gap: '24px',
      maxWidth: '900px',
      margin: '0 auto',
    },
    mainCard: {
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '16px',
      padding: '28px',
    },
    card: {
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '16px',
      padding: '24px',
    },
    cardHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginBottom: '24px',
    },
    cardTitle: {
      fontSize: '18px',
      fontWeight: '600',
      color: colors.text,
      margin: 0,
    },
    badge: {
      fontSize: '12px',
      fontWeight: '500',
      color: colors.accent,
      backgroundColor: colors.accentLight,
      padding: '4px 10px',
      borderRadius: '20px',
    },
    countBadge: {
      fontSize: '13px',
      fontWeight: '500',
      color: colors.textSec,
      backgroundColor: colors.chipBg,
      padding: '4px 12px',
      borderRadius: '20px',
    },
    smilesSection: {
      marginBottom: '0',
    },
    sectionHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginBottom: '12px',
    },
    sectionTitle: {
      fontSize: '14px',
      fontWeight: '600',
      color: colors.text,
    },
    sectionCount: {
      fontSize: '12px',
      fontWeight: '500',
      color: colors.accent,
      backgroundColor: colors.accentLight,
      padding: '2px 8px',
      borderRadius: '10px',
    },
    identifierTypeSection: {
      display: 'flex',
      alignItems: 'center',
      gap: '12px',
      marginBottom: '12px',
    },
    dropdownLabel: {
      fontSize: '13px',
      fontWeight: '500',
      color: colors.text,
    },
    dropdown: {
      padding: '8px 12px',
      fontSize: '14px',
      fontFamily: "'DM Sans', sans-serif",
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      backgroundColor: colors.card,
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      color: colors.text,
      fontWeight: '500',
    },
    smilesInputFieldFull: {
      width: '100%',
      padding: '12px 16px',
      fontSize: '14px',
      fontFamily: "'JetBrains Mono', monospace",
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      backgroundColor: colors.inputBg,
      transition: 'all 0.15s ease',
      marginBottom: '12px',
    },
    quickAdd: {
      display: 'flex',
      alignItems: 'center',
      gap: '8px',
      flexWrap: 'wrap',
    },
    quickAddLabel: {
      fontSize: '12px',
      color: colors.textSec,
    },
    chip: {
      padding: '5px 12px',
      fontSize: '12px',
      fontWeight: '500',
      backgroundColor: colors.accentLight,
      border: `1px solid ${colors.chipBorder}`,
      borderRadius: '16px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      color: colors.accent,
      fontFamily: "'DM Sans', sans-serif",
    },
    moleculesList: {
      marginTop: '16px',
      display: 'flex',
      flexDirection: 'column',
      gap: '8px',
    },
    moleculeTag: {
      display: 'flex',
      alignItems: 'center',
      gap: '12px',
      padding: '12px 14px',
      backgroundColor: colors.inputBg,
      borderRadius: '8px',
      border: `1px solid ${colors.divider}`,
    },
    moleculeTagIndex: {
      width: '24px',
      height: '24px',
      backgroundColor: colors.accent,
      color: '#ffffff',
      borderRadius: '6px',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      fontSize: '12px',
      fontWeight: '600',
      flexShrink: 0,
    },
    moleculeTagContent: {
      flex: 1,
      minWidth: 0,
    },
    moleculeTagName: {
      fontSize: '13px',
      fontWeight: '600',
      color: colors.text,
      display: 'block',
      marginBottom: '2px',
    },
    moleculeTagSmiles: {
      fontSize: '11px',
      fontFamily: "'JetBrains Mono', monospace",
      color: colors.textSec,
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
    },
    moleculeTagIdentifier: {
      display: 'flex',
      gap: '6px',
      alignItems: 'center',
    },
    moleculeTagType: {
      fontSize: '10px',
      fontWeight: '600',
      color: colors.accent,
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
    },
    moleculeTagRemove: {
      width: '24px',
      height: '24px',
      backgroundColor: 'transparent',
      border: `1px solid ${colors.border}`,
      borderRadius: '6px',
      cursor: 'pointer',
      fontSize: '16px',
      color: colors.textSec,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      transition: 'all 0.15s ease',
      flexShrink: 0,
      fontFamily: "'DM Sans', sans-serif",
    },
    divider: {
      height: '1px',
      backgroundColor: colors.progressBg,
      margin: '24px 0',
    },
    promptSection: {
      marginBottom: '24px',
    },
    promptReady: {
      fontSize: '11px',
      fontWeight: '600',
      color: '#ffffff',
      backgroundColor: colors.accent,
      padding: '3px 10px',
      borderRadius: '10px',
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
    },
    promptTextarea: {
      width: '100%',
      padding: '14px 16px',
      fontSize: '14px',
      fontFamily: "'DM Sans', sans-serif",
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      backgroundColor: colors.inputBg,
      resize: 'vertical',
      minHeight: '100px',
      transition: 'all 0.15s ease',
      lineHeight: '1.6',
    },
    promptSuggestions: {
      marginTop: '12px',
      display: 'flex',
      alignItems: 'flex-start',
      gap: '8px',
      flexWrap: 'wrap',
    },
    suggestionChips: {
      display: 'flex',
      gap: '6px',
      flexWrap: 'wrap',
    },
    suggestionChip: {
      padding: '5px 12px',
      fontSize: '12px',
      fontWeight: '500',
      backgroundColor: colors.chipBg,
      border: `1px solid ${colors.border}`,
      borderRadius: '16px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      color: colors.navText,
      fontFamily: "'DM Sans', sans-serif",
    },
    submitSection: {
      paddingTop: '20px',
      borderTop: '1px solid #f3f4f6',
    },
    submitButton: {
      width: '100%',
      padding: '14px 24px',
      fontSize: '15px',
      fontWeight: '600',
      backgroundColor: colors.accent,
      color: '#ffffff',
      border: 'none',
      borderRadius: '10px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    buttonDisabled: {
      backgroundColor: colors.buttonDisabled,
      cursor: 'not-allowed',
    },
    submitHint: {
      fontSize: '12px',
      color: colors.textSec,
      textAlign: 'center',
      marginTop: '10px',
      margin: '10px 0 0 0',
    },
    queueSection: {
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '16px',
      padding: '24px',
    },
    queueSectionHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginBottom: '20px',
      paddingBottom: '16px',
      borderBottom: `1px solid ${colors.divider}`,
    },
    queueTitleRow: {
      display: 'flex',
      alignItems: 'center',
      gap: '10px',
    },
    queueSectionTitle: {
      fontSize: '16px',
      fontWeight: '600',
      color: colors.text,
      margin: 0,
    },
    queueCount: {
      fontSize: '12px',
      fontWeight: '600',
      color: colors.accent,
      backgroundColor: colors.accentLight,
      padding: '2px 8px',
      borderRadius: '10px',
    },
    queueActions: {
      display: 'flex',
      gap: '10px',
    },
    runAllButtonInline: {
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '600',
      backgroundColor: colors.accent,
      color: '#ffffff',
      border: 'none',
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    clearAllButton: {
      padding: '10px 20px',
      fontSize: '14px',
      fontWeight: '500',
      backgroundColor: colors.card,
      color: colors.textSec,
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    queueEmptyHorizontal: {
      display: 'flex',
      alignItems: 'center',
      gap: '16px',
      padding: '24px',
      backgroundColor: colors.inputBg,
      borderRadius: '10px',
      border: '1px dashed #e5e7eb',
    },
    queueEmptyTextInline: {
      fontSize: '14px',
      fontWeight: '500',
      color: colors.textSec,
      margin: '0 0 2px 0',
    },
    queueEmptyHintInline: {
      fontSize: '12px',
      color: colors.textSec,
      margin: 0,
    },
    queueGrid: {
      display: 'grid',
      gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))',
      gap: '16px',
    },
    queueCard: {
      backgroundColor: colors.inputBg,
      border: `1px solid ${colors.divider}`,
      borderRadius: '10px',
      padding: '16px',
    },
    queueCardHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'flex-start',
      marginBottom: '10px',
    },
    queueCardTitleRow: {
      display: 'flex',
      alignItems: 'center',
      gap: '8px',
    },
    queueCardNumber: {
      fontSize: '13px',
      fontWeight: '700',
      color: colors.text,
    },
    queueCardMolCount: {
      fontSize: '11px',
      fontWeight: '500',
      color: colors.accent,
      backgroundColor: colors.accentLight,
      padding: '2px 6px',
      borderRadius: '8px',
    },
    queueCardRemove: {
      width: '22px',
      height: '22px',
      backgroundColor: 'transparent',
      border: `1px solid ${colors.border}`,
      borderRadius: '4px',
      cursor: 'pointer',
      fontSize: '14px',
      color: colors.textSec,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    queueCardMolecules: {
      display: 'flex',
      gap: '4px',
      flexWrap: 'wrap',
      marginBottom: '8px',
    },
    queueCardMolChip: {
      padding: '3px 8px',
      fontSize: '11px',
      fontWeight: '500',
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '10px',
      color: colors.text,
    },
    queueCardPrompt: {
      fontSize: '12px',
      lineHeight: '1.4',
      color: colors.textSec,
      margin: '0 0 12px 0',
    },
    queueCardRunButton: {
      width: '100%',
      padding: '8px 12px',
      fontSize: '12px',
      fontWeight: '600',
      backgroundColor: colors.accent,
      color: '#ffffff',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    libraryContainer: {
      maxWidth: '1000px',
      margin: '0 auto',
    },
    emptyState: {
      textAlign: 'center',
      padding: '48px 24px',
    },
    emptyIcon: {
      marginBottom: '16px',
    },
    emptyText: {
      fontSize: '16px',
      fontWeight: '500',
      color: colors.textSec,
      margin: '0 0 16px 0',
    },
    primaryButton: {
      padding: '12px 24px',
      fontSize: '14px',
      fontWeight: '600',
      backgroundColor: colors.accent,
      color: '#ffffff',
      border: 'none',
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    secondaryButton: {
      padding: '12px 24px',
      fontSize: '14px',
      fontWeight: '500',
      backgroundColor: colors.card,
      color: colors.navText,
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    requestsList: {
      display: 'flex',
      flexDirection: 'column',
      gap: '16px',
    },
    requestCard: {
      backgroundColor: colors.inputBg,
      border: `1px solid ${colors.divider}`,
      borderRadius: '12px',
      padding: '20px',
    },
    requestHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginBottom: '16px',
      paddingBottom: '12px',
      borderBottom: `1px solid ${colors.border}`,
    },
    requestTitleRow: {
      display: 'flex',
      alignItems: 'center',
      gap: '12px',
    },
    requestNumber: {
      fontSize: '14px',
      fontWeight: '700',
      color: colors.text,
    },
    requestTime: {
      fontSize: '12px',
      color: colors.textSec,
    },
    statusBadge: {
      fontSize: '11px',
      fontWeight: '600',
      padding: '4px 10px',
      borderRadius: '12px',
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
    },
    statusPending: {
      backgroundColor: '#fef3c7',
      color: '#d97706',
    },
    statusCompleted: {
      backgroundColor: colors.accentLight,
      color: colors.accent,
    },
    requestContent: {
      display: 'grid',
      gridTemplateColumns: '1fr 1fr',
      gap: '20px',
      marginBottom: '16px',
    },
    requestMolecules: {},
    requestLabel: {
      fontSize: '12px',
      fontWeight: '600',
      color: colors.textSec,
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
      display: 'block',
      marginBottom: '8px',
    },
    requestMoleculesList: {
      display: 'flex',
      flexDirection: 'column',
      gap: '6px',
    },
    requestMoleculeItem: {
      padding: '8px 10px',
      backgroundColor: colors.card,
      borderRadius: '6px',
      border: `1px solid ${colors.border}`,
    },
    requestMolName: {
      fontSize: '13px',
      fontWeight: '600',
      color: colors.text,
      display: 'block',
      marginBottom: '2px',
    },
    requestMolSmiles: {
      fontSize: '10px',
      fontFamily: "'JetBrains Mono', monospace",
      color: colors.textSec,
      display: 'block',
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
    },
    requestPrompt: {},
    requestPromptText: {
      fontSize: '14px',
      lineHeight: '1.6',
      color: colors.text,
      margin: 0,
      padding: '10px 12px',
      backgroundColor: colors.card,
      borderRadius: '6px',
      border: `1px solid ${colors.border}`,
    },
    requestActions: {
      display: 'flex',
      gap: '8px',
      paddingTop: '12px',
      borderTop: `1px solid ${colors.border}`,
    },
    actionButton: {
      padding: '8px 16px',
      fontSize: '13px',
      fontWeight: '500',
      backgroundColor: colors.accent,
      color: '#ffffff',
      border: 'none',
      borderRadius: '6px',
      cursor: 'pointer',
      transition: 'all 0.15s ease',
      fontFamily: "'DM Sans', sans-serif",
    },
    dangerButton: {
      backgroundColor: colors.card,
      color: '#dc2626',
      border: '1px solid #fecaca',
    },
    bulkActions: {
      display: 'flex',
      gap: '12px',
      marginTop: '24px',
      paddingTop: '20px',
      borderTop: `1px solid ${colors.border}`,
    },
    resultsContainer: {
      maxWidth: '900px',
      margin: '0 auto',
    },
    resultsPlaceholder: {
      textAlign: 'center',
      padding: '48px 24px',
    },
    resultsIcon: {
      marginBottom: '24px',
    },
    resultsTitle: {
      fontSize: '20px',
      fontWeight: '600',
      color: colors.text,
      margin: '0 0 8px 0',
    },
    resultsText: {
      fontSize: '15px',
      color: colors.textSec,
      margin: '0 0 24px 0',
      maxWidth: '400px',
      marginLeft: 'auto',
      marginRight: 'auto',
      lineHeight: '1.6',
    },
    completedResults: {
      display: 'flex',
      flexDirection: 'column',
      gap: '20px',
    },
    resultCard: {
      backgroundColor: colors.inputBg,
      border: `1px solid ${colors.divider}`,
      borderRadius: '12px',
      padding: '24px',
    },
    resultHeader: {
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginBottom: '20px',
      paddingBottom: '12px',
      borderBottom: `1px solid ${colors.border}`,
    },
    resultTitle: {
      fontSize: '16px',
      fontWeight: '600',
      color: colors.text,
    },
    statusCompletedBadge: {
      fontSize: '11px',
      fontWeight: '600',
      padding: '4px 10px',
      borderRadius: '12px',
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
      backgroundColor: colors.accentLight,
      color: colors.accent,
    },
    requestSummary: {
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '10px',
      marginBottom: '20px',
      overflow: 'hidden',
    },
    summaryHeader: {
      display: 'flex',
      alignItems: 'center',
      gap: '8px',
      padding: '12px 16px',
      backgroundColor: colors.secondaryBg,
      borderBottom: `1px solid ${colors.border}`,
    },
    summaryTitle: {
      fontSize: '13px',
      fontWeight: '600',
      color: colors.text,
    },
    summaryContent: {
      padding: '16px',
    },
    summaryMolecules: {
      marginBottom: '12px',
    },
    summaryLabel: {
      fontSize: '11px',
      fontWeight: '600',
      color: colors.textSec,
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
      display: 'block',
      marginBottom: '6px',
    },
    summaryMolList: {
      display: 'flex',
      gap: '6px',
      flexWrap: 'wrap',
    },
    summaryMolChip: {
      padding: '4px 10px',
      fontSize: '12px',
      fontWeight: '500',
      backgroundColor: colors.accentLight,
      border: `1px solid ${colors.chipBorder}`,
      borderRadius: '14px',
      color: colors.accent,
    },
    summaryPromptSection: {},
    summaryPromptText: {
      fontSize: '13px',
      lineHeight: '1.5',
      color: colors.text,
      margin: 0,
      padding: '10px 12px',
      backgroundColor: colors.secondaryBg,
      borderRadius: '6px',
      border: `1px solid ${colors.divider}`,
    },
    resultMetrics: {
      display: 'grid',
      gridTemplateColumns: 'repeat(3, 1fr)',
      gap: '12px',
      marginBottom: '20px',
    },
    metricCard: {
      backgroundColor: colors.card,
      border: `1px solid ${colors.border}`,
      borderRadius: '8px',
      padding: '14px',
    },
    metricLabel: {
      fontSize: '12px',
      fontWeight: '500',
      color: colors.textSec,
      display: 'block',
      marginBottom: '6px',
    },
    metricValue: {
      fontSize: '18px',
      fontWeight: '700',
      color: colors.text,
      display: 'block',
      marginBottom: '10px',
    },
    progressBar: {
      height: '6px',
      backgroundColor: colors.progressBg,
      borderRadius: '3px',
      overflow: 'hidden',
    },
    progressFill: {
      height: '100%',
      backgroundColor: colors.accent,
      borderRadius: '3px',
      transition: 'width 0.3s ease',
    },
    resultInsight: {
      paddingTop: '16px',
      borderTop: `1px solid ${colors.border}`,
    },
    insightLabel: {
      fontSize: '11px',
      fontWeight: '600',
      color: colors.textSec,
      textTransform: 'uppercase',
      letterSpacing: '0.5px',
      display: 'block',
      marginBottom: '8px',
    },
    insightText: {
      fontSize: '14px',
      lineHeight: '1.6',
      color: colors.text,
      margin: 0,
      fontStyle: 'italic',
    },
    footer: {
      padding: '20px 32px',
      borderTop: `1px solid ${colors.border}`,
      textAlign: 'center',
      backgroundColor: colors.card,
    },
    footerText: {
      fontSize: '13px',
      color: colors.footerText,
      margin: 0,
    },
    themeToggle: {
      padding: '8px',
      fontSize: '14px',
      backgroundColor: 'transparent',
      border: `1px solid ${themeColors[theme].navBorder}`,
      borderRadius: '8px',
      cursor: 'pointer',
      color: themeColors[theme].navText,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      marginLeft: '8px',
      transition: 'all 0.15s ease',
    },
  };
};

// Intro/Landing Page Styles
const introStyles = {
  container: {
    minHeight: '100vh',
    backgroundColor: '#0a0a0f',
    fontFamily: "'DM Sans', sans-serif",
    color: '#ffffff',
    position: 'relative',
    overflow: 'hidden',
  },
  backgroundGradient: {
    position: 'absolute',
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    background: 'radial-gradient(ellipse at 50% 0%, rgba(22, 163, 74, 0.15) 0%, transparent 50%), radial-gradient(ellipse at 80% 50%, rgba(59, 130, 246, 0.1) 0%, transparent 40%), radial-gradient(ellipse at 20% 80%, rgba(236, 72, 153, 0.08) 0%, transparent 40%)',
    pointerEvents: 'none',
  },
  backgroundOrbs: {
    position: 'absolute',
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    pointerEvents: 'none',
    overflow: 'hidden',
  },
  orb1: {
    position: 'absolute',
    width: '600px',
    height: '600px',
    borderRadius: '50%',
    background: 'radial-gradient(circle, rgba(22, 163, 74, 0.12) 0%, transparent 70%)',
    top: '-200px',
    right: '-100px',
    animation: 'float 20s ease-in-out infinite',
  },
  orb2: {
    position: 'absolute',
    width: '400px',
    height: '400px',
    borderRadius: '50%',
    background: 'radial-gradient(circle, rgba(59, 130, 246, 0.1) 0%, transparent 70%)',
    bottom: '10%',
    left: '-100px',
    animation: 'float 25s ease-in-out infinite reverse',
  },
  orb3: {
    position: 'absolute',
    width: '300px',
    height: '300px',
    borderRadius: '50%',
    background: 'radial-gradient(circle, rgba(236, 72, 153, 0.08) 0%, transparent 70%)',
    top: '50%',
    right: '10%',
    animation: 'float 15s ease-in-out infinite',
  },
  nav: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    padding: '24px 48px',
    position: 'relative',
    zIndex: 10,
  },
  navLogo: {
    display: 'flex',
    alignItems: 'center',
    gap: '12px',
  },
  navLogoIcon: {
    width: '44px',
    height: '44px',
    backgroundColor: 'rgba(22, 163, 74, 0.15)',
    borderRadius: '12px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    border: '1px solid rgba(22, 163, 74, 0.3)',
  },
  navLogoText: {
    fontSize: '22px',
    fontWeight: '700',
    color: '#ffffff',
    letterSpacing: '-0.02em',
  },
  navCta: {
    padding: '12px 24px',
    fontSize: '14px',
    fontWeight: '600',
    backgroundColor: 'rgba(255, 255, 255, 0.1)',
    border: '1px solid rgba(255, 255, 255, 0.2)',
    borderRadius: '10px',
    color: '#ffffff',
    cursor: 'pointer',
    transition: 'all 0.2s ease',
    fontFamily: "'DM Sans', sans-serif",
    backdropFilter: 'blur(10px)',
  },
  hero: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    padding: '60px 48px 80px',
    position: 'relative',
    zIndex: 10,
    maxWidth: '1400px',
    margin: '0 auto',
    gap: '60px',
  },
  heroContent: {
    flex: 1,
    maxWidth: '600px',
  },
  heroBadge: {
    display: 'inline-flex',
    alignItems: 'center',
    gap: '8px',
    padding: '8px 16px',
    backgroundColor: 'rgba(22, 163, 74, 0.15)',
    border: '1px solid rgba(22, 163, 74, 0.3)',
    borderRadius: '30px',
    fontSize: '13px',
    fontWeight: '500',
    color: '#4ade80',
    marginBottom: '24px',
  },
  heroBadgeIcon: {
    fontSize: '14px',
  },
  heroTitle: {
    fontSize: '56px',
    fontWeight: '700',
    lineHeight: '1.1',
    color: '#ffffff',
    margin: '0 0 24px 0',
    letterSpacing: '-0.03em',
  },
  heroTitleGradient: {
    background: 'linear-gradient(135deg, #4ade80 0%, #22c55e 50%, #16a34a 100%)',
    WebkitBackgroundClip: 'text',
    WebkitTextFillColor: 'transparent',
    backgroundClip: 'text',
  },
  heroSubtitle: {
    fontSize: '18px',
    lineHeight: '1.7',
    color: 'rgba(255, 255, 255, 0.7)',
    margin: '0 0 36px 0',
  },
  heroButtons: {
    display: 'flex',
    gap: '16px',
    flexWrap: 'wrap',
  },
  heroPrimaryBtn: {
    display: 'inline-flex',
    alignItems: 'center',
    padding: '16px 32px',
    fontSize: '16px',
    fontWeight: '600',
    background: 'linear-gradient(135deg, #22c55e 0%, #16a34a 100%)',
    border: 'none',
    borderRadius: '12px',
    color: '#ffffff',
    cursor: 'pointer',
    transition: 'all 0.2s ease',
    fontFamily: "'DM Sans', sans-serif",
    boxShadow: '0 8px 32px rgba(22, 163, 74, 0.3)',
  },
  heroSecondaryBtn: {
    display: 'inline-flex',
    alignItems: 'center',
    padding: '16px 28px',
    fontSize: '16px',
    fontWeight: '600',
    backgroundColor: 'rgba(255, 255, 255, 0.08)',
    border: '1px solid rgba(255, 255, 255, 0.15)',
    borderRadius: '12px',
    color: '#ffffff',
    cursor: 'pointer',
    transition: 'all 0.2s ease',
    fontFamily: "'DM Sans', sans-serif",
    backdropFilter: 'blur(10px)',
  },
  heroVisual: {
    flex: 1,
    position: 'relative',
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    minHeight: '400px',
  },
  moleculeContainer: {
    position: 'relative',
  },
  moleculeCard: {
    position: 'relative',
    padding: '32px',
    backgroundColor: 'rgba(255, 255, 255, 0.05)',
    borderRadius: '24px',
    border: '1px solid rgba(255, 255, 255, 0.1)',
    backdropFilter: 'blur(20px)',
  },
  moleculeGlow: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
    width: '200px',
    height: '200px',
    background: 'radial-gradient(circle, rgba(22, 163, 74, 0.3) 0%, transparent 70%)',
    borderRadius: '50%',
    pointerEvents: 'none',
  },
  moleculeSvg: {
    position: 'relative',
    zIndex: 1,
  },
  moleculeLabel: {
    marginTop: '20px',
    textAlign: 'center',
  },
  moleculeName: {
    display: 'block',
    fontSize: '18px',
    fontWeight: '600',
    color: '#ffffff',
    marginBottom: '4px',
  },
  moleculeSmiles: {
    fontSize: '14px',
    color: 'rgba(255, 255, 255, 0.5)',
    fontFamily: "'JetBrains Mono', monospace",
  },
  floatingCard1: {
    position: 'absolute',
    top: '20px',
    left: '-40px',
    display: 'flex',
    alignItems: 'center',
    gap: '12px',
    padding: '14px 18px',
    backgroundColor: 'rgba(255, 255, 255, 0.1)',
    borderRadius: '14px',
    border: '1px solid rgba(255, 255, 255, 0.15)',
    backdropFilter: 'blur(20px)',
    animation: 'floatCard 6s ease-in-out infinite',
  },
  floatingCard2: {
    position: 'absolute',
    top: '50%',
    right: '-60px',
    transform: 'translateY(-50%)',
    display: 'flex',
    alignItems: 'center',
    gap: '12px',
    padding: '14px 18px',
    backgroundColor: 'rgba(255, 255, 255, 0.1)',
    borderRadius: '14px',
    border: '1px solid rgba(255, 255, 255, 0.15)',
    backdropFilter: 'blur(20px)',
    animation: 'floatCard 5s ease-in-out infinite 1s',
  },
  floatingCard3: {
    position: 'absolute',
    bottom: '20px',
    left: '-20px',
    display: 'flex',
    alignItems: 'center',
    gap: '12px',
    padding: '14px 18px',
    backgroundColor: 'rgba(255, 255, 255, 0.1)',
    borderRadius: '14px',
    border: '1px solid rgba(255, 255, 255, 0.15)',
    backdropFilter: 'blur(20px)',
    animation: 'floatCard 7s ease-in-out infinite 0.5s',
  },
  floatingCardIcon: {
    fontSize: '24px',
  },
  floatingCardValue: {
    fontSize: '18px',
    fontWeight: '700',
    color: '#ffffff',
  },
  floatingCardLabel: {
    fontSize: '12px',
    color: 'rgba(255, 255, 255, 0.6)',
  },
  features: {
    padding: '80px 48px',
    position: 'relative',
    zIndex: 10,
    maxWidth: '1200px',
    margin: '0 auto',
  },
  featuresHeader: {
    textAlign: 'center',
    marginBottom: '60px',
  },
  featuresTitle: {
    fontSize: '40px',
    fontWeight: '700',
    color: '#ffffff',
    margin: '0 0 16px 0',
    letterSpacing: '-0.02em',
  },
  featuresSubtitle: {
    fontSize: '18px',
    color: 'rgba(255, 255, 255, 0.6)',
    margin: 0,
  },
  featuresGrid: {
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))',
    gap: '24px',
  },
  featureCard: {
    padding: '32px',
    backgroundColor: 'rgba(255, 255, 255, 0.04)',
    borderRadius: '20px',
    border: '1px solid rgba(255, 255, 255, 0.08)',
    transition: 'all 0.3s ease',
  },
  featureIcon: {
    width: '56px',
    height: '56px',
    backgroundColor: 'rgba(255, 255, 255, 0.08)',
    borderRadius: '14px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    marginBottom: '20px',
  },
  featureTitle: {
    fontSize: '20px',
    fontWeight: '600',
    color: '#ffffff',
    margin: '0 0 12px 0',
  },
  featureDesc: {
    fontSize: '15px',
    lineHeight: '1.6',
    color: 'rgba(255, 255, 255, 0.6)',
    margin: 0,
  },
  cta: {
    padding: '60px 48px 80px',
    position: 'relative',
    zIndex: 10,
  },
  ctaCard: {
    maxWidth: '700px',
    margin: '0 auto',
    padding: '60px',
    backgroundColor: 'rgba(22, 163, 74, 0.1)',
    borderRadius: '28px',
    border: '1px solid rgba(22, 163, 74, 0.25)',
    textAlign: 'center',
    backdropFilter: 'blur(20px)',
  },
  ctaTitle: {
    fontSize: '36px',
    fontWeight: '700',
    color: '#ffffff',
    margin: '0 0 16px 0',
    letterSpacing: '-0.02em',
  },
  ctaText: {
    fontSize: '18px',
    color: 'rgba(255, 255, 255, 0.7)',
    margin: '0 0 32px 0',
  },
  ctaButton: {
    display: 'inline-flex',
    alignItems: 'center',
    padding: '18px 40px',
    fontSize: '17px',
    fontWeight: '600',
    background: 'linear-gradient(135deg, #22c55e 0%, #16a34a 100%)',
    border: 'none',
    borderRadius: '14px',
    color: '#ffffff',
    cursor: 'pointer',
    transition: 'all 0.2s ease',
    fontFamily: "'DM Sans', sans-serif",
    boxShadow: '0 12px 40px rgba(22, 163, 74, 0.35)',
  },
  footer: {
    padding: '32px 48px',
    borderTop: '1px solid rgba(255, 255, 255, 0.1)',
    position: 'relative',
    zIndex: 10,
  },
  footerContent: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    maxWidth: '1200px',
    margin: '0 auto',
  },
  footerLogo: {
    display: 'flex',
    alignItems: 'center',
    gap: '10px',
  },
  footerLogoText: {
    fontSize: '18px',
    fontWeight: '600',
    color: '#ffffff',
  },
  footerCopyright: {
    fontSize: '14px',
    color: 'rgba(255, 255, 255, 0.5)',
    margin: 0,
  },
};

export default App;

