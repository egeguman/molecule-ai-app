import React, { useState } from 'react';
import './index.css';

const App = () => {
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

const styles = {
  container: {
    minHeight: '100vh',
    backgroundColor: '#f8fafc',
    fontFamily: "'DM Sans', sans-serif",
    color: '#1f2937',
    display: 'flex',
    flexDirection: 'column',
  },
  header: {
    backgroundColor: '#ffffff',
    borderBottom: '1px solid #e5e7eb',
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
    backgroundColor: '#f0fdf4',
    borderRadius: '12px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
  },
  logoTitle: {
    fontSize: '20px',
    fontWeight: '700',
    color: '#111827',
    margin: 0,
    letterSpacing: '-0.02em',
  },
  logoSubtitle: {
    fontSize: '13px',
    color: '#6b7280',
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
    border: '1px solid #e5e7eb',
    borderRadius: '8px',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    color: '#4b5563',
    fontFamily: "'DM Sans', sans-serif",
  },
  navButtonActive: {
    backgroundColor: '#16a34a',
    borderColor: '#16a34a',
    color: '#ffffff',
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
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
    borderRadius: '16px',
    padding: '28px',
  },
  card: {
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
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
    color: '#111827',
    margin: 0,
  },
  badge: {
    fontSize: '12px',
    fontWeight: '500',
    color: '#16a34a',
    backgroundColor: '#f0fdf4',
    padding: '4px 10px',
    borderRadius: '20px',
  },
  countBadge: {
    fontSize: '13px',
    fontWeight: '500',
    color: '#6b7280',
    backgroundColor: '#f3f4f6',
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
    color: '#374151',
  },
  sectionCount: {
    fontSize: '12px',
    fontWeight: '500',
    color: '#16a34a',
    backgroundColor: '#f0fdf4',
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
    color: '#374151',
  },
  dropdown: {
    padding: '8px 12px',
    fontSize: '14px',
    fontFamily: "'DM Sans', sans-serif",
    border: '1px solid #e5e7eb',
    borderRadius: '8px',
    backgroundColor: '#ffffff',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    color: '#111827',
    fontWeight: '500',
  },
  smilesInputFieldFull: {
    width: '100%',
    padding: '12px 16px',
    fontSize: '14px',
    fontFamily: "'JetBrains Mono', monospace",
    border: '1px solid #e5e7eb',
    borderRadius: '8px',
    backgroundColor: '#fafafa',
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
    color: '#6b7280',
  },
  chip: {
    padding: '5px 12px',
    fontSize: '12px',
    fontWeight: '500',
    backgroundColor: '#f0fdf4',
    border: '1px solid #bbf7d0',
    borderRadius: '16px',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    color: '#16a34a',
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
    backgroundColor: '#fafafa',
    borderRadius: '8px',
    border: '1px solid #f3f4f6',
  },
  moleculeTagIndex: {
    width: '24px',
    height: '24px',
    backgroundColor: '#16a34a',
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
    color: '#111827',
    display: 'block',
    marginBottom: '2px',
  },
  moleculeTagSmiles: {
    fontSize: '11px',
    fontFamily: "'JetBrains Mono', monospace",
    color: '#6b7280',
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
    color: '#16a34a',
    textTransform: 'uppercase',
    letterSpacing: '0.5px',
  },
  moleculeTagRemove: {
    width: '24px',
    height: '24px',
    backgroundColor: 'transparent',
    border: '1px solid #e5e7eb',
    borderRadius: '6px',
    cursor: 'pointer',
    fontSize: '16px',
    color: '#9ca3af',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    transition: 'all 0.15s ease',
    flexShrink: 0,
    fontFamily: "'DM Sans', sans-serif",
  },
  divider: {
    height: '1px',
    backgroundColor: '#e5e7eb',
    margin: '24px 0',
  },
  promptSection: {
    marginBottom: '24px',
  },
  promptReady: {
    fontSize: '11px',
    fontWeight: '600',
    color: '#ffffff',
    backgroundColor: '#16a34a',
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
    border: '1px solid #e5e7eb',
    borderRadius: '8px',
    backgroundColor: '#fafafa',
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
    backgroundColor: '#f3f4f6',
    border: '1px solid #e5e7eb',
    borderRadius: '16px',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    color: '#4b5563',
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
    backgroundColor: '#16a34a',
    color: '#ffffff',
    border: 'none',
    borderRadius: '10px',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    fontFamily: "'DM Sans', sans-serif",
  },
  buttonDisabled: {
    backgroundColor: '#d1d5db',
    cursor: 'not-allowed',
  },
  submitHint: {
    fontSize: '12px',
    color: '#9ca3af',
    textAlign: 'center',
    marginTop: '10px',
    margin: '10px 0 0 0',
  },
  queueSection: {
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
    borderRadius: '16px',
    padding: '24px',
  },
  queueSectionHeader: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '20px',
    paddingBottom: '16px',
    borderBottom: '1px solid #f3f4f6',
  },
  queueTitleRow: {
    display: 'flex',
    alignItems: 'center',
    gap: '10px',
  },
  queueSectionTitle: {
    fontSize: '16px',
    fontWeight: '600',
    color: '#111827',
    margin: 0,
  },
  queueCount: {
    fontSize: '12px',
    fontWeight: '600',
    color: '#16a34a',
    backgroundColor: '#f0fdf4',
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
    backgroundColor: '#16a34a',
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
    backgroundColor: '#ffffff',
    color: '#6b7280',
    border: '1px solid #e5e7eb',
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
    backgroundColor: '#fafafa',
    borderRadius: '10px',
    border: '1px dashed #e5e7eb',
  },
  queueEmptyTextInline: {
    fontSize: '14px',
    fontWeight: '500',
    color: '#6b7280',
    margin: '0 0 2px 0',
  },
  queueEmptyHintInline: {
    fontSize: '12px',
    color: '#9ca3af',
    margin: 0,
  },
  queueGrid: {
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))',
    gap: '16px',
  },
  queueCard: {
    backgroundColor: '#fafafa',
    border: '1px solid #f3f4f6',
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
    color: '#111827',
  },
  queueCardMolCount: {
    fontSize: '11px',
    fontWeight: '500',
    color: '#16a34a',
    backgroundColor: '#f0fdf4',
    padding: '2px 6px',
    borderRadius: '8px',
  },
  queueCardRemove: {
    width: '22px',
    height: '22px',
    backgroundColor: 'transparent',
    border: '1px solid #e5e7eb',
    borderRadius: '4px',
    cursor: 'pointer',
    fontSize: '14px',
    color: '#9ca3af',
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
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
    borderRadius: '10px',
    color: '#374151',
  },
  queueCardPrompt: {
    fontSize: '12px',
    lineHeight: '1.4',
    color: '#6b7280',
    margin: '0 0 12px 0',
  },
  queueCardRunButton: {
    width: '100%',
    padding: '8px 12px',
    fontSize: '12px',
    fontWeight: '600',
    backgroundColor: '#16a34a',
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
    color: '#6b7280',
    margin: '0 0 16px 0',
  },
  primaryButton: {
    padding: '12px 24px',
    fontSize: '14px',
    fontWeight: '600',
    backgroundColor: '#16a34a',
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
    backgroundColor: '#ffffff',
    color: '#4b5563',
    border: '1px solid #e5e7eb',
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
    backgroundColor: '#fafafa',
    border: '1px solid #f3f4f6',
    borderRadius: '12px',
    padding: '20px',
  },
  requestHeader: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '16px',
    paddingBottom: '12px',
    borderBottom: '1px solid #e5e7eb',
  },
  requestTitleRow: {
    display: 'flex',
    alignItems: 'center',
    gap: '12px',
  },
  requestNumber: {
    fontSize: '14px',
    fontWeight: '700',
    color: '#111827',
  },
  requestTime: {
    fontSize: '12px',
    color: '#9ca3af',
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
    backgroundColor: '#f0fdf4',
    color: '#16a34a',
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
    color: '#6b7280',
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
    backgroundColor: '#ffffff',
    borderRadius: '6px',
    border: '1px solid #e5e7eb',
  },
  requestMolName: {
    fontSize: '13px',
    fontWeight: '600',
    color: '#111827',
    display: 'block',
    marginBottom: '2px',
  },
  requestMolSmiles: {
    fontSize: '10px',
    fontFamily: "'JetBrains Mono', monospace",
    color: '#6b7280',
    display: 'block',
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap',
  },
  requestPrompt: {},
  requestPromptText: {
    fontSize: '14px',
    lineHeight: '1.6',
    color: '#374151',
    margin: 0,
    padding: '10px 12px',
    backgroundColor: '#ffffff',
    borderRadius: '6px',
    border: '1px solid #e5e7eb',
  },
  requestActions: {
    display: 'flex',
    gap: '8px',
    paddingTop: '12px',
    borderTop: '1px solid #e5e7eb',
  },
  actionButton: {
    padding: '8px 16px',
    fontSize: '13px',
    fontWeight: '500',
    backgroundColor: '#16a34a',
    color: '#ffffff',
    border: 'none',
    borderRadius: '6px',
    cursor: 'pointer',
    transition: 'all 0.15s ease',
    fontFamily: "'DM Sans', sans-serif",
  },
  dangerButton: {
    backgroundColor: '#ffffff',
    color: '#dc2626',
    border: '1px solid #fecaca',
  },
  bulkActions: {
    display: 'flex',
    gap: '12px',
    marginTop: '24px',
    paddingTop: '20px',
    borderTop: '1px solid #e5e7eb',
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
    color: '#111827',
    margin: '0 0 8px 0',
  },
  resultsText: {
    fontSize: '15px',
    color: '#6b7280',
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
    backgroundColor: '#fafafa',
    border: '1px solid #f3f4f6',
    borderRadius: '12px',
    padding: '24px',
  },
  resultHeader: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '20px',
    paddingBottom: '12px',
    borderBottom: '1px solid #e5e7eb',
  },
  resultTitle: {
    fontSize: '16px',
    fontWeight: '600',
    color: '#111827',
  },
  statusCompletedBadge: {
    fontSize: '11px',
    fontWeight: '600',
    padding: '4px 10px',
    borderRadius: '12px',
    textTransform: 'uppercase',
    letterSpacing: '0.5px',
    backgroundColor: '#f0fdf4',
    color: '#16a34a',
  },
  requestSummary: {
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
    borderRadius: '10px',
    marginBottom: '20px',
    overflow: 'hidden',
  },
  summaryHeader: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    padding: '12px 16px',
    backgroundColor: '#f9fafb',
    borderBottom: '1px solid #e5e7eb',
  },
  summaryTitle: {
    fontSize: '13px',
    fontWeight: '600',
    color: '#374151',
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
    color: '#6b7280',
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
    backgroundColor: '#f0fdf4',
    border: '1px solid #bbf7d0',
    borderRadius: '14px',
    color: '#16a34a',
  },
  summaryPromptSection: {},
  summaryPromptText: {
    fontSize: '13px',
    lineHeight: '1.5',
    color: '#374151',
    margin: 0,
    padding: '10px 12px',
    backgroundColor: '#f9fafb',
    borderRadius: '6px',
    border: '1px solid #f3f4f6',
  },
  resultMetrics: {
    display: 'grid',
    gridTemplateColumns: 'repeat(3, 1fr)',
    gap: '12px',
    marginBottom: '20px',
  },
  metricCard: {
    backgroundColor: '#ffffff',
    border: '1px solid #e5e7eb',
    borderRadius: '8px',
    padding: '14px',
  },
  metricLabel: {
    fontSize: '12px',
    fontWeight: '500',
    color: '#6b7280',
    display: 'block',
    marginBottom: '6px',
  },
  metricValue: {
    fontSize: '18px',
    fontWeight: '700',
    color: '#111827',
    display: 'block',
    marginBottom: '10px',
  },
  progressBar: {
    height: '6px',
    backgroundColor: '#e5e7eb',
    borderRadius: '3px',
    overflow: 'hidden',
  },
  progressFill: {
    height: '100%',
    backgroundColor: '#16a34a',
    borderRadius: '3px',
    transition: 'width 0.3s ease',
  },
  resultInsight: {
    paddingTop: '16px',
    borderTop: '1px solid #e5e7eb',
  },
  insightLabel: {
    fontSize: '11px',
    fontWeight: '600',
    color: '#6b7280',
    textTransform: 'uppercase',
    letterSpacing: '0.5px',
    display: 'block',
    marginBottom: '8px',
  },
  insightText: {
    fontSize: '14px',
    lineHeight: '1.6',
    color: '#374151',
    margin: 0,
    fontStyle: 'italic',
  },
  footer: {
    padding: '20px 32px',
    borderTop: '1px solid #e5e7eb',
    textAlign: 'center',
    backgroundColor: '#ffffff',
  },
  footerText: {
    fontSize: '13px',
    color: '#9ca3af',
    margin: 0,
  },
};

export default App;
