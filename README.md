# MoleculeAI - Drug Discovery Platform

A modern React frontend for drug discovery analysis using SMILES molecular notation.

## Features

- **SMILES Input**: Enter molecular structures using SMILES notation
- **Analysis Prompts**: Describe your analysis goals (e.g., compare drug properties)
- **Request Queue**: Bundle multiple molecules with prompts for batch analysis
- **Results View**: View analysis results with metrics and AI insights
- **Modern UI**: Clean, sleek design with green and white color scheme

## Getting Started

### Prerequisites

- Node.js (v16 or higher)
- npm or yarn

### Installation

1. Extract the zip file
2. Navigate to the project directory:
   ```bash
   cd molecule-ai-app
   ```
3. Install dependencies:
   ```bash
   npm install
   ```
4. Start the development server:
   ```bash
   npm start
   ```
5. Open [http://localhost:3000](http://localhost:3000) in your browser

## Usage

1. **Add Molecules**: Enter SMILES codes in the input field and press Enter, or use the quick-add buttons
2. **Write Prompt**: Describe what you want to analyze (e.g., "Compare the drug properties of these molecules")
3. **Submit Request**: Click "Add to Analysis Queue" to bundle molecules + prompt together
4. **Run Analysis**: Click "Run All Analyses" or run individual requests
5. **View Results**: Navigate to the Results tab to see analysis output

## Project Structure

```
molecule-ai-app/
├── public/
│   └── index.html
├── src/
│   ├── App.jsx        # Main application component
│   ├── index.css      # Global styles
│   └── index.js       # Entry point
├── package.json
└── README.md
```

## Backend Integration

This is a frontend-only application. To connect to a backend:

1. Update the `handleAnalyzeRequest` and `handleAnalyzeAll` functions in `App.jsx`
2. Replace the setTimeout mock with actual API calls
3. The request payload will include:
   - `molecules`: Array of {id, name, smiles}
   - `prompt`: Analysis prompt string

## Tech Stack

- React 18
- Create React App
- Inline styles (no external CSS framework)

## License

MIT
