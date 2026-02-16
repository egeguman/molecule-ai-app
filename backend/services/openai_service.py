import json
import logging

from config.settings import OPENAI_API_KEY, OPENAI_MODEL

logger = logging.getLogger(__name__)


class OpenAIService:
    def __init__(self):
        self.client = None
        if OPENAI_API_KEY:
            try:
                from openai import OpenAI
                self.client = OpenAI(api_key=OPENAI_API_KEY)
                logger.info("OpenAI client initialized (model: %s)", OPENAI_MODEL)
            except ImportError:
                logger.warning("openai package not installed — using mock interpretation")

    def interpret(self, predictions: dict, user_prompt: str) -> dict:
        """Get GPT-4o interpretation of ADMET predictions.

        Args:
            predictions: dict keyed by SMILES, each value is a dict of endpoint results
            user_prompt: the user's analysis request text
        """
        if not self.client:
            return self._mock_interpretation(predictions)

        try:
            return self._call_openai(predictions, user_prompt)
        except Exception as e:
            logger.error("OpenAI API call failed: %s", e)
            return self._mock_interpretation(predictions)

    def _call_openai(self, predictions: dict, user_prompt: str) -> dict:
        system_prompt = """You are an expert medicinal chemist and pharmacologist specializing in ADMET property analysis for drug discovery.

Given ADMET prediction data for one or more molecules, provide a structured drug discovery interpretation.

Respond ONLY with valid JSON matching this exact schema:
{
  "ai_overview": "A comprehensive paragraph (5-8 sentences) summarizing the full analysis report: cover overall drug-likeness assessment, key ADMET strengths and weaknesses across absorption/distribution/metabolism/excretion/toxicity, highlight the most critical findings, and end with the overall recommendation. This should read as a standalone executive briefing that captures the essence of the entire report.",
  "executive_summary": "2-3 sentence overview of the molecule(s) drug-likeness and key ADMET findings",
  "high_risk_flags": ["list of concerning findings that could limit clinical development"],
  "admet_interpretation": [
    {
      "endpoint": "endpoint display name",
      "value": "predicted value with unit",
      "assessment": "favorable|moderate|concerning",
      "explanation": "brief pharmacological explanation"
    }
  ],
  "recommended_next_steps": ["actionable recommendations for further evaluation"],
  "disclaimer": "Standard disclaimer about computational predictions"
}

Focus on clinically relevant insights. Flag any values outside normal drug-like ranges.
For classification endpoints, values >0.7 indicate positive (active/toxic), <0.3 indicate negative."""

        user_message = f"""ADMET Predictions:
{json.dumps(predictions, indent=2)}

User's analysis request: {user_prompt}

Provide your structured interpretation."""

        response = self.client.chat.completions.create(
            model=OPENAI_MODEL,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_message},
            ],
            response_format={"type": "json_object"},
            temperature=0.3,
            max_tokens=3000,
        )

        return json.loads(response.choices[0].message.content)

    def _mock_interpretation(self, predictions: dict) -> dict:
        """Fallback when no OpenAI API key is configured."""
        # Build basic interpretation from predictions
        interpretations = []
        high_risk = []

        for smiles, endpoints in predictions.items():
            for name, data in endpoints.items():
                assessment = "favorable"
                if data["type"] == "classification" and data["value"] > 0.7:
                    assessment = "concerning"
                    high_risk.append(f"{data['display_name']}: high probability ({data['value']:.2f}) for {smiles}")
                elif data["type"] == "classification" and data["value"] > 0.4:
                    assessment = "moderate"

                interpretations.append({
                    "endpoint": data["display_name"],
                    "value": f"{data['value']} {data['unit'] or ''}".strip(),
                    "assessment": assessment,
                    "explanation": f"Mock assessment for {smiles}",
                })

        molecule_count = len(predictions)
        smiles_str = ", ".join(predictions.keys())

        return {
            "ai_overview": (
                f"This analysis covers {molecule_count} molecule(s) ({smiles_str}) evaluated across 19 ADMET endpoints "
                "spanning absorption, distribution, metabolism, excretion, and toxicity properties. "
                "Predictions were generated in mock mode using deterministic fallback values rather than trained ML models. "
                "As such, these results should not be used for decision-making. "
                "To obtain real AI-powered interpretation with clinically relevant insights, "
                "configure the OPENAI_API_KEY environment variable and optionally train Chemprop models. "
                "Once configured, this overview will provide a comprehensive drug-likeness assessment with actionable recommendations."
            ),
            "executive_summary": (
                f"Analysis of {molecule_count} molecule(s) ({smiles_str}). "
                "ADMET predictions generated in mock mode. "
                "Configure OPENAI_API_KEY in .env for AI-powered interpretation."
            ),
            "high_risk_flags": high_risk if high_risk else ["No high-risk flags detected (mock mode)"],
            "admet_interpretation": interpretations[:10],
            "recommended_next_steps": [
                "Configure OpenAI API key for detailed AI interpretation",
                "Train Chemprop models for real ADMET predictions (POST /api/train)",
                "Validate predictions experimentally",
            ],
            "disclaimer": (
                "These are computational predictions generated in mock mode and do not "
                "represent real ADMET data. All predictions should be validated through "
                "experimental assays before making any drug development decisions."
            ),
        }
