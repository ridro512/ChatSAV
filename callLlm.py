#!/usr/bin/env python3
"""
callLlm.py - File for calling LLM analysis of SpliceAI, Pangolin, AlphaGenome and GTEx results
Integrates with OPENAI to provide genetic analysis

Input: SpliceAI results, Pangolin results, AlphaGenome results, GTEx results, variant coordinates
Output: LLM interpretation of all prediction results
"""

import json
import os
from typing import Dict, Any, Optional, List


# configure api keys from environmental variables
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("ERROR: OPENAI_API_KEY environment variable not set")
    print("Please set your OpenAI API key: export OPENAI_API_KEY='your_api_key_here'")
    client = None
else:
    try:
        # import and initialise only after check for api key
        from openai import OpenAI
        client = OpenAI(api_key=api_key)
        print("OpenAI client initialized successfully")
    except ImportError:
        print("ERROR: openai package not installed. Run: pip install openai")
        client = None

def format_transcript(splice_scores: List[Dict[str, Any]]) -> str:
    if not splice_scores:
        return "No transcript data available"
    
    formatted_data = []
    for i, transcript in enumerate(splice_scores, 1):
        transcript_info = f"""
        Transcript {i}:
        - Transcript ID: {transcript.get('transcript_id', 'Unknown')}
        - Gene Name: {transcript.get('g_name', 'Unknown')}
        - Ensembl Gene ID: {transcript.get('ensembl_id', 'Unknown')}
        - Transcript Type: {transcript.get('t_type', 'Unknown')}
        - Transcript Priority: {transcript.get('t_priority', 'Unknown')}
        - Strand: {transcript.get('t_strand', 'Unknown')}
        - RefSeq IDs: {', '.join(transcript.get('t_refseq_ids', [])) if transcript.get('t_refseq_ids') else 'None'}

        Splice Disruption Scores (DS):
        * Acceptor Gain (DS_AG): {transcript.get('DS_AG', 0)}
            - Reference: {transcript.get('DS_AG_REF', 0)}
            - Alternate: {transcript.get('DS_AG_ALT', 0)}
        * Acceptor Loss (DS_AL): {transcript.get('DS_AL', 0)}
            - Reference: {transcript.get('DS_AL_REF', 0)}
            - Alternate: {transcript.get('DS_AL_ALT', 0)}
        * Donor Gain (DS_DG): {transcript.get('DS_DG', 0)}
            - Reference: {transcript.get('DS_DG_REF', 0)}
            - Alternate: {transcript.get('DS_DG_ALT', 0)}
        * Donor Loss (DS_DL): {transcript.get('DS_DL', 0)}
            - Reference: {transcript.get('DS_DL_REF', 0)}
            - Alternate: {transcript.get('DS_DL_ALT', 0)}

        Position Changes (DP in bp from variant site):
        * Acceptor Gain (DP_AG): {transcript.get('DP_AG', 0)}
        * Acceptor Loss (DP_AL): {transcript.get('DP_AL', 0)}
        * Donor Gain (DP_DG): {transcript.get('DP_DG', 0)}
        * Donor Loss (DP_DL): {transcript.get('DP_DL', 0)}
        """
        formatted_data.append(transcript_info.strip())
    return "\n\n".join(formatted_data)

def format_pangolin_data(pangolin_results: Optional[Dict[str, Any]]) -> str:
    """Format Pangolin results for LLM prompt."""
    if not pangolin_results or 'error' in pangolin_results:
        return "No Pangolin data available"
    
    formatted_data = []
    for i, transcript in enumerate(pangolin_results.get('pangolin_scores', []), 1):
        transcript_info = f"""
        Pangolin Transcript {i}:
        - Transcript ID: {transcript.get('transcript_id', 'Unknown')}
        - Gene Name: {transcript.get('g_name', 'Unknown')}
        - Ensembl Gene ID: {transcript.get('ensembl_id', 'Unknown')}
        - Transcript Type: {transcript.get('t_type', 'Unknown')}
        - Transcript Priority: {transcript.get('t_priority', 'Unknown')}
        - Strand: {transcript.get('t_strand', 'Unknown')}
        - RefSeq IDs: {', '.join(transcript.get('t_refseq_ids', [])) if transcript.get('t_refseq_ids') else 'None'}

        Splice Disruption Scores (DS):
        * Splice Gain (DS_SG): {transcript.get('DS_SG', 0)}
            - Reference: {transcript.get('SG_REF', 0)}
            - Alternate: {transcript.get('SG_ALT', 0)}
        * Splice Loss (DS_SL): {transcript.get('DS_SL', 0)}
            - Reference: {transcript.get('SL_REF', 0)}
            - Alternate: {transcript.get('SL_ALT', 0)}

        Position Changes (DP in bp from variant site):
        * Splice Gain (DP_SG): {transcript.get('DP_SG', 0)}
        * Splice Loss (DP_SL): {transcript.get('DP_SL', 0)}
        """
        formatted_data.append(transcript_info.strip())
    return "\n\n".join(formatted_data)

def format_alphagenome_data(alphagenome_results: Optional[Dict[str, Any]]) -> str:
    """Format AlphaGenome results for LLM prompt."""
    if not alphagenome_results or 'error' in alphagenome_results:
        return "No AlphaGenome data available"
    
    scores = alphagenome_results.get('alphagenome_scores', {})
    formatted_parts = []
    
    # summary statistics
    summary = scores.get('summary_stats', {})
    if summary and 'error' not in summary:
        formatted_parts.append(f"""
AlphaGenome Summary:
- Total predictions: {summary.get('total_predictions', 0)}
- Significant effects (|score| > 0.5): {summary.get('significant_effects', 0)}
- Mean quantile score: {summary.get('mean_quantile_score', 0):.3f}
- Max absolute effect: {summary.get('max_absolute_effect', 0):.3f}
- Unique tissues analyzed: {summary.get('unique_tissues', 0)}""")
    
    # splice predictions
    splice_preds = scores.get('splice_predictions', {})
    if splice_preds and 'error' not in splice_preds:
        formatted_parts.append("\nAlphaGenome Splice Predictions:")
        for splice_type, data in splice_preds.items():
            if isinstance(data, dict) and 'total_predictions' in data:
                formatted_parts.append(f"""
- {splice_type.replace('_', ' ').title()}:
  * Total predictions: {data.get('total_predictions', 0)}
  * Significant predictions: {data.get('significant_predictions', 0)}
  * Mean quantile score: {data.get('mean_quantile', 0):.3f}""")
    
    # expression predictions
    expr_preds = scores.get('expression_predictions', {})
    if expr_preds and 'error' not in expr_preds:
        formatted_parts.append("\nAlphaGenome Expression Predictions:")
        for expr_type, data in expr_preds.items():
            if isinstance(data, dict) and 'total_predictions' in data:
                formatted_parts.append(f"""
- {expr_type.upper()}:
  * Total predictions: {data.get('total_predictions', 0)}
  * Upregulated: {data.get('upregulated_count', 0)}
  * Downregulated: {data.get('downregulated_count', 0)}
  * Mean effect: {data.get('mean_effect', 0):.3f}""")
    
    # chromatin predictions
    chromatin_preds = scores.get('chromatin_predictions', {})
    if chromatin_preds and 'error' not in chromatin_preds:
        formatted_parts.append("\nAlphaGenome Chromatin Accessibility:")
        for chrom_type, data in chromatin_preds.items():
            if isinstance(data, dict) and 'total_predictions' in data:
                formatted_parts.append(f"""
- {chrom_type.upper()}:
  * Total predictions: {data.get('total_predictions', 0)}
  * Accessibility increase: {data.get('accessibility_increase', 0)}
  * Accessibility decrease: {data.get('accessibility_decrease', 0)}""")
    
    # top predictions
    top_preds = alphagenome_results.get('top_predictions', [])
    if top_preds:
        formatted_parts.append("\nTop AlphaGenome Predictions:")
        for i, pred in enumerate(top_preds[:5], 1):  # Top 5
            formatted_parts.append(f"""
{i}. {pred.get('output_type', 'Unknown')} in {pred.get('biosample_name', 'Unknown')}:
   - Quantile score: {pred.get('quantile_score', 0):.3f}
   - Effect: {pred.get('effect_direction', 'unknown')}
   - Significance: {pred.get('significance', 'unknown')}""")
    
    return "\n".join(formatted_parts) if formatted_parts else "AlphaGenome analysis completed but no interpretable results available"

# construct a comprehensive prompt for the LLM
# the structure follows: provides context for the variant, includes spliceAI, pangolin, alphagenome and gtex data
# gives LLM a clear analysis framework, request specific output
# returning formatted prompt string for LLM
def construct_prompt(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
                     pangolin_results: Optional[Dict[Any, Any]] = None,
                     alphagenome_results: Optional[Dict[Any, Any]] = None,
                     chrom: str = None, pos: int = None, ref: str = None, alt: str = None, 
                     context: Optional[str] = None) -> str:
    
    transcript_data = format_transcript(spliceai_results.get('splice_scores', []))
    pangolin_data = format_pangolin_data(pangolin_results)
    alphagenome_data = format_alphagenome_data(alphagenome_results)

    # include context if provided
    context_section = ""
    if context and context.strip():
        context_section = f"""
USER CONTEXT:
{context.strip()}

Please consider this context when providing your analysis.
"""
        
    prompt = f"""You are an expert genetic analyst specializing in splice-altering variants (SAVs). 
Please analyze the following genetic variant using multiple prediction tools and provide a focused, concise evaluation for laboratory follow-up.

CONSTRAINTS:
- Provide ONLY the most relevant and accurate analysis
- Focus on the most significant findings across all prediction tools
- Give ONE clear recommendation per section, not multiple options
- Be specific and actionable in your recommendations
- Limit responses to essential information only
- For experimental recommendations, suggest only the SINGLE most important/relevant test
- When tools disagree, explain which prediction is more reliable and why

VARIANT INFORMATION:
- Chromosome: {chrom}
- Position: {pos}
- Reference allele: {ref}
- Alternative allele: {alt}
- Variant notation: {chrom}:{pos}:{ref}:{alt}

SPLICEAI PREDICTION DATA:
{transcript_data}

PANGOLIN PREDICTION DATA:
{pangolin_data}

ALPHAGENOME PREDICTION DATA:
{alphagenome_data}

GTEX EXPRESSION DATA:
{json.dumps(gtex_results, indent=2)}
{context_section}

ANALYSIS FRAMEWORK:
Provide a structured analysis with the following sections. Each section should contain ONLY the most relevant point:

1. COMPARATIVE SPLICE IMPACT ANALYSIS:
   - Compare SpliceAI, Pangolin, and AlphaGenome predictions
   - Identify the most reliable prediction and explain why
   - State the consensus prediction or highlight key disagreements
   - Focus on the transcript/gene with strongest evidence

2. MULTI-MODAL EVIDENCE INTEGRATION:
   - Integrate splice predictions with AlphaGenome expression and chromatin data
   - Assess whether expression changes support splice predictions
   - Consider tissue-specific effects from AlphaGenome analysis
   - One sentence summary of integrated evidence

3. CLINICAL SIGNIFICANCE:
   - State whether this variant is likely pathogenic based on ALL evidence
   - Provide confidence level (High/Medium/Low) considering tool agreement
   - Key reason supporting assessment using strongest evidence

4. PATHOGENICITY ASSESSMENT:
   - Overall likelihood of being disease-causing
   - Confidence level based on tool consensus
   - Most compelling evidence for or against pathogenicity

5. EXPERIMENTAL PRIORITY:
   - Assign priority level: High/Medium/Low
   - Recommend the SINGLE most appropriate testing method
   - Specify optimal tissue/cell type considering AlphaGenome tissue predictions
   - Consider which experimental approach would best validate the predictions

6. SUMMARY RECOMMENDATION:
   - One clear action item for the curator
   - Key consideration for experimental validation
   - Note any important limitations or uncertainties

SCORING INTERPRETATION:
- SpliceAI/Pangolin scores >0.5 are significant, >0.8 are high confidence
- AlphaGenome quantile scores >0.5 indicate significant tissue effects
- When tools disagree, prioritize: (1) Tool agreement, (2) AlphaGenome multi-modal evidence, (3) Higher confidence scores

For experimental recommendations, choose only the most critical test based on the strongest prediction evidence.
Keep each section to 2-3 sentences maximum. Focus on actionable insights rather than general explanations."""

    return prompt

# main function that calls OPENAI
# construct a prompt with variant and prediction data from all tools
# extract information from the LLM response
# return results in standardised format
def call_llm(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
             pangolin_results: Optional[Dict[Any, Any]] = None,
             alphagenome_results: Optional[Dict[Any, Any]] = None,
             chrom: str = None, pos: int = None, ref: str = None, alt: str = None,
             context: Optional[str] = None,
             model: str = "gpt-4o", temperature: float = 0.3) -> Dict[str, Any]:
    if not client:
        return {
            "llm_interpretation": "OpenAI client not initialized. Please check your API key.",
            "priority_level": "unknown",
            "pathogenicity_assessment": "unable_to_assess",
            "experimental_recommendations": "manual_review_required",
            "error": "OpenAI client not available"
        }
    try:
        prompt = construct_prompt(spliceai_results, gtex_results, pangolin_results, 
                                alphagenome_results, chrom, pos, ref, alt, context)
        
        # determine system message based on available tools
        available_tools = []
        if spliceai_results and 'error' not in spliceai_results:
            available_tools.append("SpliceAI")
        if pangolin_results and 'error' not in pangolin_results:
            available_tools.append("Pangolin")
        if alphagenome_results and 'error' not in alphagenome_results:
            available_tools.append("AlphaGenome")
        if gtex_results and 'error' not in gtex_results:
            available_tools.append("GTEx")
        
        system_content = f"""You are an expert genetic analyst with access to multiple prediction tools: {', '.join(available_tools)}. 

Provide focused, concise analysis with clear recommendations. When multiple tools are available:
1. Compare and contrast their predictions
2. Identify consensus or explain disagreements
3. Weight evidence based on tool reliability and confidence scores
4. Integrate multi-modal evidence (splice + expression + chromatin when available)

Avoid redundancy and focus on the most relevant findings that inform clinical decision-making."""
        
        # call openai
        response = client.chat.completions.create(
            model=model,
            messages=[
                {
                    "role": "system",
                    "content": system_content
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ],
            max_tokens=1200,
            temperature=temperature
        )
        # extract text response
        analysis_text = response.choices[0].message.content

        # parse the response, helper functions scan the text for specific keywords and patterns
        priority = extract_priority_level(analysis_text)
        pathogenicity = extract_pathogenicity_assessment(analysis_text)
        
        return {
            "llm_interpretation": analysis_text,
            "priority_level": priority,
            "pathogenicity_assessment": pathogenicity,
            "experimental_recommendations": extract_experimental_recommendations(analysis_text),
            "model_used": model,
            "tools_analyzed": available_tools,
            "prompt_used": prompt,
            "context_provided": context if context else "None"
        }        
    except Exception as e:
        return {
            "llm_interpretation": f"Error calling OpenAI API: {str(e)}",
            "priority_level": "unknown",
            "pathogenicity_assessment": "unable_to_assess",
            "experimental_recommendations": "manual_review_required",
            "error": str(e)
        }


# these helpers extract information from openai's natural language response

# scan LLM response for keywords indicating recommend priority for lab follow-up
# medium by default
def extract_priority_level(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    
    # check for explicit priority statements
    if any(phrase in text_lower for phrase in ["high priority", "priority: high", "priority level: high"]):
        return "high"
    elif any(phrase in text_lower for phrase in ["medium priority", "priority: medium", "priority level: medium"]):
        return "medium"
    elif any(phrase in text_lower for phrase in ["low priority", "priority: low", "priority level: low"]):
        return "low"
    
    # check for implicit high priority indicators
    high_priority_indicators = [
        "significant disruption", "likely pathogenic", "high confidence", "strong evidence",
        "consensus prediction", "multiple tools agree", "tissue-specific", "functional impact"
    ]
    if any(phrase in text_lower for phrase in high_priority_indicators):
        return "high"
    
    # check for implicit low priority indicators
    low_priority_indicators = [
        "conflicting predictions", "low confidence", "minimal effect", "unlikely pathogenic",
        "benign", "no significant"
    ]
    if any(phrase in text_lower for phrase in low_priority_indicators):
        return "low"
    
    # check for uncertainty indicators
    if any(phrase in text_lower for phrase in ["uncertain", "moderate", "possible", "disagreement"]):
        return "medium"
    
    return "medium"

# extract pathogenicity, looks for standard pathogenicity terms following ACMG guideline
def extract_pathogenicity_assessment(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    
    # check for explicit pathogenicity statements
    if "likely pathogenic" in text_lower:
        return "likely_pathogenic"
    elif "possibly pathogenic" in text_lower or "probably pathogenic" in text_lower:
        return "possibly_pathogenic"
    elif "pathogenic" in text_lower and "not" not in text_lower.split("pathogenic")[0][-20:]:
        return "likely_pathogenic"
    elif "likely benign" in text_lower:
        return "likely_benign"
    elif "benign" in text_lower and "not" not in text_lower.split("benign")[0][-20:]:
        return "likely_benign"
    
    # check for uncertainty indicators
    uncertainty_indicators = [
        "uncertain significance", "uncertain", "unclear", "conflicting evidence",
        "tool disagreement", "moderate confidence"
    ]
    if any(phrase in text_lower for phrase in uncertainty_indicators):
        return "uncertain_significance"
    
    # check for pathogenic indicators
    pathogenic_indicators = [
        "disease-causing", "functional impact", "splice disruption", "expression change",
        "significant effect", "consensus pathogenic"
    ]
    if any(phrase in text_lower for phrase in pathogenic_indicators):
        return "likely_pathogenic"
    
    # check for benign indicators
    benign_indicators = [
        "no functional impact", "minimal effect", "unlikely disease", "normal expression"
    ]
    if any(phrase in text_lower for phrase in benign_indicators):
        return "likely_benign"
    
    return "uncertain_significance"

# identifies which in vitro testing methods openai recommends
def extract_experimental_recommendations(analysis_text: str) -> str:
    """Identify the most important experimental recommendation from the LLM response."""
    text_lower = analysis_text.lower()
    
    # priority order for experimental recommendations
    recommendation_priority = [
        ("minigene_assay", ["minigene", "splicing reporter"]),
        ("rt_pcr", ["rt-pcr", "rt pcr", "reverse transcription"]),
        ("rna_seq", ["rna-seq", "rna seq", "transcriptome", "rna sequencing"]),
        ("patient_sample_testing", ["patient sample", "patient tissue", "clinical sample"]),
        ("functional_assay", ["functional assay", "functional testing", "in vitro assay"]),
        ("cell_culture_testing", ["cell culture", "tissue culture", "cultured cells"]),
        ("whole_blood_testing", ["whole blood", "blood sample"]),
        ("tissue_specific_testing", ["tissue-specific", "tissue specific"]),
        ("alphagenome_validation", ["alphagenome", "multi-modal", "chromatin"]),
    ]
    
    for rec_code, keywords in recommendation_priority:
        for keyword in keywords:
            if keyword in text_lower:
                return rec_code
    
    # check for general experimental terms
    if any(term in text_lower for term in ["experiment", "testing", "validation", "assay"]):
        return "functional_assay"
    
    return "manual_review_required"