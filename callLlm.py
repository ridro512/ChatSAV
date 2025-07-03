#!/usr/bin/env python3
"""
callLlm.py - File for calling LLM analysis of SpliceAI and GTEx results
Integrates with OPENAI to provide genetic analysis

Input: SpliceAI results, GTEx results, variant coordinates
Output: LLM interpretation of SpliceAI and GTEx results
"""

import json
import os
from typing import Dict, Any, Optional, List

# configure api keys from environmental variables
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("ERROR: OPENAI_API_KEY environment variable not set")
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
- Gene ID: {transcript.get('ensembl_id', 'Unknown')}
- Splice Disruption Scores:
  * Acceptor Gain (DS_AG): {transcript.get('DS_AG', 0)}
  * Acceptor Loss (DS_AL): {transcript.get('DS_AL', 0)}
  * Donor Gain (DS_DG): {transcript.get('DS_DG', 0)}
  * Donor Loss (DS_DL): {transcript.get('DS_DL', 0)}
- Position Changes:
  * Acceptor Gain (DP_AG): {transcript.get('DP_AG', 0)}
  * Acceptor Loss (DP_AL): {transcript.get('DP_AL', 0)}
  * Donor Gain (DP_DG): {transcript.get('DP_DG', 0)}
  * Donor Loss (DP_DL): {transcript.get('DP_DL', 0)}"""
        formatted_data.append(transcript_info)
    return "\n".join(formatted_data)

# construct a comprehensive prompt for the LLM
# the structure follows: provides context for the variant, includes apliceAI and gtex data in json
# gives LLM a clear analysis framework, request specific output
# returning formatted prompt string for LLM
def construct_prompt(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
                     chrom: str, pos: int, ref: str, alt: str, context: Optional[str] = None) -> str:
    
    transcript_data = format_transcript(spliceai_results.get('splice_scores', []))

    # include context if provided
    context_section = ""
    if context and context.strip():
        context_section = f"""
USER CONTEXT:
{context.strip()}

Please consider this context when providing your analysis.
"""
        
    prompt = f"""You are an expert genetic analyst specializing in splice-altering variants (SAVs). 
Please analyze the following genetic variant and provide a focused, concise evaluation for laboratory follow-up.

CONSTRAINTS:
- Provide ONLY the most relevant and accurate analysis
- Focus on the transcript with the highest splice disruption scores
- Give ONE clear recommendation per section, not multiple options
- Be specific and actionable in your recommendations
- Limit responses to essential information only
- For experimental recommendations, suggest only the SINGLE most important/relevant test

VARIANT INFORMATION:
- Chromosome: {chrom}
- Position: {pos}
- Reference allele: {ref}
- Alternative allele: {alt}
- Variant notation: {chrom}:{pos}:{ref}:{alt}

SPLICEAI PREDICTION DATA:
{transcript_data}

GTEX EXPRESSION DATA:
{json.dumps(gtex_results, indent=2)}
{context_section}

ANALYSIS FRAMEWORK:
Provide a structured analysis with the following sections. Each section should contain ONLY the most relevant point:

1. PRIMARY SPLICING IMPACT:
   - Identify the transcript with highest disruption score (>0.5 is significant)
   - State the most likely splicing mechanism affected
   - One sentence summary of predicted impact

2. CLINICAL SIGNIFICANCE:
   - State whether this variant is likely pathogenic based on the evidence
   - Provide confidence level (High/Medium/Low)
   - One key reason supporting your assessment

3. PATHOGENICITY ASSESSMENT:
   - Overall likelihood of being disease-causing
   - Confidence level in prediction
   - Key evidence supporting or refuting pathogenicity

4. EXPERIMENTAL PRIORITY:
   - Assign priority level: High/Medium/Low
   - Recommend the SINGLE most appropriate and important testing method
   - Specify optimal tissue/cell type for testing

5. SUMMARY RECOMMENDATION:
   - One clear action item for the curator
   - Key consideration for experimental validation

For experimental recommendations, choose only the most critical test - do not list multiple options.
Keep each section to 2-3 sentences maximum. Focus on actionable insights rather than general explanations."""

    return prompt

# main function that calls OPENAI
# construct a prompt with variant and prediction data
# extract information from the LLM response
# return results in standardised format
def call_llm(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
             chrom: str = None, pos: int = None, ref: str = None, alt: str = None,
             context: Optional[str] = None,
             model: str = "gpt-4.1-nano", temperature: float = 0.3) -> Dict[str, Any]:
    if not client:
        return {
            "llm_interpretation": "OpenAI client not initialized. Please check your API key.",
            "priority_level": "unknown",
            "pathogenicity_assessment": "unable_to_assess",
            "experimental_recommendations": "manual_review_required",
            "error": "OpenAI client not available"
        }
    try:
        prompt = construct_prompt(spliceai_results, gtex_results, chrom, pos, ref, alt, context)
        # call openai
        response = client.chat.completions.create(
            model=model,
            messages=[
                {
                    "role": "system",
                    "content": "You are an expert genetic analyst. Provide focused, concise analysis with clear recommendations. Avoid redundancy and focus on the most relevant findings only."
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ],
            max_tokens=1024,
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
    if any(phrase in text_lower for phrase in ["high priority", "priority: high", "priority level: high"]):
        return "high"
    elif any(phrase in text_lower for phrase in ["medium priority", "priority: medium", "priority level: medium"]):
        return "medium"
    elif any(phrase in text_lower for phrase in ["low priority", "priority: low", "priority level: low"]):
        return "low"
    
    if any(phrase in text_lower for phrase in ["significant disruption", "likely pathogenic", "high confidence"]):
        return "high"
    elif any(phrase in text_lower for phrase in ["uncertain", "moderate", "possible"]):
        return "medium"
    else:
        return "low"

# extract pathogenicity, looks for standard pathogenicity terms following ACMG guideline
def extract_pathogenicity_assessment(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    if "likely pathogenic" in text_lower:
        return "likely_pathogenic"
    elif "possibly pathogenic" in text_lower or "probably pathogenic" in text_lower:
        return "possibly_pathogenic"
    elif "pathogenic" in text_lower and "not" not in text_lower.split("pathogenic")[0][-20:]:
        return "likely_pathogenic"
    elif "likely benign" in text_lower:
        return "likely_benign"
    elif "benign" in text_lower:
        return "likely_benign"
    elif any(phrase in text_lower for phrase in ["uncertain significance", "uncertain", "unclear"]):
        return "uncertain_significance"
    else:
        return "uncertain_significance"

# identifies which in vitro testing methods openai recommends
def extract_experimental_recommendations(analysis_text: str) -> str:
    """Identify the most important experimental recommendation from the LLM response."""
    text_lower = analysis_text.lower()
    
    recommendation_priority = [
        ("minigene_assay", ["minigene"]),
        ("rt_pcr", ["rt-pcr", "rt pcr"]),
        ("rna_seq", ["rna-seq", "rna seq", "transcriptome"]),
        ("patient_sample_testing", ["patient sample", "patient tissue"]),
        ("functional_assay", ["functional assay", "functional testing"]),
        ("whole_blood_testing", ["whole blood"]),
    ]
    
    for rec_code, keywords in recommendation_priority:
        for keyword in keywords:
            if keyword in text_lower:
                return rec_code
    
    return "manual_review_required"