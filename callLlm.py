#!/usr/bin/env python3
"""
callLlm.py - File for calling LLM analysis of SpliceAI and GTEx results
Integrates with Claude AI to provide genetic analysis

Input: SpliceAI results, GTEx results, variant coordinates
Output: LLM interpretation of SpliceAI and GTEx results
"""

import json
import os
from typing import Dict, Any, Optional

# configure api keys from environmental variables
api_key = os.getenv("ANTHROPIC_API_KEY")
if not api_key:
    print("ERROR: ANTHROPIC_API_KEY environment variable not set")
    client = None
else:
    try:
        # import and initialise only after check for api key
        import anthropic
        client = anthropic.Anthropic(api_key=api_key)
        print("Claude client initialized successfully")
    except ImportError:
        print("ERROR: anthropic package not installed. Run: pip install anthropic")
        client = None

# construct a comprehensive prompt for the LLM
# the structure follows: provides context for the variant, includes apliceAI and gtex data in json
# gives LLM a clear analysis framework, request specific output
# returning formatted prompt string for LLM
def construct_prompt(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
                     chrom: str, pos: int, ref: str, alt: str) -> str:
    
    prompt = f"""You are an expert genetic analyst specializing in splice-altering variants (SAVs). 
Please analyze the following genetic variant and provide a comprehensive evaluation for laboratory follow-up.

VARIANT INFORMATION:
- Chromosome: {chrom}
- Position: {pos}
- Reference allele: {ref}
- Alternative allele: {alt}
- Variant notation: {chrom}:{pos}:{ref}:{alt}

SPLICEAI PREDICTION DATA:
{json.dumps(spliceai_results, indent=2)}

GTEX EXPRESSION DATA:
{json.dumps(gtex_results, indent=2)}

ANALYSIS FRAMEWORK:
Please evaluate this variant based on the following criteria:

1. SPLICING DISRUPTION LIKELIHOOD:
   - Interpret SpliceAI prediction scores (note: scores >0.5 are considered significant)
   - Assess donor/acceptor site gain/loss probabilities
   - Determine likely splicing mechanism affected

2. CLINICAL RELEVANCE:
   - Evaluate tissue-specific expression relevance
   - Consider whether expression pattern matches potential disease phenotypes
   - Assess if whole blood testing would be informative

3. PATHOGENICITY ASSESSMENT:
   - Overall likelihood of being disease-causing
   - Confidence level in prediction
   - Key evidence supporting or refuting pathogenicity

4. EXPERIMENTAL RECOMMENDATIONS:
   - Recommend appropriate in vitro testing methods (minigene assay, RT-PCR, etc.)
   - Suggest optimal tissue/cell types for testing
   - Priority level for laboratory follow-up (High/Medium/Low)

5. SUMMARY:
   - Concise recommendation for curator action
   - Key points to consider for experimental validation

Please provide a structured analysis addressing each of these points. Be specific about the biological implications and practical next steps."""

    return prompt

# main function that calls Claude AI
# construct a prompt with variant and prediction data
# extract information from the LLM response
# return results in standardised format
def call_llm(spliceai_results: Dict[Any, Any], gtex_results: Dict[Any, Any],
             chrom: str = None, pos: int = None, ref: str = None, alt: str = None,
             model: str = "claude-3-5-sonnet-20241022", temperature: float = 0.3) -> Dict[str, Any]:
    try:
        prompt = construct_prompt(spliceai_results, gtex_results, chrom, pos, ref, alt)
        # call claude
        response = client.messages.create(
            model=model,
            max_tokens=1024,
            temperature=temperature,
            messages=[
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )
        # extract text response
        analysis_text = response.content[0].text

        # parse the response, helper functions scan the text for specific keywords and patterns
        priority = extract_priority_level(analysis_text)
        pathogenicity = extract_pathogenicity_assessment(analysis_text)
        
        return {
            "llm_interpretation": analysis_text,
            "priority_level": priority,
            "pathogenicity_assessment": pathogenicity,
            "experimental_recommendations": extract_experimental_recommendations(analysis_text),
            "model_used": model,
            "prompt_used": prompt
        }
        
    except Exception as e:
        return {
            "llm_interpretation": f"Error calling Claude API: {str(e)}",
            "priority_level": "unknown",
            "pathogenicity_assessment": "unable_to_assess",
            "experimental_recommendations": "manual_review_required",
            "error": str(e)
        }
        
    except Exception as e:
        return {
            "llm_interpretation": f"Error calling LLM API: {str(e)}",
            "priority_level": "unknown",
            "pathogenicity_assessment": "unable_to_assess",
            "experimental_recommendations": "manual_review_required",
            "error": str(e)
        }


# these helpers extract information from claude's natural language response

# scan LLM response for keywords indicating recommend priority for lab follow-up
# medium by default
def extract_priority_level(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    if "high priority" in text_lower or "priority: high" in text_lower:
        return "high"
    elif "medium priority" in text_lower or "priority: medium" in text_lower:
        return "medium"
    elif "low priority" in text_lower or "priority: low" in text_lower:
        return "low"
    else:
        return "medium"

# extract pathogenicity, looks for standard pathogenicity terms following ACMG guideline
def extract_pathogenicity_assessment(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    if "likely pathogenic" in text_lower:
        return "likely_pathogenic"
    elif "possibly pathogenic" in text_lower:
        return "possibly_pathogenic"
    elif "likely benign" in text_lower:
        return "likely_benign"
    elif "uncertain significance" in text_lower:
        return "uncertain_significance"
    else:
        return "uncertain_significance"

# identifies which in vitro testing methods claude recommends
def extract_experimental_recommendations(analysis_text: str) -> str:
    text_lower = analysis_text.lower()
    recommendations = []
    
    if "minigene" in text_lower:
        recommendations.append("minigene_assay")
    if "rt-pcr" in text_lower:
        recommendations.append("rt_pcr")
    if "whole blood" in text_lower:
        recommendations.append("whole_blood_testing")
    return ",".join(recommendations) if recommendations else "manual_review_required"