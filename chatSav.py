#!/usr/bin/env python3
"""
chatSav.py - Main file for LLM-powered SAV detecting assistant

Input: Coordinates of a variant (chr:pos:ref:alt format)
Output: SpliceAI predictions and whole blood expression status

VARIANT FORMAT INFORMATION:
- Chromosome: chr
- Position: pos
- Reference allele: ref
- Alternate allele: alt
"""

from callSplice import call_splice  # results of SpliceAI analysis
#from callGtex import call_gtex  # results of GTEx whole blood expression
from callLlm import call_llm  # results of LLM analysis

def parse_variant_coordinates(variant_coord: str):
    chrom, pos, ref, alt = variant_coord.split(":")
    return chrom, int(pos), ref, alt

def main():
    print("Welcome to ChatSAV")
    print("This tool analyzes splice-altering variants (SAVs) using SpliceAI and GTEx data.")
    # will change to plain text parsed through LLM into these coordinates
    variant_coord = input("Please enter the variant coordinates in the format chr:pos:ref:alt (e.g., chr1:123456:A:T):\n")
    
    # Parse the variant coordinates
    chrom, pos, ref, alt = parse_variant_coordinates(variant_coord)
    
    # Call SpliceAI, GTEx, and LLM functions (to be implemented)
    spliceai_results = call_splice(variant_coord)
    #gtex_results = call_gtex(variant_coord)
    # THIS IS TEST THIS IS TEMPORARY
    gtex_results = {
        "whole_blood_tpm": 0.5,  # Placeholder value
        "expressed_in_whole_blood": False  # Placeholder value
    }
    llm_results = call_llm(spliceai_results, gtex_results, chrom, pos, ref, alt, model="claude-3-5-sonnet-20241022")
    
    print("Result of a mutation of", ref, ">", alt, "at position", pos, "on chromosome", chrom)
    print(f"Priority: {llm_results['priority_level']}")
    print(f"Assessment: {llm_results['pathogenicity_assessment']}")
    print(f"Recommendations: {llm_results['experimental_recommendations']}")
    if 'error' in llm_results:
        print(f"\nError: {llm_results['error']}")
    else:
        print("\nSuccess! Full analysis:")
        print("=" * 50)
        print(llm_results['llm_interpretation'])

if __name__ == "__main__":
    main()