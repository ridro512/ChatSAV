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

from typing import Tuple, Dict, Any

from callSplice import call_splice
# from callGtex import call_gtex  # results of GTEx whole blood expression
from callLlm import call_llm


def parse_variant_coordinates(variant_coord: str) -> Tuple[str, int, str, str]:
    """
    Parse variant coordinates from string format.
    
    Args:
        variant_coord: Variant in format "chr:pos:ref:alt"
        
    Returns:
        Tuple of (chromosome, position, reference_allele, alternate_allele)
    """
    chrom, pos, ref, alt = variant_coord.split(":")
    return chrom, int(pos), ref, alt


def print_results(llm_results: Dict[str, Any], chrom: str, pos: int, ref: str, alt: str) -> None:
    """
    Print formatted analysis results.
    
    Args:
        llm_results: Dictionary containing LLM analysis results
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
    """
    print(f"\nResult of a mutation of {ref} > {alt} at position {pos} on chromosome {chrom}")
    print("=" * 60)
    
    print(f"Priority: {llm_results['priority_level']}")
    print(f"Assessment: {llm_results['pathogenicity_assessment']}")
    print(f"Recommendations: {llm_results['experimental_recommendations']}")
    
    if 'error' in llm_results:
        print(f"\nError: {llm_results['error']}")
    else:
        print("\nFull Analysis:")
        print("-" * 40)
        print(llm_results['llm_interpretation'])


def main() -> None:
    """Main function to run the ChatSAV analysis."""
    print("Welcome to ChatSAV")
    print("This assistant analyses splice-altering variants (SAVs) using SpliceAI and GTEx data.")
    
    # Get user input
    variant_coord = input(
        "Please enter the variant coordinates in the format chr:pos:ref:alt "
        "(e.g., chr1:123456:A:T):\n"
    )
    hg = input(
        "Please enter the genome build [37 for GRCh37/hg19] or [38 for GRCh38/hg38]:\n"
    )
    tissue = input("Please enter the tissue type (e.g., Whole_Blood):\n")
    
    # Parse variant coordinates
    chrom, pos, ref, alt = parse_variant_coordinates(variant_coord)
    
    # Call analysis functions
    spliceai_results = call_splice(variant_coord, hg)
    
    # TODO: Replace with actual GTEx call when implemented
    gtex_results = {
        "whole_blood_tpm": 0.5,
        "expressed_in_whole_blood": False
    }
    
    llm_results = call_llm(
        spliceai_results, 
        gtex_results, 
        chrom, 
        pos, 
        ref, 
        alt, 
        model="claude-3-5-sonnet-20241022"
    )
    
    # Display results
    print_results(llm_results, chrom, pos, ref, alt)
    
if __name__ == "__main__":
    main()
