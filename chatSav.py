#!/usr/bin/env python3
"""
chatSav.py - Main file for LLM-powered SAV detecting assistant

Input: Coordinates of a variant (chr:pos:ref:alt format)
Output: SpliceAI predictions and whole blood expression status
"""
"""
from callSplice import #results of spliceAI
from callGtex import #results of GTEx whole blood expression
from callLlm import #results of LLM analysis
"""
def parse_variant_coordinates(variant_coord: str):
    chrom, pos, ref, alt = variant_coord.split(":")
    return chrom, int(pos), ref, alt

def main():
    # Example variant coordinate
    variant_coord = "chr1:123456:A:T"
    
    # Parse the variant coordinates
    chrom, pos, ref, alt = parse_variant_coordinates(variant_coord)
    
    # Call SpliceAI, GTEx, and LLM functions (to be implemented)
    spliceai_results = call_splice(chrom, pos, ref, alt)
    gtex_results = call_gtex(chrom, pos, ref, alt)
    llm_results = call_llm(chrom, pos, ref, alt)
    
    # Print or process results
    print("SpliceAI Results:", spliceai_results)
    print("GTEx Results:", gtex_results)
    print("LLM Results:", llm_results)