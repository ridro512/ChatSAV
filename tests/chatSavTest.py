#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

from callLlm import call_llm
from callSplice import call_splice

spliceai_data = call_splice("chr8:140300616:T:G", 38)

gtex_data = {
    "whole_blood_tpm": 0.5,
    "expressed_in_whole_blood": False
}

result = call_llm(
    spliceai_data, 
    gtex_data, 
    "chr8", 
    140300616, 
    "T", 
    "G",
    model="claude-3-5-sonnet-20241022"
)


print("Result:")
print(f"Priority: {result['priority_level']}")
print(f"Assessment: {result['pathogenicity_assessment']}")
print(f"Recommendations: {result['experimental_recommendations']}")

if 'error' in result:
    print(f"\nError: {result['error']}")
else:
    print("\nSuccess! Full analysis:")
    print("=" * 50)
    print(result['llm_interpretation'])