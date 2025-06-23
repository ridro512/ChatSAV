#!/usr/bin/env python3

from callLlm import call_llm
from callSplice import call_splice

spliceai_data = call_splice("chr8:140300616:T:G")

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