#!/usr/bin/env python3
"""
callGtexTest.py - Test script for calling GTEx API for gene expression data
"""
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

from callGtex import call_gtex


def test_gtex_api():
    """Test the GTEx API with known genes and tissues"""
    
    print("Testing GTEx API...")
    print("=" * 50)
    
    # Test cases: (gene_id, tissue, expected_result_description)
    test_cases = [
        ("ENSG00000141510.16", "Whole_Blood", "TP53 - should have expression")
    ]
    
    for gene_id, tissue, description in test_cases:
        print(f"\nTest: {description}")
        print(f"Query: {gene_id} in {tissue}")
        
        result = call_gtex(gene_id, tissue)
        
        if result is not None:
            print(f"✅ Result: {result} TPM")
            if result > 1.0:
                print("   High expression detected")
            elif result > 0.1:
                print("   Moderate expression detected")
            else:
                print("   Low expression detected")
        else:
            print("❌ Result: No data found")
        
        print("-" * 30)


def test_specific_gene():
    """Test a specific gene across multiple tissues"""
    
    gene_id = "ENSG00000141510.16"  # TP53
    tissues = ["Whole_Blood", "Brain_Cortex", "Heart_Left_Ventricle", "Liver", "Lung"]
    
    print(f"\nTesting {gene_id} (TP53) across multiple tissues:")
    print("=" * 50)
    
    for tissue in tissues:
        result = call_gtex(gene_id, tissue)
        if result is not None:
            print(f"{tissue:20}: {result:6.2f} TPM")
        else:
            print(f"{tissue:20}: No data")


if __name__ == "__main__":
    test_gtex_api()
    test_specific_gene()