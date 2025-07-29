#!/usr/bin/env python3
"""
test_alphagenome.py - Testing script for AlphaGenome integration

This script tests the AlphaGenome integration with various variants and scenarios.
Run this to verify that the AlphaGenome API is working correctly with your setup.
"""

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

import json
from typing import Dict, Any, List

# Import the AlphaGenome integration module
try:
    from callAlphaGenome import call_alphagenome, ALPHAGENOME_AVAILABLE
except ImportError:
    print("Error: Cannot import callAlphaGenome module. Make sure it's in the same directory.")
    sys.exit(1)

def test_single_variant(variant: str, hg: str, sequence_length: str = "100KB") -> Dict[str, Any]:
    """Test AlphaGenome with a single variant."""
    print(f"\n🧬 Testing variant: {variant} (GRCh{hg})")
    print(f"   Sequence length: {sequence_length}")
    print("-" * 60)
    
    try:
        results = call_alphagenome(variant, hg, sequence_length)
        
        if "error" not in results:
            print("✅ AlphaGenome Analysis Successful!")
            
            # Summary statistics
            summary = results.get("alphagenome_scores", {}).get("summary_stats", {})
            print(f"\n📊 Summary Statistics:")
            print(f"   Total predictions: {summary.get('total_predictions', 'N/A')}")
            print(f"   Output modalities: {', '.join(summary.get('output_types', []))}")
            print(f"   Unique tissues: {summary.get('unique_tissues', 'N/A')}")
            print(f"   Significant effects: {summary.get('significant_effects', 'N/A')}")
            print(f"   Mean quantile score: {summary.get('mean_quantile_score', 'N/A'):.3f}")
            print(f"   Max absolute effect: {summary.get('max_absolute_effect', 'N/A'):.3f}")
            
            # Splice predictions
            splice_preds = results.get("alphagenome_scores", {}).get("splice_predictions", {})
            if splice_preds and 'message' not in splice_preds:
                print(f"\n🧬 Splice Predictions:")
                for splice_type, data in splice_preds.items():
                    if isinstance(data, dict) and 'total_predictions' in data:
                        print(f"   {splice_type.upper()}:")
                        print(f"     - Total: {data.get('total_predictions', 0)}")
                        print(f"     - Significant: {data.get('significant_predictions', 0)}")
                        if data.get('mean_quantile') is not None:
                            print(f"     - Mean score: {data.get('mean_quantile'):.3f}")
            
            # Expression predictions
            expr_preds = results.get("alphagenome_scores", {}).get("expression_predictions", {})
            if expr_preds and 'message' not in expr_preds:
                print(f"\n📈 Expression Predictions:")
                for expr_type, data in expr_preds.items():
                    if isinstance(data, dict) and 'total_predictions' in data:
                        upregulated = data.get('upregulated_count', 0)
                        downregulated = data.get('downregulated_count', 0)
                        print(f"   {expr_type.upper()}: ↑{upregulated} ↓{downregulated}")
            
            # Top predictions
            top_preds = results.get("top_predictions", [])[:5]
            if top_preds:
                print(f"\n🏆 Top 5 Predicted Effects:")
                for i, pred in enumerate(top_preds, 1):
                    direction = "↑" if pred['effect_direction'] == 'increase' else "↓"
                    print(f"   {i}. {pred['output_type']} in {pred['biosample_name']}")
                    print(f"      {direction} {pred['quantile_score']:.3f} ({pred['significance']})")
            
            return results
            
        else:
            print(f"❌ Error: {results['error']}")
            return results
            
    except Exception as e:
        print(f"❌ Test failed with exception: {e}")
        return {"error": str(e)}

def test_sequence_lengths(variant: str, hg: str) -> None:
    """Test different sequence lengths for the same variant."""
    print(f"\n🔬 Testing Different Sequence Lengths for {variant}")
    print("=" * 70)
    
    sequence_lengths = ["2KB", "16KB", "100KB", "500KB"]
    
    for seq_len in sequence_lengths:
        print(f"\n📏 Testing sequence length: {seq_len}")
        try:
            results = call_alphagenome(variant, hg, seq_len)
            
            if "error" not in results:
                summary = results.get("alphagenome_scores", {}).get("summary_stats", {})
                total_preds = summary.get('total_predictions', 0)
                significant = summary.get('significant_effects', 0)
                max_effect = summary.get('max_absolute_effect', 0)
                
                print(f"   ✅ Success: {total_preds} predictions, {significant} significant, max effect: {max_effect:.3f}")
            else:
                print(f"   ❌ Error: {results['error']}")
                
        except Exception as e:
            print(f"   ❌ Exception: {e}")

def test_multiple_variants() -> None:
    """Test multiple variants to ensure robustness."""
    print(f"\n🧪 Testing Multiple Variants")
    print("=" * 70)
    
    # Test variants (some known to affect splicing)
    test_variants = [
        ("chr22:36201698:A:C", "38", "Known variant in APOB gene region"),
        ("chr17:43124096:G:A", "38", "BRCA1 region variant"),
        ("chr7:117559593:G:A", "38", "CFTR region variant"),
        ("chr11:5247992:G:A", "38", "HBB region variant"),
    ]
    
    results_summary = []
    
    for variant, hg, description in test_variants:
        print(f"\n🔍 Testing: {variant} ({description})")
        try:
            results = call_alphagenome(variant, hg, "100KB")
            
            if "error" not in results:
                summary = results.get("alphagenome_scores", {}).get("summary_stats", {})
                results_summary.append({
                    "variant": variant,
                    "description": description,
                    "status": "success",
                    "total_predictions": summary.get('total_predictions', 0),
                    "significant_effects": summary.get('significant_effects', 0),
                    "max_effect": summary.get('max_absolute_effect', 0)
                })
                print(f"   ✅ Success: {summary.get('total_predictions', 0)} predictions")
            else:
                results_summary.append({
                    "variant": variant,
                    "description": description,
                    "status": "error",
                    "error": results['error']
                })
                print(f"   ❌ Error: {results['error']}")
                
        except Exception as e:
            results_summary.append({
                "variant": variant,
                "description": description,
                "status": "exception",
                "error": str(e)
            })
            print(f"   ❌ Exception: {e}")
    
    # Summary report
    print(f"\n📋 Multi-Variant Test Summary:")
    print("-" * 50)
    successful = sum(1 for r in results_summary if r['status'] == 'success')
    total = len(results_summary)
    print(f"Success rate: {successful}/{total} ({successful/total*100:.1f}%)")
    
    for result in results_summary:
        status_emoji = "✅" if result['status'] == 'success' else "❌"
        print(f"{status_emoji} {result['variant']}: {result['description']}")
        if result['status'] == 'success':
            print(f"   Predictions: {result['total_predictions']}, Significant: {result['significant_effects']}")

def test_error_handling() -> None:
    """Test error handling with invalid inputs."""
    print(f"\n🛡️  Testing Error Handling")
    print("=" * 70)
    
    error_test_cases = [
        ("invalid:format", "38", "Invalid variant format"),
        ("chr1:123:A:T", "99", "Invalid genome build"),
        ("chr99:123456:A:T", "38", "Invalid chromosome"),
        ("chr1:abc:A:T", "38", "Non-numeric position"),
    ]
    
    for variant, hg, description in error_test_cases:
        print(f"\n🔍 Testing: {variant} (GRCh{hg}) - {description}")
        try:
            results = call_alphagenome(variant, hg, "100KB")
            
            if "error" in results:
                print(f"   ✅ Error correctly handled: {results['error']}")
            else:
                print(f"   ⚠️  Unexpected success - error should have been caught")
                
        except Exception as e:
            print(f"   ✅ Exception correctly caught: {e}")

def test_api_connectivity() -> None:
    """Test basic API connectivity and authentication."""
    print(f"\n🌐 Testing API Connectivity")
    print("=" * 70)
    
    if not ALPHAGENOME_AVAILABLE:
        print("❌ AlphaGenome package not available")
        print("   To install: pip install alphagenome")
        return False
    
    # Test with a simple, known variant
    test_variant = "chr22:36201698:A:C"
    print(f"Testing connectivity with variant: {test_variant}")
    
    try:
        # First, test a very basic call with minimal parameters
        print("Step 1: Testing basic AlphaGenome client initialization...")
        
        from alphagenome.models import dna_client
        from alphagenome.data import genome
        
        # Try to create the client
        try:
            API_KEY = "AIzaSyDk6NJZUysBCXV5oywr_8_gtrutcKFwpRg"
            dna_model = dna_client.create(API_KEY)
            print("✅ Client initialized successfully")
        except Exception as client_error:
            print(f"❌ Client initialization failed: {client_error}")
            return False
        
        print("Step 2: Testing variant object creation...")
        
        # Test variant creation
        try:
            chrom, pos, ref, alt = test_variant.split(":")
            variant = genome.Variant(
                chromosome=chrom,
                position=int(pos),
                reference_bases=ref,
                alternate_bases=alt,
            )
            print(f"✅ Variant object created: {variant}")
        except Exception as variant_error:
            print(f"❌ Variant creation failed: {variant_error}")
            return False
        
        print("Step 3: Testing basic API call...")
        
        # Test with a very small sequence and minimal scorers
        try:
            # Use the smallest possible sequence length
            seq_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_2KB']
            interval = variant.reference_interval.resize(seq_length)
            
            print(f"   Interval: {interval}")
            print(f"   Sequence length: 2KB")
            
            # Try with just one basic scorer
            from alphagenome.models import variant_scorers
            
            # Get available scorers
            available_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys())
            print(f"   Available scorers: {available_scorers}")
            
            # Try with the simplest scorer first
            if 'RNA_SEQ' in variant_scorers.RECOMMENDED_VARIANT_SCORERS:
                test_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
                print("   Using RNA_SEQ scorer for test...")
            else:
                # Fallback to first available scorer
                scorer_name = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys())[0]
                test_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[scorer_name]
                print(f"   Using {scorer_name} scorer for test...")
            
            # Make the API call
            variant_scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=[test_scorer]
            )
            
            print("✅ Basic API call successful!")
            
            # Try to process results
            df_scores = variant_scorers.tidy_scores(variant_scores)
            print(f"✅ Results processed: {len(df_scores)} predictions")
            
            return True
            
        except Exception as api_error:
            print(f"❌ API call failed: {api_error}")
            print(f"   Error type: {type(api_error).__name__}")
            import traceback
            print(f"   Traceback: {traceback.format_exc()}")
            return False
            
    except Exception as e:
        print(f"❌ Overall test failed: {e}")
        import traceback
        print(f"   Traceback: {traceback.format_exc()}")
        return False

def run_comprehensive_test() -> None:
    """Run all tests in sequence."""
    print("🚀 AlphaGenome Integration Test Suite")
    print("=" * 70)
    print("This script will test the AlphaGenome integration comprehensively.")
    print("Make sure you have:")
    print("  1. Installed the alphagenome package: pip install alphagenome")
    print("  2. Set your API key: export ALPHA_GENOME_API_KEY='your_key_here'")
    print("  3. Have an internet connection")
    print("=" * 70)
    
    # Test 1: Check package availability
    print(f"\n1️⃣  Package Availability Check")
    if not ALPHAGENOME_AVAILABLE:
        print("❌ AlphaGenome package not available. Install with: pip install alphagenome")
        return
    else:
        print("✅ AlphaGenome package loaded successfully")
    
    # Test 2: API connectivity
    print(f"\n2️⃣  API Connectivity Test")
    if not test_api_connectivity():
        print("❌ API connectivity failed. Check your setup and try again.")
        return
    
    # Test 3: Single variant test
    print(f"\n3️⃣  Single Variant Test")
    test_single_variant("chr22:36201698:A:C", "38", "100KB")
    
    # Test 4: Sequence length tests
    print(f"\n4️⃣  Sequence Length Tests")
    test_sequence_lengths("chr22:36201698:A:C", "38")
    
    # Test 5: Multiple variants
    print(f"\n5️⃣  Multiple Variant Tests")
    test_multiple_variants()
    
    # Test 6: Error handling
    print(f"\n6️⃣  Error Handling Tests")
    test_error_handling()
    
    print(f"\n🎉 Test Suite Complete!")
    print("=" * 70)
    print("If all tests passed, AlphaGenome integration is ready for use in ChatSAV!")

def main():
    """Main function with command line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Test AlphaGenome integration")
    parser.add_argument("--variant", help="Test specific variant (format: chr:pos:ref:alt)")
    parser.add_argument("--hg", choices=["37", "38"], default="38", help="Genome build")
    parser.add_argument("--length", choices=["2KB", "16KB", "100KB", "500KB", "1MB"], 
                       default="100KB", help="Sequence length")
    parser.add_argument("--comprehensive", action="store_true", 
                       help="Run comprehensive test suite")
    
    args = parser.parse_args()
    
    if args.comprehensive:
        run_comprehensive_test()
    elif args.variant:
        test_single_variant(args.variant, args.hg, args.length)
    else:
        # Default: run comprehensive tests
        run_comprehensive_test()

if __name__ == "__main__":
    main()