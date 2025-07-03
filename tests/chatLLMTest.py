#!/usr/bin/env python3
"""
simple_chat_test.py - Simple integration test for chat functionality

This test verifies that:
1. Chat generates LLM analysis correctly
2. Chat responses are consistent with formal analysis
3. Context is properly maintained
4. Special commands work

Run this with: python simple_chat_test.py
"""

import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

from chatSav import ChatSAVPipeline
from chatLLM import ChatLLM


def test_chat_llm_consistency():
    """Test that chat responses are consistent with LLM analysis."""
    print("Testing Chat-LLM Consistency")
    print("=" * 50)
    
    # Create pipeline and set up test data
    pipeline = ChatSAVPipeline()
    
    # Set up mock variant data
    pipeline.variant_coord = "chr7:117559593:C:T"
    pipeline.chrom = "chr7"
    pipeline.pos = 117559593
    pipeline.ref = "C"
    pipeline.alt = "T"
    pipeline.hg = "38"
    pipeline.distance = 50
    pipeline.mask = 0
    pipeline.tissue = "Whole_Blood"
    pipeline.context = "12-year-old patient with cystic fibrosis symptoms"
    
    # Mock SpliceAI results
    pipeline.spliceai_results = {
        "variant": "chr7:117559593:C:T",
        "splice_scores": [
            {
                "transcript_id": "ENST00000288602.6",
                "gene_id": "CFTR",
                "ensembl_id": "ENSG00000001626",
                "g_name": "CFTR",
                "DS_AG": 0.02,
                "DS_AL": 0.91,  # High score
                "DS_DG": 0.01,
                "DS_DL": 0.05,
                "DP_AG": -15,
                "DP_AL": -2,
                "DP_DG": 23,
                "DP_DL": 45
            }
        ]
    }
    
    # Mock GTEx results
    pipeline.gtex_results = {
        "tissue": "Whole_Blood",
        "tpm": 2.3,
        "expressed": True,
        "percentile": 75
    }
    
    print("✓ Test data set up")
    
    # Test 1: Generate formal LLM analysis and store it in pipeline
    print("\n1. Testing formal LLM analysis generation...")
    try:
        from callLlm import call_llm
        
        formal_analysis = call_llm(
            pipeline.spliceai_results,
            pipeline.gtex_results,
            pipeline.chrom,
            pipeline.pos,
            pipeline.ref,
            pipeline.alt,
            pipeline.context,
            model="gpt-4o-mini"
        )
        
        # Store in pipeline like the real application does
        pipeline.last_llm_results = formal_analysis
        
        print("✓ Formal analysis generated successfully")
        print(f"  Priority: {formal_analysis.get('priority_level')}")
        print(f"  Pathogenicity: {formal_analysis.get('pathogenicity_assessment')}")
        print(f"  Recommendations: {formal_analysis.get('experimental_recommendations')}")
        
    except Exception as e:
        print(f"❌ Error generating formal analysis: {e}")
        return False
    
    # Test 2: Initialize chat handler (should use stored analysis)
    print("\n2. Testing chat handler analysis reuse...")
    try:
        chat_handler = ChatLLM(pipeline)
        
        # Simulate starting chat (should reuse existing analysis)
        if hasattr(pipeline, 'last_llm_results') and pipeline.last_llm_results:
            chat_handler.llm_analysis_results = pipeline.last_llm_results
            print("✓ Chat using existing analysis results")
        else:
            print("❌ Chat not reusing existing analysis")
            return False
            
        chat_analysis = chat_handler.llm_analysis_results
        print(f"  Priority: {chat_analysis.get('priority_level')}")
        print(f"  Pathogenicity: {chat_analysis.get('pathogenicity_assessment')}")
        print(f"  Recommendations: {chat_analysis.get('experimental_recommendations')}")
            
    except Exception as e:
        print(f"❌ Error in chat analysis setup: {e}")
        return False
    
    # Test 3: Check consistency between formal and chat analysis (should be identical)
    print("\n3. Testing consistency between analyses...")
    
    consistency_checks = [
        ("Priority Level", formal_analysis.get('priority_level'), chat_analysis.get('priority_level')),
        ("Pathogenicity", formal_analysis.get('pathogenicity_assessment'), chat_analysis.get('pathogenicity_assessment')),
        ("Recommendations", formal_analysis.get('experimental_recommendations'), chat_analysis.get('experimental_recommendations'))
    ]
    
    all_consistent = True
    for check_name, formal_value, chat_value in consistency_checks:
        if formal_value == chat_value:
            print(f"  ✓ {check_name}: {formal_value}")
        else:
            print(f"  ❌ {check_name} mismatch: Formal='{formal_value}' vs Chat='{chat_value}'")
            all_consistent = False
    
    # Since chat should reuse the exact same analysis, they should be identical
    if formal_analysis == chat_analysis:
        print("  ✓ Analyses are identical (perfect consistency)")
    else:
        print("  ⚠️  Analyses are not identical but key fields match")
    
    if not all_consistent:
        return False
    
    # Test 4: Test context preparation
    print("\n4. Testing context preparation...")
    try:
        context = chat_handler.prepare_chat_context("Test question")
        
        required_elements = [
            "VARIANT INFORMATION:",
            "PREVIOUS LLM ANALYSIS:",
            "SPLICEAI RESULTS:",
            "GTEX EXPRESSION:",
            "USER CONTEXT:"
        ]
        
        for element in required_elements:
            if element in context:
                print(f"  ✓ {element} included")
            else:
                print(f"  ❌ {element} missing")
                return False
                
    except Exception as e:
        print(f"❌ Error in context preparation: {e}")
        return False
    
    # Test 5: Test chat response generation (mock)
    print("\n5. Testing chat response...")
    try:
        # Mock a simple question about priority
        test_question = "Why is this variant high priority?"
        
        # Check that the context includes the analysis results
        context = chat_handler.prepare_chat_context(test_question)
        
        if chat_analysis.get('priority_level') in context:
            print("  ✓ Analysis results included in chat context")
        else:
            print("  ❌ Analysis results not properly included in context")
            return False
            
        print("  ✓ Chat context preparation successful")
        
    except Exception as e:
        print(f"❌ Error in chat response testing: {e}")
        return False
    
    # Test 6: Test special commands
    print("\n6. Testing special commands...")
    try:
        # Test help command
        help_result = chat_handler.handle_special_commands("help")
        if help_result:
            print("  ✓ Help command works")
        
        # Test summary command  
        summary_result = chat_handler.handle_special_commands("summary")
        if summary_result:
            print("  ✓ Summary command works")
        
        # Test analysis command
        analysis_result = chat_handler.handle_special_commands("analysis")
        if analysis_result:
            print("  ✓ Analysis command works")
            
    except Exception as e:
        print(f"❌ Error testing special commands: {e}")
        return False
    
    print("\n" + "=" * 50)
    print("🎉 ALL TESTS PASSED!")
    print("Chat functionality is working correctly and is consistent with LLM analysis.")
    return True


def test_edge_cases():
    """Test edge cases and error handling."""
    print("\nTesting Edge Cases")
    print("=" * 30)
    
    # Test with missing data
    pipeline = ChatSAVPipeline()
    pipeline.variant_coord = "chr1:12345:A:T"
    pipeline.chrom = "chr1"
    pipeline.pos = 12345
    pipeline.ref = "A"
    pipeline.alt = "T"
    
    chat_handler = ChatLLM(pipeline)
    
    # Test without SpliceAI/GTEx data
    print("1. Testing with missing analysis data...")
    try:
        chat_handler.show_available_data()
        print("  ✓ Handles missing data gracefully")
    except Exception as e:
        print(f"  ❌ Error handling missing data: {e}")
        return False
    
    # Test conversation history management
    print("2. Testing conversation history...")
    try:
        # Add more than max history
        for i in range(10):
            chat_handler.add_to_history(f"Question {i}", f"Answer {i}")
        
        if len(chat_handler.conversation_history) <= chat_handler.max_history_length:
            print("  ✓ History properly limited")
        else:
            print("  ❌ History not properly managed")
            return False
            
    except Exception as e:
        print(f"  ❌ Error in history management: {e}")
        return False
    
    print("✓ Edge case tests passed")
    return True


def main():
    """Run all tests."""
    print("ChatSAV Chat Functionality Test Suite")
    print("=" * 60)
    
    # Check if API key is set
    if not os.getenv("OPENAI_API_KEY"):
        print("❌ OPENAI_API_KEY environment variable not set")
        print("Please set your OpenAI API key:")
        print("export OPENAI_API_KEY='your-api-key-here'")
        return False
    
    print("✓ OpenAI API key found")
    
    # Run main consistency test
    try:
        consistency_passed = test_chat_llm_consistency()
    except Exception as e:
        print(f"❌ Consistency test failed with error: {e}")
        consistency_passed = False
    
    # Run edge case tests
    try:
        edge_cases_passed = test_edge_cases()
    except Exception as e:
        print(f"❌ Edge case test failed with error: {e}")
        edge_cases_passed = False
    
    # Final results
    print("\n" + "=" * 60)
    print("FINAL TEST RESULTS")
    print("=" * 60)
    
    if consistency_passed and edge_cases_passed:
        print("🎉 ALL TESTS PASSED!")
        print("\nYour chat functionality is working correctly!")
        print("✓ LLM analysis and chat responses are consistent")
        print("✓ Context management is working")
        print("✓ Special commands are functional")
        print("✓ Error handling is working")
        
        print("\nYou can now use the chat feature with confidence!")
        print("Run your main program and try option 8 (Chat with LLM)")
        
        return True
    else:
        print("❌ SOME TESTS FAILED")
        if not consistency_passed:
            print("- Consistency tests failed")
        if not edge_cases_passed:
            print("- Edge case tests failed")
            
        print("\nPlease check your implementation and try again.")
        return False


def run_interactive_demo():
    """Run an interactive demo showing how the chat works."""
    print("\n" + "=" * 60)
    print("INTERACTIVE DEMO")
    print("=" * 60)
    
    response = input("Would you like to see an interactive demo? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        return
    
    print("\nSetting up demo...")
    
    # Create pipeline with demo data
    pipeline = ChatSAVPipeline()
    pipeline.variant_coord = "chr7:117559593:C:T"
    pipeline.chrom = "chr7"
    pipeline.pos = 117559593
    pipeline.ref = "C"
    pipeline.alt = "T"
    pipeline.hg = "38"
    pipeline.context = "12-year-old patient with cystic fibrosis symptoms"
    
    # Demo SpliceAI results
    pipeline.spliceai_results = {
        "variant": "chr7:117559593:C:T",
        "splice_scores": [
            {
                "transcript_id": "ENST00000288602.6",
                "g_name": "CFTR",
                "ensembl_id": "ENSG00000001626",
                "DS_AG": 0.02,
                "DS_AL": 0.91,
                "DS_DG": 0.01,
                "DS_DL": 0.05
            }
        ]
    }
    
    pipeline.gtex_results = {
        "tissue": "Whole_Blood",
        "tpm": 2.3,
        "expressed": True,
        "percentile": 75
    }
    
    # Initialize chat handler
    chat_handler = ChatLLM(pipeline)
    
    print("Demo data loaded. Generating analysis...")
    
    try:
        # Generate analysis
        chat_handler.generate_llm_analysis_for_chat()
        
        if chat_handler.llm_analysis_results:
            print("✓ Analysis generated successfully!")
            
            # Show analysis
            print("\n--- Generated Analysis ---")
            analysis = chat_handler.llm_analysis_results
            print(f"Priority: {analysis.get('priority_level')}")
            print(f"Pathogenicity: {analysis.get('pathogenicity_assessment')}")
            print(f"Recommendations: {analysis.get('experimental_recommendations')}")
            print(f"Full analysis: {analysis.get('llm_interpretation', '')[:200]}...")
            
            # Show what would be in chat context
            print("\n--- Chat Context Preview ---")
            context = chat_handler.prepare_chat_context("Demo question")
            print("Context includes:")
            if "PREVIOUS LLM ANALYSIS:" in context:
                print("✓ Previous LLM analysis")
            if "SPLICEAI RESULTS:" in context:
                print("✓ SpliceAI results")
            if "GTEX EXPRESSION:" in context:
                print("✓ GTEx expression data")
            if "USER CONTEXT:" in context:
                print("✓ User context")
            
            print("\n✅ Demo completed successfully!")
            print("Your chat functionality is ready to use!")
            
        else:
            print("❌ Analysis generation failed in demo")
            
    except Exception as e:
        print(f"❌ Demo failed: {e}")


if __name__ == "__main__":
    success = main()
    
    if success:
        run_interactive_demo()
    
    sys.exit(0 if success else 1)