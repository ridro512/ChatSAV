#!/usr/bin/env python3
"""
chatLLM.py - Interactive chat functionality

Handles all chat-related interactions with the LLM, including:
- Interactive Q&A about variants
- Context management
- Conversation history
- Help and summary functions
- Support for SpliceAI, Pangolin, AlphaGenome, and GTEx data
"""

from typing import List, Dict, Any, Optional
from callLlm import call_llm

class ChatLLM:
    def __init__(self, pipeline_instance):
        self.pipeline = pipeline_instance
        self.conversation_history: List[Dict[str, str]] = []
        self.max_history_length = 5
        self.max_response_tokens = 300
        self.llm_analysis_results = None
    
    def start_chat(self) -> None:
        # start interactive chat session
        print("\n--- Chat with LLM ---")
        
        if not self.pipeline.variant_coord:
            print("Error: Please input variant coordinates first (Option 1)")
            self.pipeline.wait_for_user()
            return

        if not self.pipeline.spliceai_results:
            print("Warning: No SpliceAI results available. Please run SpliceAI analysis first (Option 2) for better chat experience.")
        
        if not self.pipeline.pangolin_results:
            print("Warning: No Pangolin results available. Consider running Pangolin analysis (Option 3) for comparative analysis.")
        
        if not self.pipeline.alphagenome_results:
            print("Warning: No AlphaGenome results available. Consider running AlphaGenome analysis (Option 5) for multi-modal predictions.")
        
        if not self.pipeline.gtex_results:
            print("Warning: No GTEx results available. Please run GTEx analysis first (Option 8) for better chat experience.")
        
        # generate LLM analysis if we have the required data
        if hasattr(self.pipeline, 'last_llm_results') and self.pipeline.last_llm_results:
            print("Using existing LLM analysis results for chat context...")
            self.llm_analysis_results = self.pipeline.last_llm_results
        elif self.pipeline.spliceai_results:  # Minimum requirement
            print("Generating comprehensive analysis for chat context...")
            self.generate_llm_analysis_for_chat()
        
        print(f"Starting interactive chat about variant: {self.pipeline.variant_coord}")
        print("=" * 60)
        
        self.show_available_data()
        self.show_chat_instructions()

        # main chat loop
        while True:
            user_input = input("\nYour question: ").strip()
            
            if user_input.lower() in ['exit', 'quit', 'back']:
                print("Ending chat session...")
                break
            
            if not user_input:
                continue
            
            if self.handle_special_commands(user_input):
                continue

            self.process_chat_question(user_input)
    
    # show what analysis data is available for discussion
    def show_available_data(self) -> None:
        """Show what analysis data is available for discussion."""
        available_data = []
        if self.pipeline.spliceai_results:
            available_data.append("SpliceAI predictions")
        if self.pipeline.pangolin_results:
            available_data.append("Pangolin predictions")
        if self.pipeline.alphagenome_results:
            available_data.append("AlphaGenome multi-modal predictions")
        if self.pipeline.gtex_results:
            available_data.append("GTEx expression data")
        if self.pipeline.context:
            available_data.append("User context")
        if self.llm_analysis_results:
            available_data.append("LLM comparative analysis")
        
        if available_data:
            print(f"Available data for discussion: {', '.join(available_data)}")
        else:
            print("Note: Limited analysis data available. Consider running SpliceAI, Pangolin, AlphaGenome, and GTEx analyses first.")
    
    def generate_llm_analysis_for_chat(self) -> None:
        """Generate LLM analysis results to use as context for chat."""
        try:
            self.llm_analysis_results = call_llm(
                self.pipeline.spliceai_results or {}, 
                self.pipeline.gtex_results or {},
                self.pipeline.pangolin_results,
                self.pipeline.alphagenome_results, 
                self.pipeline.chrom, 
                self.pipeline.pos, 
                self.pipeline.ref, 
                self.pipeline.alt, 
                self.pipeline.context,
                model="gpt-4.1-nano"
            )
            print("✓ Multi-tool analysis generated for chat context")
        except Exception as e:
            print(f"Warning: Could not generate analysis for chat context: {e}")
            self.llm_analysis_results = None
    
    def show_chat_instructions(self) -> None:
        print("\nChat Commands:")
        print("- Type your question and press Enter")
        print("- Type 'exit', 'quit', or 'back' to return to main menu")
        print("- Type 'summary' to get a quick overview of current findings")
        print("- Type 'analysis' to see the full LLM analysis")
        print("- Type 'comparison' to see tool-by-tool comparison")
        print("- Type 'help' to see example questions")
        print("- Type 'clear' to clear conversation history")
        print("-" * 60)
    
    def handle_special_commands(self, user_input: str) -> bool:
        command = user_input.lower()
        
        if command == 'help':
            self.show_chat_help()
            return True
        elif command == 'summary':
            self.show_quick_summary()
            return True
        elif command == 'analysis':
            self.show_full_analysis()
            return True
        elif command == 'comparison':
            self.show_tool_comparison()
            return True
        elif command == 'clear':
            self.conversation_history.clear()
            print("Conversation history cleared.")
            return True
        
        return False
    
    def show_tool_comparison(self) -> None:
        """Show a comparison of predictions across different tools."""
        print("\n--- Tool-by-Tool Comparison ---")
        print("=" * 50)
        
        # SpliceAI results
        if self.pipeline.spliceai_results:
            max_score, max_type, gene_name = self.get_max_spliceai_score()
            print(f"SpliceAI: Highest score = {max_score} ({max_type}) in {gene_name}")
        else:
            print("SpliceAI: No results available")
        
        # Pangolin results
        if self.pipeline.pangolin_results:
            max_score, max_type, gene_name = self.get_max_pangolin_score()
            print(f"Pangolin: Highest score = {max_score} ({max_type}) in {gene_name}")
        else:
            print("Pangolin: No results available")
        
        # AlphaGenome results
        if self.pipeline.alphagenome_results:
            alpha_summary = self.get_alphagenome_summary()
            print(f"AlphaGenome: {alpha_summary}")
        else:
            print("AlphaGenome: No results available")
        
        # GTEx results
        if self.pipeline.gtex_results:
            gtex_summary = self.get_gtex_summary()
            print(f"GTEx: {gtex_summary}")
        else:
            print("GTEx: No results available")
        
        # Tool agreement analysis
        print("\n--- Tool Agreement Analysis ---")
        self.analyze_tool_agreement()
    
    def get_alphagenome_summary(self) -> str:
        """Get a summary of AlphaGenome predictions."""
        if not self.pipeline.alphagenome_results:
            return "No data"
        
        scores = self.pipeline.alphagenome_results.get('alphagenome_scores', {})
        summary_stats = scores.get('summary_stats', {})
        
        if 'error' in summary_stats:
            return "Analysis failed"
        
        total_preds = summary_stats.get('total_predictions', 0)
        sig_effects = summary_stats.get('significant_effects', 0)
        max_effect = summary_stats.get('max_absolute_effect', 0)
        
        return f"{total_preds} predictions, {sig_effects} significant (max effect: {max_effect:.3f})"
    
    def get_gtex_summary(self) -> str:
        """Get a summary of GTEx expression data."""
        if not self.pipeline.gtex_results:
            return "No data"
        
        # actual implementation depends on GTEx result structure
        if isinstance(self.pipeline.gtex_results, dict):
            gene_count = len(self.pipeline.gtex_results)
            tissues = self.pipeline.tissues or ["unknown"]
            return f"{gene_count} genes analyzed in {len(tissues)} tissues"
        
        return "Expression data available"
    
    def analyze_tool_agreement(self) -> None:
        """Analyze agreement between prediction tools."""
        splice_tools = []
        
        if self.pipeline.spliceai_results:
            max_score, _, _ = self.get_max_spliceai_score()
            splice_tools.append(("SpliceAI", max_score))
        
        if self.pipeline.pangolin_results:
            max_score, _, _ = self.get_max_pangolin_score()
            splice_tools.append(("Pangolin", max_score))
        
        if len(splice_tools) >= 2:
            scores = [score for _, score in splice_tools]
            if all(s > 0.5 for s in scores):
                print("✓ Strong agreement: Both tools predict significant splice disruption")
            elif all(s <= 0.2 for s in scores):
                print("✓ Strong agreement: Both tools predict minimal splice impact")
            else:
                print("⚠️ Disagreement: Tools show conflicting predictions")
                for tool, score in splice_tools:
                    print(f"  {tool}: {score}")
        else:
            print("Insufficient tools for agreement analysis")
        
        # AlphaGenome integration
        if self.pipeline.alphagenome_results:
            alpha_scores = self.pipeline.alphagenome_results.get('alphagenome_scores', {})
            splice_preds = alpha_scores.get('splice_predictions', {})
            if splice_preds and 'error' not in splice_preds:
                print("AlphaGenome provides additional multi-modal evidence for validation")
    
    def show_full_analysis(self) -> None:
        """Show the full LLM analysis if available."""
        if self.llm_analysis_results:
            print("\n--- Full LLM Analysis ---")
            print("=" * 50)
            print(f"Priority: {self.llm_analysis_results.get('priority_level', 'N/A')}")
            print(f"Pathogenicity: {self.llm_analysis_results.get('pathogenicity_assessment', 'N/A')}")
            print(f"Recommendations: {self.llm_analysis_results.get('experimental_recommendations', 'N/A')}")
            
            # Show which tools were analyzed
            tools_analyzed = self.llm_analysis_results.get('tools_analyzed', [])
            if tools_analyzed:
                print(f"Tools analyzed: {', '.join(tools_analyzed)}")
            
            print("\nDetailed Analysis:")
            print("-" * 30)
            print(self.llm_analysis_results.get('llm_interpretation', 'No analysis available'))
        else:
            print("No LLM analysis available. Please run the full analysis first (Option 9) or ensure prediction data is available.")
    
    def process_chat_question(self, user_input: str) -> None:
        """Process a regular chat question."""
        try:
            # prepare context for the chat
            chat_context = self.prepare_chat_context(user_input)

            response = self.get_chat_response(chat_context, user_input)
            print(f"\nLLM: {response}")
            
            # add to conversation history
            self.add_to_history(user_input, response)
            
        except Exception as e:
            print(f"Error: {e}")
            print("Please try rephrasing your question or check your API connection.")
    
    def show_chat_help(self) -> None:
        """Show example questions for the chat interface."""
        print("\n--- Example Questions ---")
        print("General variant questions:")
        print("• What type of variant is this?")
        print("• Which gene does this affect?")
        print("• What are the potential consequences?")
        
        if self.pipeline.spliceai_results:
            print("\nSpliceAI-specific questions:")
            print("• Which transcript has the highest disruption score?")
            print("• What type of splicing change is predicted?")
            print("• How confident should I be in these predictions?")
            print("• What does a DS_AL score of X mean?")
        
        if self.pipeline.pangolin_results:
            print("\nPangolin-specific questions:")
            print("• How do Pangolin and SpliceAI predictions compare?")
            print("• What does the Pangolin score tell us?")
            print("• Which model is more reliable for this variant?")
        
        if self.pipeline.alphagenome_results:
            print("\nAlphaGenome-specific questions:")
            print("• What does AlphaGenome predict beyond splicing?")
            print("• Which tissues show the strongest effects?")
            print("• How do expression and chromatin predictions relate?")
            print("• Does AlphaGenome support the splice predictions?")
            print("• What multi-modal evidence is available?")
        
        if self.pipeline.gtex_results:
            print("\nExpression-related questions:")
            print("• Is this gene expressed in the selected tissue?")
            print("• Would expression testing be informative?")
            print("• What tissues should I test?")
        
        if self.llm_analysis_results:
            print("\nComparative analysis questions:")
            print("• Why did you assign this priority level?")
            print("• Which tool provides the strongest evidence?")
            print("• How do the different predictions compare?")
            print("• What makes you recommend these experiments?")
            print("• What additional evidence would change your assessment?")
        
        print("\nClinical questions:")
        print("• What experiments should I prioritize?")
        print("• How would you classify this variant's pathogenicity?")
        print("• What additional evidence would be helpful?")
        print("• Is this variant likely to cause disease?")
        print("• Which prediction tool should I trust most?")
    
    # show a quick summary of current analysis results
    def show_quick_summary(self) -> None:
        print("\n--- Quick Summary ---")
        print(f"Variant: {self.pipeline.variant_coord} ({self.pipeline.chrom}:{self.pipeline.pos} {self.pipeline.ref}>{self.pipeline.alt})")
        
        if self.pipeline.spliceai_results:
            max_score, max_type, gene_name = self.get_max_spliceai_score()
            print(f"SpliceAI: Highest score = {max_score} ({max_type}) in {gene_name}")
        
        if self.pipeline.pangolin_results:
            max_score, max_type, gene_name = self.get_max_pangolin_score()
            print(f"Pangolin: Highest score = {max_score} ({max_type}) in {gene_name}")
        
        if self.pipeline.alphagenome_results:
            alpha_summary = self.get_alphagenome_summary()
            print(f"AlphaGenome: {alpha_summary}")
        
        if self.pipeline.gtex_results:
            gtex_summary = self.get_gtex_summary()
            print(f"GTEx: {gtex_summary}")
        
        if self.llm_analysis_results:
            print(f"LLM Assessment: {self.llm_analysis_results.get('pathogenicity_assessment', 'N/A')}")
            print(f"Priority: {self.llm_analysis_results.get('priority_level', 'N/A')}")
            tools_analyzed = self.llm_analysis_results.get('tools_analyzed', [])
            if tools_analyzed:
                print(f"Tools compared: {', '.join(tools_analyzed)}")
        
        if self.pipeline.context:
            context_preview = self.pipeline.context[:100] + ('...' if len(self.pipeline.context) > 100 else '')
            print(f"User context: {context_preview}")
        
        print(f"Conversation history: {len(self.conversation_history)} exchanges")
    
    def get_max_spliceai_score(self) -> tuple:
        max_score = 0
        max_type = "None"
        gene_name = "Unknown"
        
        for transcript in self.pipeline.spliceai_results.get('splice_scores', []):
            current_gene = transcript.get('g_name', transcript.get('ensembl_id', 'Unknown'))
            for score_type in ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']:
                try:
                    score = float(transcript.get(score_type, 0))
                    if score > max_score:
                        max_score = score
                        max_type = score_type
                        gene_name = current_gene
                except (ValueError, TypeError):
                    continue
        
        return max_score, max_type, gene_name
    
    def get_max_pangolin_score(self) -> tuple:
        max_score = 0
        max_type = "None"
        gene_name = "Unknown"
        
        for transcript in self.pipeline.pangolin_results.get('pangolin_scores', []):
            current_gene = transcript.get('g_name', transcript.get('ensembl_id', 'Unknown'))
            for score_type in ['DS_SG', 'DS_SL']:
                score_val = transcript.get(score_type)
                if score_val is not None:
                    try:
                        score = abs(float(score_val))
                        if score > max_score:
                            max_score = score
                            max_type = score_type
                            gene_name = current_gene
                    except (ValueError, TypeError):
                        continue
        
        return max_score, max_type, gene_name
    
    def prepare_chat_context(self, user_question: str) -> str:
        context_parts = []
        
        # add variant information
        context_parts.append("VARIANT INFORMATION:")
        context_parts.append(f"- Coordinates: {self.pipeline.variant_coord}")
        context_parts.append(f"- Chromosome: {self.pipeline.chrom}, Position: {self.pipeline.pos}")
        context_parts.append(f"- Change: {self.pipeline.ref} → {self.pipeline.alt}")
        context_parts.append(f"- Genome build: GRCh{self.pipeline.hg}")

        # add LLM analysis results if available (priority)
        if self.llm_analysis_results:
            context_parts.append("\nPREVIOUS LLM ANALYSIS:")
            context_parts.append(f"- Priority Level: {self.llm_analysis_results.get('priority_level')}")
            context_parts.append(f"- Pathogenicity Assessment: {self.llm_analysis_results.get('pathogenicity_assessment')}")
            context_parts.append(f"- Experimental Recommendations: {self.llm_analysis_results.get('experimental_recommendations')}")
            tools_analyzed = self.llm_analysis_results.get('tools_analyzed', [])
            if tools_analyzed:
                context_parts.append(f"- Tools Analyzed: {', '.join(tools_analyzed)}")
            context_parts.append(f"- Full Analysis: {self.llm_analysis_results.get('llm_interpretation', '')[:500]}...")
        
        # add SpliceAI results if available
        if self.pipeline.spliceai_results:
            context_parts.append("\nSPLICEAI RESULTS:")
            context_parts.append(f"- Distance parameter: {self.pipeline.distance}")
            context_parts.append(f"- Mask parameter: {self.pipeline.mask}")
            
            for transcript in self.pipeline.spliceai_results.get('splice_scores', []):
                gene_info = transcript.get('g_name', transcript.get('ensembl_id', 'Unknown'))
                context_parts.append(
                    f"- {gene_info} ({transcript.get('transcript_id')}): "
                    f"AG={transcript.get('DS_AG')}, AL={transcript.get('DS_AL')}, "
                    f"DG={transcript.get('DS_DG')}, DL={transcript.get('DS_DL')}"
                )
        
        # add Pangolin results if available
        if self.pipeline.pangolin_results:
            context_parts.append("\nPANGOLIN RESULTS:")
            for transcript in self.pipeline.pangolin_results.get('pangolin_scores', []):
                gene_info = transcript.get('g_name', transcript.get('ensembl_id', 'Unknown'))
                context_parts.append(
                    f"- {gene_info} ({transcript.get('transcript_id')}): "
                    f"SG={transcript.get('DS_SG')}, SL={transcript.get('DS_SL')}"
                )
        
        # add AlphaGenome results if available
        if self.pipeline.alphagenome_results:
            context_parts.append("\nALPHAGENOME RESULTS:")
            scores = self.pipeline.alphagenome_results.get('alphagenome_scores', {})
            
            # Summary stats
            summary = scores.get('summary_stats', {})
            if summary and 'error' not in summary:
                context_parts.append(
                    f"- Summary: {summary.get('total_predictions', 0)} predictions, "
                    f"{summary.get('significant_effects', 0)} significant effects"
                )
            
            # Top predictions
            top_preds = self.pipeline.alphagenome_results.get('top_predictions', [])
            if top_preds:
                context_parts.append("- Top effects:")
                for pred in top_preds[:3]:  # Top 3
                    context_parts.append(
                        f"  * {pred.get('output_type')} in {pred.get('biosample_name')}: "
                        f"{pred.get('quantile_score', 0):.3f} ({pred.get('significance')})"
                    )
        
        # add GTEx results if available
        if self.pipeline.gtex_results:
            context_parts.append("\nGTEX EXPRESSION:")
            # Simplified representation - adapt based on actual GTEx result structure
            context_parts.append(f"- Tissues analyzed: {', '.join(self.pipeline.tissues or ['Unknown'])}")
            context_parts.append(f"- Genes analyzed: {len(self.pipeline.gtex_results) if isinstance(self.pipeline.gtex_results, dict) else 1}")
        
        # add user context if available
        if self.pipeline.context:
            context_parts.append("\nUSER CONTEXT:")
            context_parts.append(self.pipeline.context)
        
        # add recent conversation history
        if self.conversation_history:
            context_parts.append("\nRECENT CONVERSATION:")
            for exchange in self.conversation_history[-3:]:
                context_parts.append(f"User: {exchange['user']}")
                assistant_preview = exchange['assistant'][:200] + ('...' if len(exchange['assistant']) > 200 else '')
                context_parts.append(f"Assistant: {assistant_preview}")
        
        return "\n".join(context_parts)
    
    def get_chat_response(self, context: str, user_question: str) -> str:
        """Get response from LLM for chat interaction."""
        try:
            from callLlm import client
            if not client:
                return "OpenAI client not available. Please check your API key configuration."
        except ImportError:
            return "OpenAI client not properly configured. Please check callLlm.py setup."
        
        try:
            # Determine available tools for system prompt
            available_tools = []
            if self.pipeline.spliceai_results:
                available_tools.append("SpliceAI")
            if self.pipeline.pangolin_results:
                available_tools.append("Pangolin")
            if self.pipeline.alphagenome_results:
                available_tools.append("AlphaGenome")
            if self.pipeline.gtex_results:
                available_tools.append("GTEx")
            
            system_prompt = f"""You are an expert genetic analyst having a conversation about a specific genetic variant.

You have access to comprehensive analysis results from multiple prediction tools: {', '.join(available_tools) if available_tools else 'SpliceAI only'}.

Use this information to provide informed, consistent responses.

CRITICAL: Maintain strict consistency with the previous analysis, especially:
- Use the EXACT same priority level from the previous analysis
- Use the EXACT same pathogenicity assessment from the previous analysis  
- When discussing experimental recommendations, focus on the MOST IMPORTANT recommendation from the previous analysis
- Do not suggest new experiments beyond what was already recommended in the formal analysis

Provide conversational, helpful responses that:
- Reference and build upon the previous LLM analysis when relevant
- Compare predictions across different tools when appropriate
- Answer the user's specific question clearly and accurately
- Explain your reasoning based on the available multi-tool data
- Maintain strict consistency with previous assessments 
- Explain complex genetic concepts in accessible terms
- Maintain a helpful, professional tone

For splice prediction scores:
- SpliceAI/Pangolin scores >0.5 are considered significant, >0.8 are high confidence
- AlphaGenome quantile scores >0.5 indicate significant tissue effects
- When tools disagree, explain which prediction is more reliable and why

When referencing previous analysis, be specific about what led to those conclusions and maintain the same conclusions."""
            
            full_prompt = f"""{context}

USER QUESTION: {user_question}

Please provide a helpful response based on the available multi-tool information about this variant."""
            
            # call OpenAI
            response = client.chat.completions.create(
                model="gpt-4.1-nano",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": full_prompt}
                ],
                max_tokens=self.max_response_tokens,
                temperature=0.3
            )
            
            return response.choices[0].message.content.strip()
            
        except Exception as e:
            return f"Error getting response from LLM: {str(e)}"
    
    def add_to_history(self, user_input: str, response: str) -> None:
        self.conversation_history.append({
            "user": user_input,
            "assistant": response
        })
        
        # keep only recent history to avoid token limits
        if len(self.conversation_history) > self.max_history_length:
            self.conversation_history = self.conversation_history[-self.max_history_length:]
    
    # get summary of current convo for external use
    def get_conversation_summary(self) -> str:
        if not self.conversation_history:
            return "No conversation history available."
        
        summary_parts = [f"Conversation about variant {self.pipeline.variant_coord}:"]
        for i, exchange in enumerate(self.conversation_history, 1):
            summary_parts.append(f"{i}. Q: {exchange['user'][:100]}...")
            summary_parts.append(f"   A: {exchange['assistant'][:150]}...")
        
        return "\n".join(summary_parts)
    
    def clear_history(self) -> None:
        self.conversation_history.clear()