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

from typing import Tuple, Dict, Any, Optional

from callSplice import call_splice
from callPangolin import call_pangolin
from callGtex import (
    medianGeneExpression, 
    geneExpression, 
    medianExonExpression, 
    medianJunctionExpression,
    topExpressedGenes,
    call_gtex  # Keep for backward compatibility
)
from callLlm import call_llm
from chatLLM import ChatLLM


class ChatSAVPipeline:
    """Main pipeline class for ChatSAV analysis."""
    
    def __init__(self):
        self.variant_coord: Optional[str] = None
        self.hg: Optional[str] = None
        self.distance: int = 50
        self.mask: int = 0
        self.tissues: Optional[List[str]] = []
        self.spliceai_results: Optional[Dict[str, Any]] = None
        self.pangolin_results: Optional[Dict[str, Any]] = None
        self.gtex_results: Optional[Dict[str, Any]] = None
        self.context: Optional[str] = None
        self.chrom: Optional[str] = None
        self.pos: Optional[int] = None
        self.ref: Optional[str] = None
        self.alt: Optional[str] = None
        self.gtex_endpoint: str = "medianGeneExpression"
        self.chat_handler = ChatLLM(self)

    def wait_for_user(self) -> None:
        """Wait for user to press a key before continuing."""
        input("\nPress Enter to return to main menu...")

    def parse_variant_coordinates(self, variant_coord: str) -> Tuple[str, int, str, str]:
        """
        Parse variant coordinates from string format.
        
        Args:
            variant_coord: Variant in format "chr:pos:ref:alt"
            
        Returns:
            Tuple of (chromosome, position, reference_allele, alternate_allele)
        """
        try:
            chrom, pos, ref, alt = variant_coord.split(":")
            return chrom, int(pos), ref, alt
        except ValueError as e:
            raise ValueError(f"Invalid variant format. Expected 'chr:pos:ref:alt', got '{variant_coord}'") from e

    def input_variant_coordinates(self) -> None:
        """Get variant coordinates and genome build from user."""
        print("\n--- Input Variant Coordinates ---")
        
        while True:
            variant_coord = input(
                "Please enter the variant coordinates in the format chr:pos:ref:alt "
                "(e.g., chr1:123456:A:T):\n"
            ).strip()
            
            try:
                chrom, pos, ref, alt = self.parse_variant_coordinates(variant_coord)
                self.variant_coord = variant_coord
                self.chrom, self.pos, self.ref, self.alt = chrom, pos, ref, alt
                print(f"✓ Variant coordinates set: {variant_coord}")
                break
            except ValueError as e:
                print(f"Error: {e}")
                print("Please try again.")
        
        while True:
            hg = input(
                "Please enter the genome build [37 for GRCh37/hg19] or [38 for GRCh38/hg38]:\n"
            ).strip()
            
            if hg in ["37", "38"]:
                self.hg = hg
                print(f"✓ Genome build set: GRCh{hg}")
                break
            else:
                print("Error: Please enter either '37' or '38'")
            
        while True:
            dist = input(
                "(Optional) Please enter the distance parameter of SpliceAI model (default: 50)\n"
            ).strip()

            try:
                distance = int(dist)
                self.distance = distance
                print(f"✓ SpliceAI distance set: {distance}")
                break
            except ValueError:
                print(f"{dist} is not a valid integer, defaulting to 50\n")
                break
            
        while True:
            mask = input(
                "(Optional) Please enter if you would like raw (0) or masked scores (1) (default: 0)\n"
            ).strip()

            try:
                mask = int(mask)
                if mask in [0,1]:
                    self.mask = mask
                    print(f"✓ SpliceAI mask set: {mask}")
                    break
                else:
                    print(f"{mask} is not a valid option, defaulting to raw scores 0\n")
                    break
            except ValueError:
                print(f"{mask} is not a valid option, defaulting to raw scores 0\n")
                break

        self.wait_for_user()

    def get_spliceai_results(self) -> None:
        """Get SpliceAI results for the current variant."""
        print("\n--- Getting SpliceAI Results ---")
        
        if not self.variant_coord or not self.hg:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        print(f"Analyzing variant: {self.variant_coord} (GRCh{self.hg})")
        print(f"Distance parameter: {self.distance}")
        print(f"Mask parameter: {self.mask} (0 = raw score, 1 = masked score)")
        
        try:
            self.spliceai_results = call_splice(self.variant_coord, self.hg, self.distance, self.mask)
            
            if 'error' in self.spliceai_results:
                print(f"Error calling SpliceAI: {self.spliceai_results['error']}")
            else:
                print("\n✓ SpliceAI Results:")
                print("=" * 60)
                
                print(f"Variant: {self.spliceai_results.get('variant')}")
                
                # Group by gene (since scores are usually the same per gene)
                gene_groups = {}
                for transcript in self.spliceai_results.get('splice_scores', []):
                    gene_id = transcript.get('ensembl_id', 'Unknown')
                    if gene_id not in gene_groups:
                        gene_groups[gene_id] = {
                            'scores': transcript,
                            'transcripts': []
                        }
                    gene_groups[gene_id]['transcripts'].append(transcript.get('transcript_id'))
                
                for gene_id, data in gene_groups.items():
                    print(f"\nGene: {gene_id}")
                    print(f"Gene name: {data['scores'].get('g_name', 'N/A')}")
                    print(f"Affected transcripts: {', '.join(data['transcripts'])}")
                    print(f"Transcript type: {data['scores'].get('t_type', 'N/A')}")
                    print(f"Transcript priority: {data['scores'].get('t_priority', 'N/A')}")
                    print(f"Strand: {data['scores'].get('t_strand', 'N/A')}")
                    if data['scores'].get('t_refseq_ids'):
                        print(f"RefSeq IDs: {', '.join(data['scores'].get('t_refseq_ids', []))}")
                    print()
                    
                    print("Splice Disruption Scores:")
                    print(f"  Acceptor Gain: {data['scores'].get('DS_AG')}")
                    print(f"  Acceptor Loss: {data['scores'].get('DS_AL')}")
                    print(f"  Donor Gain: {data['scores'].get('DS_DG')}")
                    print(f"  Donor Loss: {data['scores'].get('DS_DL')}")
                    print()
                    
                    print("Reference vs Alternative Scores:")
                    print(f"  Acceptor Gain - REF: {data['scores'].get('DS_AG_REF')} | ALT: {data['scores'].get('DS_AG_ALT')}")
                    print(f"  Acceptor Loss - REF: {data['scores'].get('DS_AL_REF')} | ALT: {data['scores'].get('DS_AL_ALT')}")
                    print(f"  Donor Gain - REF: {data['scores'].get('DS_DG_REF')} | ALT: {data['scores'].get('DS_DG_ALT')}")
                    print(f"  Donor Loss - REF: {data['scores'].get('DS_DL_REF')} | ALT: {data['scores'].get('DS_DL_ALT')}")
                    print()
                    
                    print("Position Changes:")
                    print(f"  Acceptor Gain: {data['scores'].get('DP_AG')}")
                    print(f"  Acceptor Loss: {data['scores'].get('DP_AL')}")
                    print(f"  Donor Gain: {data['scores'].get('DP_DG')}")
                    print(f"  Donor Loss: {data['scores'].get('DP_DL')}")
                    print()
                    
                    # Highlight significant scores (>0.5)
                    significant = []
                    for score_type, label in [
                        ('DS_AG', 'Acceptor Gain'),
                        ('DS_AL', 'Acceptor Loss'), 
                        ('DS_DG', 'Donor Gain'),
                        ('DS_DL', 'Donor Loss')
                    ]:
                        score = float(data['scores'].get(score_type, 0))
                        if score > 0.5:
                            significant.append(f"{label} ({score})")
                    
                    if significant:
                        print(f"Significant scores (>0.5): {', '.join(significant)}")
                    else:
                        print("No significant splice disruption predicted (all scores ≤0.5)")
                    
                    print("-" * 60)
                
        except Exception as e:
            print(f"Error: Failed to get SpliceAI results: {e}")
        
        self.wait_for_user()

    def get_pangolin_results(self) -> None:
        """Get Pangolin results for the current variant."""
        print("\n--- Getting Pangolin Results ---")
        
        if not self.variant_coord or not self.hg:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        print(f"Analyzing variant: {self.variant_coord} (GRCh{self.hg})")
        print(f"Distance parameter: {self.distance}")
        print(f"Mask parameter: {self.mask} (0 = raw score, 1 = masked score)")
        
        try:
            self.pangolin_results = call_pangolin(self.variant_coord, self.hg, self.distance, self.mask)
            
            if 'error' in self.pangolin_results:
                print(f"Error calling Pangolin: {self.pangolin_results['error']}")
            else:
                print("\n✓ Pangolin Results:")
                print("=" * 60)
                
                print(f"Variant: {self.pangolin_results.get('variant')}")
                
                # Group by gene (since scores are usually the same per gene)
                gene_groups = {}
                for transcript in self.pangolin_results.get('pangolin_scores', []):
                    gene_id = transcript.get('ensembl_id', 'Unknown')
                    if gene_id not in gene_groups:
                        gene_groups[gene_id] = {
                            'scores': transcript,
                            'transcripts': []
                        }
                    gene_groups[gene_id]['transcripts'].append(transcript.get('transcript_id'))
                
                for gene_id, data in gene_groups.items():
                    print(f"\nGene: {gene_id}")
                    print(f"Gene name: {data['scores'].get('g_name', 'N/A')}")
                    print(f"Affected transcripts: {', '.join(data['transcripts'])}")
                    print(f"Transcript type: {data['scores'].get('t_type', 'N/A')}")
                    print(f"Transcript priority: {data['scores'].get('t_priority', 'N/A')}")
                    print(f"Strand: {data['scores'].get('t_strand', 'N/A')}")
                    if data['scores'].get('t_refseq_ids'):
                        print(f"RefSeq IDs: {', '.join(data['scores'].get('t_refseq_ids', []))}")
                    print()
                    
                    print("Splice Disruption Scores:")
                    print(f"  Splice Gain: {data['scores'].get('DS_SG')}")
                    print(f"  Splice Loss: {data['scores'].get('DS_SL')}")
                    print()
                    
                    print("Reference vs Alternative Scores:")
                    print(f"  Splice Gain - REF: {data['scores'].get('SG_REF')} | ALT: {data['scores'].get('SG_ALT')}")
                    print(f"  Splice Loss - REF: {data['scores'].get('SL_REF')} | ALT: {data['scores'].get('SL_ALT')}")
                    print()
                    
                    print("Position Changes:")
                    print(f"  Splice Gain: {data['scores'].get('DP_SG')}")
                    print(f"  Splice Loss: {data['scores'].get('DP_SL')}")
                    print()
                    
                    # Highlight significant scores (>0.5)
                    significant = []
                    for score_type, label in [
                        ('DS_SG', 'Splice Gain'),
                        ('DS_SL', 'Splice Loss')
                    ]:
                        score_val = data['scores'].get(score_type)
                        if score_val is not None:
                            try:
                                score = float(score_val)
                                if abs(score) > 0.5:  # Use abs() since splice loss can be negative
                                    significant.append(f"{label} ({score})")
                            except (ValueError, TypeError):
                                pass
                    
                    if significant:
                        print(f"Significant scores (|score| > 0.5): {', '.join(significant)}")
                    else:
                        print("No significant splice disruption predicted (all |scores| ≤ 0.5)")
                    
                    # Show all non-zero scores if available
                    if self.pangolin_results.get('allNonZeroScores'):
                        print(f"\nAll affected positions within {self.distance}bp:")
                        for pos_data in self.pangolin_results['allNonZeroScores']:
                            pos = pos_data.get('pos')
                            sg_ref = pos_data.get('SG_REF')
                            sg_alt = pos_data.get('SG_ALT')
                            sl_ref = pos_data.get('SL_REF')
                            sl_alt = pos_data.get('SL_ALT')
                            print(f"  Position {pos}: SG({sg_ref}→{sg_alt}) SL({sl_ref}→{sl_alt})")
                    
                    print("-" * 60)
                
        except Exception as e:
            print(f"Error: Failed to get Pangolin results: {e}")

        self.wait_for_user()


    def select_tissue_type(self) -> None:
        """Select tissue type for GTEx analysis."""
        print("\n--- Tissue Type Selection ---")
        
        tissues_input = input("\nPlease enter the tissue types you wish to analyse, seperated by commas (e.g., Whole_Blood, Brain_Cortex)):\n").strip()
        
        if tissues_input:
            tissue_list = [tissue.strip() for tissue in tissues_input.split(',')]
            self.tissues = tissue_list
            print(f"✓ Tissue types set: {', '.join(tissue_list)}")
        else:
            print("Error: Please enter at least one valid tissue type")

        self.wait_for_user()

    def select_gtex_endpoint(self) -> None:
        """Select which GTEx API endpoint to use."""
        print("\n--- GTEx API Endpoint Selection ---")
        print("\nAvailable GTEx endpoints:")
        print("1. Median Gene Expression (default)")
        print("2. Gene Expression (with sample data)")
        print("3. Median Exon Expression")
        print("4. Median Junction Expression")
        print("5. Top Expressed Genes (by tissue)")
        
        choice = input("\nSelect endpoint (1-5, default=1): ").strip()
        
        endpoint_map = {
            "1": "medianGeneExpression",
            "2": "geneExpression",
            "3": "medianExonExpression",
            "4": "medianJunctionExpression",
            "5": "topExpressedGenes"
        }
        
        if choice in endpoint_map:
            self.gtex_endpoint = endpoint_map[choice]
            print(f"✓ GTEx endpoint set to: {self.gtex_endpoint}")
        else:
            print("Invalid choice, keeping default: medianGeneExpression")
        
        self.wait_for_user()


    def get_gtex_results(self) -> None:
        """Get GTEx expression results for the current variant and tissue."""
        print("\n--- Getting GTEx Results ---")
        
        if not self.tissues:
            print("Error: Please select tissue types first (Option 4)")
            return
        
        if not self.variant_coord:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        if self.gtex_endpoint != "topExpressedGenes" and not self.spliceai_results:
            print("Error: Please get SpliceAI results first (Option 2) to identify the gene")
            return
        
        print(f"Analyzing expression in {', '.join(self.tissues)} for variant: {self.variant_coord}")
        print(f"Using GTEx endpoint: {self.gtex_endpoint}")
        
        try:
            self.gtex_results = {}
            
            # Handle topExpressedGenes differently (tissue-based, not gene-based)
            if self.gtex_endpoint == "topExpressedGenes":
                for tissue in self.tissues:
                    print(f"\nQuerying top expressed genes in {tissue}")
                    result = topExpressedGenes(tissue)
                    
                    if 'error' not in result:
                        self.gtex_results[tissue] = result
                        print(f"✓ Found {len(result.get('top_genes', []))} top expressed genes")
                        # Show top 5
                        for i, gene in enumerate(result.get('top_genes', [])[:5]):
                            print(f"  {i+1}. {gene['gene_symbol']} - {gene['median_expression']} TPM")
                    else:
                        print(f"Error: {result['error']}")
            else:
                # Extract GENCODE IDs from SpliceAI results
                gencode_ids = set()
                for transcript in self.spliceai_results.get('splice_scores', []):
                    gene_id = transcript.get('ensembl_id')
                    if gene_id:
                        gencode_ids.add(gene_id)
                
                if not gencode_ids:
                    print("Error: No GENCODE IDs found in SpliceAI results")
                    return
                
                # Get GTEx results for each gene
                for gencode_id in gencode_ids:
                    print(f"\nQuerying GTEx for gene: {gencode_id}")
                    
                    # Call appropriate endpoint
                    if self.gtex_endpoint == "medianGeneExpression":
                        result = medianGeneExpression(gencode_id, self.tissues)
                    elif self.gtex_endpoint == "geneExpression":
                        result = geneExpression(gencode_id, self.tissues)
                    elif self.gtex_endpoint == "medianExonExpression":
                        result = medianExonExpression(gencode_id, self.tissues)
                    elif self.gtex_endpoint == "medianJunctionExpression":
                        result = medianJunctionExpression(gencode_id, self.tissues)
                    
                    if 'error' not in result:
                        self.gtex_results[gencode_id] = result
                        
                        # Display results based on endpoint
                        print(f"\n✓ GTEx Results for {gencode_id}:")
                        print("=" * 50)
                        
                        if self.gtex_endpoint == "medianGeneExpression":
                            for tissue, data in result.get('expression_data', {}).items():
                                print(f"\nTissue: {tissue}")
                                print(f"  Median TPM: {data.get('median_tpm', 'N/A')}")
                                print(f"  Expressed: {data.get('expressed', 'N/A')}")
                                print(f"  Expression Level: {data.get('expression_level', 'N/A')}")
                        
                        elif self.gtex_endpoint == "geneExpression":
                            for tissue, data in result.get('expression_data', {}).items():
                                print(f"\nTissue: {tissue}")
                                print(f"  Sample count: {len(data.get('data', []))}")
                                print(f"  Gene symbol: {data.get('gene_symbol', 'N/A')}")
                        
                        elif self.gtex_endpoint == "medianExonExpression":
                            for tissue, data in result.get('exon_data', {}).items():
                                print(f"\nTissue: {tissue}")
                                print(f"  Exon count: {data.get('exon_count', 0)}")
                                # Show first few exons
                                for exon in data.get('exons', [])[:3]:
                                    print(f"  - {exon['exon_id']}: {exon['median_expression']}")
                        
                        elif self.gtex_endpoint == "medianJunctionExpression":
                            for tissue, data in result.get('junction_data', {}).items():
                                print(f"\nTissue: {tissue}")
                                print(f"  Junction count: {data.get('junction_count', 0)}")
                    else:
                        print(f"Error getting GTEx data: {result['error']}")
            
        except Exception as e:
            print(f"Error: Failed to get GTEx results: {e}")
        
        self.wait_for_user()

    def input_context(self) -> None:
        "Provide context to the LLM to tweak query/output"
        print("\n--- Context for LLM ---")

        context = input("\nPlease enter context you would like to provide to the LLM\n")
        if context:
            self.context = context
            print(f"✓ Context")
        else:
            print("Error: No context detected")
        

    def get_llm_results(self) -> None:
        """Get LLM analysis in natural language."""
        print("\n--- Getting LLM Analysis ---")
        
        if not self.spliceai_results:
            print("Error: Please get SpliceAI results first (Option 2)")
            return
        
        if not self.gtex_results:
            print("Error: Please get GTEx results first (Option 4)")
            return
        
        try:
            llm_results = call_llm(
                self.spliceai_results, 
                self.gtex_results, 
                self.chrom, 
                self.pos, 
                self.ref, 
                self.alt, 
                self.context,
                model="gpt-4.1-nano"
            )
            
            self.print_results(llm_results)
            
        except Exception as e:
            print(f"Error: Failed to get LLM analysis: {e}")

        self.wait_for_user()

    def chat_with_llm(self) -> None:
        """Interactive chat with LLM about the variant."""
        self.chat_handler.start_chat()

    def print_results(self, llm_results: Dict[str, Any]) -> None:
        """
        Print formatted analysis results.
        
        Args:
            llm_results: Dictionary containing LLM analysis results
        """
        print(f"\nResult of a mutation of {self.ref} > {self.alt} at position {self.pos} on chromosome {self.chrom}")
        print("=" * 80)
        
        print(f"Priority: {llm_results.get('priority_level', 'N/A')}")
        print(f"Assessment: {llm_results.get('pathogenicity_assessment', 'N/A')}")
        print(f"Recommendations: {llm_results.get('experimental_recommendations', 'N/A')}")
        
        if 'error' in llm_results:
            print(f"\nError: {llm_results['error']}")
        else:
            print("\nFull Analysis:")
            print("-" * 40)
            print(llm_results.get('llm_interpretation', 'No interpretation available'))

    def display_menu(self) -> str:
        """Display the main menu and get user choice."""
        print("\n" + "=" * 60)
        print("Welcome to ChatSAV")
        print("This assistant analyses splice-altering variants (SAVs) using SpliceAI and GTEx data.")
        print("=" * 60)
        
        # Show current state
        status = []
        if self.variant_coord:
            status.append(f"✓ Variant: {self.variant_coord} (GRCh{self.hg})")
        if self.tissues:
            status.append(f"✓ Tissues: {', '.join(self.tissues)}")
        if self.spliceai_results:
            status.append("✓ SpliceAI results available")
        if self.pangolin_results:
            status.append("✓ Pangolin results available")
        if self.gtex_results:
            status.append("✓ GTEx results available")
        
        if status:
            print("\nCurrent Status:")
            for s in status:
                print(f"  {s}")
        
        print("\nPlease select an option:")
        print("1. Input variant coordinates")
        print("2. Get SpliceAI results")
        print("3. Get Pangolin results")
        print("4. Select tissue types")
        print("5. Select GTEx endpoint")
        print("6. Get GTEx results")
        print("7. Get LLM analysis in natural language")
        print("8. Input context for LLM")
        print("9. Chat with LLM")
        print("0. Exit")
        print("-" * 60)
        
        return input("Your choice: ").strip()

    def run(self) -> None:
        """Main execution loop."""
        while True:
            choice = self.display_menu()
            
            try:
                match choice:
                    case "1":
                        self.input_variant_coordinates()
                    case "2":
                        self.get_spliceai_results()
                    case "3":
                        self.get_pangolin_results()
                    case "4":
                        self.select_tissue_type()
                    case "5":
                        self.select_gtex_endpoint()
                    case "6":
                        self.get_gtex_results()
                    case "7":
                        self.get_llm_results()
                    case "8":
                        self.input_context()
                    case "9":
                        self.chat_with_llm()
                    case "0":
                        print("\nThank you for using ChatSAV!")
                        break
                    case _:
                        print("Invalid choice. Please select a number from 0-6.")
            
            except KeyboardInterrupt:
                print("\n\nExiting ChatSAV...")
                break
            except Exception as e:
                print(f"\nAn unexpected error occurred: {e}")
                print("Please try again or contact support if the problem persists.")


def main() -> None:
    """Main function to run the ChatSAV analysis."""
    pipeline = ChatSAVPipeline()
    pipeline.run()


if __name__ == "__main__":
    main()