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

from typing import Tuple, Dict, Any, Optional, List

from callSplice import call_splice
from callPangolin import call_pangolin
from callAlphaGenome import call_alphagenome
from callGtex import (
    medianGeneExpression, 
    geneExpression, 
    medianExonExpression, 
    medianJunctionExpression,
    topExpressedGenes,
    getExons
)
from callLlm import call_llm
from chatLLM import ChatLLM


class ChatSAVPipeline:
    """Main pipeline class for ChatSAV analysis."""

    VALID_TISSUES = {
        "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", 
        "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Bladder", 
        "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", 
        "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
        "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
        "Brain_Hippocampus", "Brain_Hypothalamus", 
        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", 
        "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", 
        "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", 
        "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", 
        "Cervix_Ectocervix", "Cervix_Endocervix", "Colon_Sigmoid", 
        "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", 
        "Esophagus_Mucosa", "Esophagus_Muscularis", "Fallopian_Tube", 
        "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", 
        "Kidney_Medulla", "Liver", "Lung", "Minor_Salivary_Gland", 
        "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", 
        "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", 
        "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", 
        "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", 
        "Whole_Blood"
    }

    ALPHAGENOME_SEQUENCE_LENGTHS = ["2KB", "16KB", "100KB", "500KB", "1MB"]
    
    def __init__(self):
        self.variant_coord: Optional[str] = None
        self.hg: Optional[str] = 38
        self.distance: int = 50
        self.mask: int = 0
        self.tissues: Optional[List[str]] = []
        self.spliceai_results: Optional[Dict[str, Any]] = None
        self.pangolin_results: Optional[Dict[str, Any]] = None
        self.gtex_results: Optional[Dict[str, Any]] = None
        self.alphagenome_results: Optional[Dict[str, Any]] = None
        self.alphagenome_sequence_length: str = "100KB"
        self.context: Optional[str] = None
        self.chrom: Optional[str] = None
        self.pos: Optional[int] = None
        self.ref: Optional[str] = None
        self.alt: Optional[str] = None
        self.gtex_endpoint: str = "medianGeneExpression"
        self.last_llm_results: Optional[Dict[str, Any]] = None
        self.chat_handler = ChatLLM(self)

    def wait_for_user(self) -> None:
        """Wait for user to press a key before continuing."""
        input("\nPress Enter to return to main menu...")

    def wait_to_continue(self) -> None:
        "Wait for user to press a key and continues to next step"
        input("\nPress Enter to continue...")

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
                "Please enter the genome build [37 for GRCh37/hg19] or [38 for GRCh38/hg38] (default: hg38):\n"
            ).strip()
            
            if hg in ["37", "38"]:
                self.hg = hg
                print(f"✓ Genome build set: GRCh{hg}")
                break
            else:
                print(f"{hg} is not a valid build, defaulting to hg38\n")
                break
            
        while True:
            dist = input(
                "(Optional) Please enter the distance parameter of SpliceAI and Pangolin model (default: 50)\n"
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
                    print(f"✓ SpliceAI and Pangolin mask set: {mask}")
                    break
                else:
                    print(f"{mask} is not a valid option, defaulting to raw scores 0\n")
                    break
            except ValueError:
                print(f"{mask} is not a valid option, defaulting to raw scores 0\n")
                break

        self.wait_for_user()

    def get_splice_results(self) -> None:
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
                    print("-" * 60)
                
        except Exception as e:
            print(f"Error: Failed to get SpliceAI results: {e}")

        """Get Pangolin results for the current variant."""
        print("\n--- Getting Pangolin Results ---")

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

    def configure_alphagenome(self) -> None:
        """Configure AlphaGenome analysis parameters."""
        print("\n--- Configure AlphaGenome Parameters ---")
        
        print(f"Current sequence length: {self.alphagenome_sequence_length}")
        print("\nAvailable sequence lengths:")
        for i, length in enumerate(self.ALPHAGENOME_SEQUENCE_LENGTHS, 1):
            print(f"{i}. {length}")
        
        choice = input(f"\nSelect sequence length (1-{len(self.ALPHAGENOME_SEQUENCE_LENGTHS)}, default=current): ").strip()
        
        try:
            if choice:
                idx = int(choice) - 1
                if 0 <= idx < len(self.ALPHAGENOME_SEQUENCE_LENGTHS):
                    self.alphagenome_sequence_length = self.ALPHAGENOME_SEQUENCE_LENGTHS[idx]
                    print(f"✓ AlphaGenome sequence length set to: {self.alphagenome_sequence_length}")
                else:
                    print("Invalid choice, keeping current setting")
        except ValueError:
            print("Invalid choice, keeping current setting")
        
        self.wait_for_user()

    def get_alphagenome_results(self) -> None:
        """Get AlphaGenome variant effect predictions."""
        print("\n--- Getting AlphaGenome Results ---")
        
        if not self.variant_coord or not self.hg:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        print(f"Analyzing variant: {self.variant_coord} (GRCh{self.hg})")
        print(f"Sequence length: {self.alphagenome_sequence_length}")
        
        try:
            # Import the AlphaGenome function
            from callAlphaGenome import call_alphagenome
            
            # Use tissues if available, otherwise use default RNA_SEQ analysis
            tissues = None
            if self.tissues:
                # Convert GTEx tissue names to ontology terms if needed
                tissue_map = {
                    "Whole_Blood": "UBERON:0000178",
                    "Brain_Cortex": "UBERON:0000955", 
                    "Heart_Left_Ventricle": "UBERON:0000948",
                    "Liver": "UBERON:0002107",
                    "Lung": "UBERON:0002048",
                    "Muscle_Skeletal": "UBERON:0000993",
                    "Colon_Sigmoid": "UBERON:0001157",
                    "Colon_Transverse": "UBERON:0001157"
                }
                tissues = [tissue_map.get(tissue, tissue) for tissue in self.tissues]
            
            self.alphagenome_results = call_alphagenome(
                variant_coord=self.variant_coord,
                hg=self.hg,
                sequence_length=self.alphagenome_sequence_length,
                output_types=["RNA_SEQ", "ATAC", "DNASE", "CAGE", "PROCAP", "SPLICE_SITES"],  # Added splice-specific output
                tissues=tissues
            )
            
            if 'error' in self.alphagenome_results:
                print(f"Error calling AlphaGenome: {self.alphagenome_results['error']}")
            else:
                print("\n✓ AlphaGenome Results:")
                print("=" * 60)
                
                # Display basic variant information
                print(f"Variant: {self.alphagenome_results['variant']}")
                print(f"Genome Build: {self.alphagenome_results['genome_build']}")
                print(f"Sequence Length: {self.alphagenome_results['sequence_length']}")
                
                # Display interval information
                interval = self.alphagenome_results['interval']
                print(f"Analysis Interval: {interval['chromosome']}:{interval['start']}-{interval['end']}")
                print(f"Interval Length: {interval['length']:,} bp")
                
                # Display predictions for each output type
                predictions = self.alphagenome_results.get('predictions', {})
                for output_type, pred_data in predictions.items():
                    print(f"\n{output_type} Predictions:")
                    print("-" * 40)
                    
                    # Reference and alternate availability
                    if pred_data.get('reference_available'):
                        ref_summary = pred_data.get('reference_summary', {})
                        print(f"✓ Reference predictions available")
                        if 'track_count' in ref_summary:
                            print(f"  Reference tracks: {ref_summary['track_count']}")
                        if 'value_stats' in ref_summary:
                            stats = ref_summary['value_stats']
                            if 'error' not in stats:
                                print(f"  Reference stats: mean={stats.get('mean', 0):.4f}, max={stats.get('max', 0):.4f}")
                            else:
                                print(f"  Reference stats: {stats['error']}")
                    
                    if pred_data.get('alternate_available'):
                        alt_summary = pred_data.get('alternate_summary', {})
                        print(f"✓ Alternate predictions available")
                        if 'track_count' in alt_summary:
                            print(f"  Alternate tracks: {alt_summary['track_count']}")
                        if 'value_stats' in alt_summary:
                            stats = alt_summary['value_stats']
                            if 'error' not in stats:
                                print(f"  Alternate stats: mean={stats.get('mean', 0):.4f}, max={stats.get('max', 0):.4f}")
                            else:
                                print(f"  Alternate stats: {stats['error']}")
                    
                    # Prediction differences (variant effect)
                    pred_diff = pred_data.get('prediction_diff', {})
                    if pred_diff:
                        print(f"  Variant Effect:")
                        print(f"    Mean difference (ALT-REF): {pred_diff.get('mean_diff', 0):.4f}")
                        print(f"    Max difference (ALT-REF): {pred_diff.get('max_diff', 0):.4f}")
                        print(f"    Effect magnitude: {pred_diff.get('effect_magnitude', 0):.4f}")
                    
                    if pred_data.get('error'):
                        print(f"  Error: {pred_data['error']}")
                    
                    print("-" * 40)
                
                # Display variant scores for each output type
                variant_scores = self.alphagenome_results.get('variant_scores', {})
                for output_type, scores_list in variant_scores.items():
                    if scores_list:  # Check if there are any scores
                        print(f"\n{output_type} Variant Scores:")
                        print("-" * 40)
                        
                        for i, score_data in enumerate(scores_list):
                            if 'error' in score_data:
                                print(f"  Scorer {i+1}: Error - {score_data['error']}")
                                continue
                            
                            scorer_type = score_data.get('scorer_type', 'Unknown')
                            print(f"  Scorer {i+1} ({scorer_type}):")
                            print(f"    Total genes analyzed: {score_data.get('total_genes', 0)}")
                            print(f"    Total tracks analyzed: {score_data.get('total_tracks', 0)}")
                            
                            # Summary statistics
                            summary_stats = score_data.get('summary_stats', {})
                            if summary_stats and 'error' not in summary_stats:
                                print(f"    Mean score: {summary_stats.get('mean_score', 0):.4f}")
                                print(f"    Max absolute score: {summary_stats.get('max_score', 0):.4f}")
                                print(f"    Significant effects: {summary_stats.get('significant_scores', 0)}/{summary_stats.get('total_scores', 0)}")
                                
                                # Show valid score ratio if available
                                if 'valid_score_ratio' in summary_stats:
                                    print(f"    Valid score ratio: {summary_stats.get('valid_score_ratio', 1.0):.2f}")
                                
                                if summary_stats.get('note'):
                                    print(f"    Note: {summary_stats['note']}")
                            elif summary_stats and 'error' in summary_stats:
                                print(f"    Summary error: {summary_stats['error']}")
                            
                            # Splice-specific results for RNA_SEQ and SPLICE_SITES
                            splice_results = score_data.get('splice_specific_results', {})
                            if splice_results and 'error' not in splice_results:
                                print(f"    Splice Analysis:")
                                print(f"      Total splice effects: {splice_results.get('total_splice_effects', 0)}")
                                print(f"      Strong splice effects: {splice_results.get('strong_splice_effects', 0)}")
                                
                                # Show splice disruption genes
                                disruption_genes = splice_results.get('splice_disruption_genes', [])
                                if disruption_genes:
                                    print(f"      Splice disruption detected in {len(disruption_genes)} genes:")
                                    for gene in disruption_genes[:3]:  # Top 3
                                        print(f"        - {gene['gene_name']}: {gene['disruption_score']:.4f} ({gene['strength']})")
                                
                                # Show splice enhancement genes
                                enhancement_genes = splice_results.get('splice_enhancement_genes', [])
                                if enhancement_genes:
                                    print(f"      Splice enhancement detected in {len(enhancement_genes)} genes:")
                                    for gene in enhancement_genes[:3]:  # Top 3
                                        print(f"        - {gene['gene_name']}: {gene['enhancement_score']:.4f} ({gene['strength']})")
                            elif splice_results and 'error' in splice_results:
                                print(f"    Splice analysis error: {splice_results['error']}")
                            
                            # Top affected genes
                            top_genes = score_data.get('top_affected_genes', [])
                            if top_genes:
                                print(f"    Top 3 Affected Genes:")
                                for j, gene in enumerate(top_genes[:3], 1):
                                    gene_name = gene.get('gene_name', 'N/A')
                                    gene_id = gene.get('gene_id', 'N/A')
                                    max_score = gene.get('max_abs_score', 0.0)
                                    mean_score = gene.get('mean_abs_score', 0.0)
                                    
                                    print(f"      {j}. {gene_name} ({gene_id})")
                                    print(f"         Max effect: {max_score:.4f}, Mean effect: {mean_score:.4f}")
                            
                            print()
                    
                    print("-" * 40)
        
        except ImportError:
            print("Error: callAlphaGenome module not found. Please ensure callAlphaGenome.py is in the same directory.")
            print("Also ensure AlphaGenome is installed: pip install alphagenome")
        except Exception as e:
            print(f"Error: Failed to get AlphaGenome results: {e}")
        
        self.wait_for_user()

    def select_tissue_type(self) -> None:
        """Select tissue type for GTEx analysis."""
        print("\n--- Tissue Type Selection ---")
        
        print("\nAvailable tissue types (examples):")
        print("Whole_Blood, Brain_Cortex, Heart_Left_Ventricle, Liver, Lung, Muscle_Skeletal...")
        print("(See GTEx documentation for full list of valid tissue IDs)")
        
        tissues_input = input(
            "\nPlease enter the tissue types you wish to analyse, separated by commas.\n"
            "Press Enter without typing to default to Whole_Blood:\n"
        ).strip()
        
        if not tissues_input:
            # User pressed enter without input - set default
            self.tissues = ["Whole_Blood"]
            print("✓ Defaulting to Whole_Blood")
        else:
            # Parse and validate tissue types
            tissue_list = [tissue.strip() for tissue in tissues_input.split(',')]
            
            # Validate tissues
            invalid_tissues = []
            valid_tissues = []
            
            for tissue in tissue_list:
                if tissue in self.VALID_TISSUES:
                    valid_tissues.append(tissue)
                else:
                    invalid_tissues.append(tissue)
            
            if invalid_tissues:
                print(f"\n⚠️  Warning: The following tissue IDs are not valid GTEx tissue types:")
                for tissue in invalid_tissues:
                    print(f"   - {tissue}")
                print("\nValid tissue IDs include:")
                # Show a sample of valid tissues
                sample_tissues = sorted(list(self.VALID_TISSUES))[:10]
                for tissue in sample_tissues:
                    print(f"   - {tissue}")
                print("   ... and more (see GTEx documentation)")
            
            if valid_tissues:
                self.tissues = valid_tissues
                print(f"\n✓ Tissue types set: {', '.join(valid_tissues)}")
            else:
                print("\nError: No valid tissue types provided. Please try again.")

        self.wait_to_continue()

    def select_gtex_endpoint(self) -> None:
        """Select which GTEx API endpoint to use."""
        print("\n--- GTEx API Endpoint Selection ---")
        print("\nAvailable GTEx endpoints:")
        print("1. Median Gene Expression (default)")
        print("2. Gene Expression (with sample data)")
        print("3. Median Exon Expression")
        print("4. Median Junction Expression")
        print("5. Top Expressed Genes (by tissue)")
        print("6. Get Exons")
        
        choice = input("\nSelect endpoint (1-6, default=1): ").strip()
        
        endpoint_map = {
            "1": "medianGeneExpression",
            "2": "geneExpression",
            "3": "medianExonExpression",
            "4": "medianJunctionExpression",
            "5": "topExpressedGenes",
            "6": "getExons"
        }
        
        if choice in endpoint_map:
            self.gtex_endpoint = endpoint_map[choice]
            print(f"✓ GTEx endpoint set to: {self.gtex_endpoint}")
        else:
            print("Invalid choice, keeping default: medianGeneExpression")
        
        self.wait_to_continue()


    def get_gtex_results(self) -> None:
        """Get GTEx expression results for the current variant and tissue."""
        print("\n--- Getting GTEx Results ---")
        
        if not self.tissues and self.gtex_endpoint != "getExons":
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
                    elif self.gtex_endpoint == "getExons":
                        result = getExons(gencode_id)
                    
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
                        elif self.gtex_endpoint == "getExons":
                            print(f"  Total exons: {result.get('exon_count', 0)}")
                            for i, exon in enumerate(result.get('exons', [])[:5]):  # Show first 5
                                exon_num = exon.get('exon_number', i+1)
                                start = exon.get('start', 'N/A')
                                end = exon.get('end', 'N/A')
                                length = exon.get('length', 'N/A')
                                print(f"  - Exon {exon_num}: {start}-{end} ({length} bp)")
                            if result.get('exon_count', 0) > 5:
                                print(f"  ... and {result.get('exon_count') - 5} more exons")
                    else:
                        print(f"Error getting GTEx data: {result['error']}")
            
        except Exception as e:
            print(f"Error: Failed to get GTEx results: {e}")
        
        self.wait_to_continue()

    def gtex_submenu(self) -> None:
        """GTEx analysis submenu."""
        while True:
            print("\n" + "=" * 60)
            print("GTEx Analysis Submenu")
            print("=" * 60)
            
            # Show current GTEx status
            if self.tissues:
                print(f"Current tissues: {', '.join(self.tissues)}")
            else:
                print("Current tissues: None selected")
            print(f"Current endpoint: {self.gtex_endpoint}")
            
            print("\nOptions:")
            print("1. Select tissue types")
            print("2. Select GTEx endpoint")
            print("3. Get GTEx results")
            print("0. Return to main menu")
            print("-" * 60)
            
            choice = input("Your choice: ").strip()
            
            if choice == "1":
                self.select_tissue_type()
            elif choice == "2":
                self.select_gtex_endpoint()
            elif choice == "3":
                self.get_gtex_results()
            elif choice == "0":
                break
            else:
                print("Invalid choice. Please select a number from 0-3.")


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
            print("Error: Please get GTEx results first (Option 6)")
            return
        
        try:
            llm_results = call_llm(
                self.spliceai_results, 
                self.gtex_results or {}, 
                self.pangolin_results,
                self.alphagenome_results,
                self.chrom, 
                self.pos, 
                self.ref, 
                self.alt, 
                self.context,
                model="gpt-4.1-nano"
            )
            self.last_llm_results = llm_results
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

    def llm_submenu(self) -> None:
        """LLM analysis submenu."""
        while True:
            print("\n" + "=" * 60)
            print("LLM Analysis Submenu")
            print("=" * 60)
            
            # Show current LLM status
            if self.context:
                print(f"Current context: {self.context[:50]}..." if len(self.context) > 50 else f"Current context: {self.context}")
            else:
                print("Current context: None")
            
            if self.last_llm_results:
                print("✓ Previous LLM analysis available")
            
            print("\nOptions:")
            print("1. Input context for LLM")
            print("2. Get LLM analysis")
            print("3. Chat with LLM")
            print("0. Return to main menu")
            print("-" * 60)
            
            choice = input("Your choice: ").strip()
            
            if choice == "1":
                self.input_context()
            elif choice == "2":
                self.get_llm_results()
                self.wait_for_user()
            elif choice == "3":
                self.chat_with_llm()
            elif choice == "0":
                break
            else:
                print("Invalid choice. Please select a number from 0-3.")

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
        if self.alphagenome_results:
            status.append("✓ AlphaGenome results available")
        if self.gtex_results:
            status.append("✓ GTEx results available")
        
        if status:
            print("\nCurrent Status:")
            for s in status:
                print(f"  {s}")
        
        print("\nPlease select an option:")
        print("1. Input variant coordinates")
        print("2. SpliceAI and Pangolin ")
        print("3. GTEx")
        print("4. LLM")
        print("5. ALPHAGENOME PARAM")
        print("6. ALPHAGENOME RESULTS")
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
                        self.get_splice_results()
                    case "3":
                        self.gtex_submenu()
                    case "4":
                        self.llm_submenu()
                    case "5":
                        self.configure_alphagenome()
                    case "6":
                        self.get_alphagenome_results()
                    case "0":
                        print("\nThank you for using ChatSAV!")
                        break
                    case _:
                        print("Invalid choice. Please select a number from 0-11.")
            
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