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
                print(f"\n Warning: The following tissue IDs are not valid GTEx tissue types:")
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
            elif choice == "3":
                self.chat_with_llm()
            elif choice == "0":
                break
            else:
                print("Invalid choice. Please select a number from 0-3.")

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
        """Get AlphaGenome results for the current variant."""
        print("\n--- Getting AlphaGenome Results ---")
        
        if not self.variant_coord or not self.hg:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        print(f"Analyzing variant: {self.variant_coord} (GRCh{self.hg})")
        print(f"Sequence length: {self.alphagenome_sequence_length}")
        
        if self.hg == "37":
            print("⚠️  Warning: AlphaGenome is optimized for GRCh38. Results for GRCh37 may be less accurate.")
        
        try:
            print("Calling AlphaGenome API (this may take a few moments)...")
            self.alphagenome_results = call_alphagenome(
                self.variant_coord, 
                self.hg, 
                self.alphagenome_sequence_length
            )
            
            if 'error' in self.alphagenome_results:
                print(f"Error calling AlphaGenome: {self.alphagenome_results['error']}")
            else:
                print("\n✓ AlphaGenome Results:")
                print("=" * 70)
                
                print(f"Variant: {self.alphagenome_results.get('variant')}")
                print(f"Genome build: {self.alphagenome_results.get('genome_build')}")
                print(f"Sequence length: {self.alphagenome_results.get('sequence_length')}")
                
                # Add general explanation
                print("\n" + "ℹ️  INTERPRETATION GUIDE".center(70))
                print("• Quantile scores range from -1.0 to +1.0")
                print("• Scores >0.5 or <-0.5 are considered significant")
                print("• Scores >0.8 or <-0.8 indicate high confidence predictions")
                print("• Positive scores = increase, Negative scores = decrease")
                print("=" * 70)
                
                scores = self.alphagenome_results.get('alphagenome_scores', {})
                
                # Display splice predictions
                splice_preds = scores.get('splice_predictions', {})
                if splice_preds and 'error' not in splice_preds:
                    print("\n🧬 SPLICE PREDICTIONS:")
                    print("   (How the variant affects RNA splicing)")
                    print("-" * 50)
                    
                    for splice_type, data in splice_preds.items():
                        if isinstance(data, dict) and 'total_predictions' in data:
                            print(f"\n{splice_type.replace('_', ' ').title()}:")
                            if splice_type == 'splice_sites':
                                print("   → Predicts creation/disruption of splice sites")
                            elif splice_type == 'splice_site_usage':
                                print("   → Predicts if splice sites will be actively used")
                            
                            print(f"  Total predictions: {data.get('total_predictions', 0)}")
                            sig_count = data.get('significant_predictions', 0)
                            total_count = data.get('total_predictions', 0)
                            print(f"  Significant effects: {sig_count} ({sig_count/total_count*100:.1f}% of tissues)" if total_count > 0 else f"  Significant effects: {sig_count}")
                            
                            if data.get('mean_quantile') is not None:
                                mean_score = data.get('mean_quantile', 0)
                                print(f"  Mean quantile score: {mean_score:.3f}", end="")
                                if abs(mean_score) > 0.8:
                                    print(" (🔴 VERY HIGH impact)")
                                elif abs(mean_score) > 0.5:
                                    print(" (🟡 HIGH impact)")
                                elif abs(mean_score) > 0.2:
                                    print(" (🟢 MODERATE impact)")
                                else:
                                    print(" (⚪ LOW impact)")
                            
                            # Show top tissues if available
                            top_tissues = data.get('top_tissues', [])
                            if top_tissues:
                                print("  Top affected tissues:")
                                for tissue in top_tissues[:3]:
                                    if isinstance(tissue, dict):
                                        effect_score = tissue.get('max_effect', 0)
                                        effect_icon = "🔴" if abs(effect_score) > 0.8 else "🟡" if abs(effect_score) > 0.5 else "🟢"
                                        print(f"    {effect_icon} {tissue.get('tissue', 'Unknown')}: {effect_score:.3f}")
                
                # Display expression predictions
                expr_preds = scores.get('expression_predictions', {})
                if expr_preds and 'error' not in expr_preds:
                    print("\n📊 EXPRESSION PREDICTIONS:")
                    print("   (How the variant affects gene expression levels)")
                    print("-" * 50)
                    
                    for expr_type, data in expr_preds.items():
                        if isinstance(data, dict) and 'total_predictions' in data:
                            print(f"\n{expr_type.replace('_', ' ').upper()}:")
                            if expr_type == 'rna_seq':
                                print("   → Measures overall gene expression changes")
                            elif expr_type == 'cage':
                                print("   → Measures transcription start site activity")
                            
                            total = data.get('total_predictions', 0)
                            up_count = data.get('upregulated_count', 0)
                            down_count = data.get('downregulated_count', 0)
                            
                            print(f"  Total predictions: {total}")
                            print(f"  Upregulated: {up_count} ({up_count/total*100:.1f}% of tissues)" if total > 0 else f"  Upregulated: {up_count}")
                            print(f"  Downregulated: {down_count} ({down_count/total*100:.1f}% of tissues)" if total > 0 else f"  Downregulated: {down_count}")
                            
                            if data.get('mean_effect') is not None:
                                mean_effect = data.get('mean_effect', 0)
                                print(f"  Mean effect: {mean_effect:.3f}", end="")
                                if mean_effect < -0.8:
                                    print(" (🔴 STRONG downregulation)")
                                elif mean_effect < -0.5:
                                    print(" (🟡 MODERATE downregulation)")
                                elif mean_effect < -0.2:
                                    print(" (🟢 MILD downregulation)")
                                elif mean_effect > 0.8:
                                    print(" (🔴 STRONG upregulation)")
                                elif mean_effect > 0.5:
                                    print(" (🟡 MODERATE upregulation)")
                                elif mean_effect > 0.2:
                                    print(" (🟢 MILD upregulation)")
                                else:
                                    print(" (⚪ MINIMAL effect)")
                
                # Display chromatin predictions
                chromatin_preds = scores.get('chromatin_predictions', {})
                if chromatin_preds and 'error' not in chromatin_preds:
                    print("\n🔬 CHROMATIN ACCESSIBILITY PREDICTIONS:")
                    print("   (How the variant affects DNA accessibility for regulation)")
                    print("-" * 50)
                    
                    for chrom_type, data in chromatin_preds.items():
                        if isinstance(data, dict) and 'total_predictions' in data:
                            print(f"\n{chrom_type.upper()}:")
                            if chrom_type == 'atac':
                                print("   → ATAC-seq: Measures chromatin accessibility")
                            elif chrom_type == 'dnase':
                                print("   → DNase-seq: Measures DNA accessibility")
                            
                            total = data.get('total_predictions', 0)
                            inc_count = data.get('accessibility_increase', 0)
                            dec_count = data.get('accessibility_decrease', 0)
                            
                            print(f"  Total predictions: {total}")
                            print(f"  Accessibility increase: {inc_count} ({inc_count/total*100:.1f}% of tissues)" if total > 0 else f"  Accessibility increase: {inc_count}")
                            print(f"  Accessibility decrease: {dec_count} ({dec_count/total*100:.1f}% of tissues)" if total > 0 else f"  Accessibility decrease: {dec_count}")
                            
                            # Add interpretation
                            if inc_count > dec_count:
                                print("   💡 Net effect: More open chromatin (easier transcription factor binding)")
                            elif dec_count > inc_count:
                                print("   💡 Net effect: More closed chromatin (harder transcription factor binding)")
                            else:
                                print("   💡 Net effect: Balanced accessibility changes")
                
                # Display summary statistics
                summary = scores.get('summary_stats', {})
                if summary and 'error' not in summary:
                    print("\n📈 SUMMARY STATISTICS:")
                    print("-" * 50)
                    total_preds = summary.get('total_predictions', 0)
                    sig_effects = summary.get('significant_effects', 0)
                    mean_score = summary.get('mean_quantile_score', 0)
                    max_effect = summary.get('max_absolute_effect', 0)
                    tissues = summary.get('unique_tissues', 0)
                    
                    print(f"Total predictions: {total_preds}")
                    print(f"Significant effects: {sig_effects} ({sig_effects/total_preds*100:.1f}% of all predictions)" if total_preds > 0 else f"Significant effects: {sig_effects}")
                    print(f"Mean quantile score: {mean_score:.3f}", end="")
                    if abs(mean_score) > 0.5:
                        print(" (Overall significant impact)")
                    else:
                        print(" (Overall moderate impact)")
                    
                    print(f"Max absolute effect: {max_effect:.3f}", end="")
                    if max_effect > 0.8:
                        print(" (Very strong tissue-specific effects)")
                    elif max_effect > 0.5:
                        print(" (Strong tissue-specific effects)")
                    else:
                        print(" (Moderate tissue-specific effects)")
                    
                    print(f"Unique tissues analyzed: {tissues}")
                
                # Display top predictions
                top_preds = self.alphagenome_results.get('top_predictions', [])
                if top_preds:
                    print("\n🎯 TOP PREDICTIONS:")
                    print("   (Strongest predicted effects across all tissues)")
                    print("-" * 50)
                    for i, pred in enumerate(top_preds[:5], 1):  # Show top 5
                        score = pred.get('quantile_score', 0)
                        significance = pred.get('significance', 'unknown')
                        effect_dir = pred.get('effect_direction', 'effect')
                        
                        # Add visual indicators
                        if significance == 'high':
                            icon = "🔴"
                        elif significance == 'moderate':
                            icon = "🟡"
                        else:
                            icon = "🟢"
                        
                        print(f"{i}. {icon} {pred.get('output_type', 'Unknown')} in {pred.get('biosample_name', 'Unknown')}")
                        print(f"   Score: {score:.3f} ({significance} confidence {effect_dir})")
                
                # Show successful scorers info
                if self.alphagenome_results.get('successful_scorers'):
                    successful = self.alphagenome_results['successful_scorers']
                    total_attempted = self.alphagenome_results.get('total_scorers_attempted', len(successful))
                    print(f"\n🔧 MODEL PERFORMANCE:")
                    print(f"Successful predictors: {len(successful)}/{total_attempted}")
                    print(f"Active modules: {', '.join(successful)}")
                    
                    if len(successful) == total_attempted:
                        print("✅ All AlphaGenome prediction modules worked successfully")
                    elif len(successful) >= total_attempted * 0.8:
                        print("✅ Most AlphaGenome prediction modules worked successfully")
                    else:
                        print("⚠️ Some AlphaGenome prediction modules failed")
                
                # Add overall interpretation
                print(f"\n💡 OVERALL INTERPRETATION:")
                print("-" * 50)
                
                # Determine overall impact level
                summary = scores.get('summary_stats', {})
                if summary and 'error' not in summary:
                    sig_ratio = summary.get('significant_effects', 0) / max(summary.get('total_predictions', 1), 1)
                    max_effect = summary.get('max_absolute_effect', 0)
                    
                    if sig_ratio > 0.8 and max_effect > 0.8:
                        impact_level = "🔴 VERY HIGH IMPACT"
                        interpretation = "This variant is predicted to have strong, widespread effects across multiple tissues and mechanisms."
                    elif sig_ratio > 0.5 or max_effect > 0.5:
                        impact_level = "🟡 HIGH IMPACT"
                        interpretation = "This variant is predicted to have significant effects in multiple tissues."
                    elif sig_ratio > 0.2 or max_effect > 0.3:
                        impact_level = "🟢 MODERATE IMPACT"
                        interpretation = "This variant is predicted to have moderate effects in some tissues."
                    else:
                        impact_level = "⚪ LOW IMPACT"
                        interpretation = "This variant is predicted to have minimal functional effects."
                    
                    print(f"Impact Level: {impact_level}")
                    print(f"Summary: {interpretation}")
                    print("\n⚠️  Note: These are computational predictions. Experimental validation is recommended for clinical interpretation.")
                
        except Exception as e:
            print(f"Error: Failed to get AlphaGenome results: {e}")
        
        self.wait_for_user()

    def alphagenome_submenu(self) -> None:
        """submenu for alphagenome options"""
        while True:
            print("\n" + "=" * 60)
            print("AlphaGenome Submenu")
            print("=" * 60)

            print("\nOptions:")
            print("1. AlphaGenome Configuration")
            print("2. Get AlphaGenome Analysis")
            print("0. Return to main menu")
            print("-" * 60)
            
            choice = input("Your choice: ").strip()
            
            if choice == "1":
                self.configure_alphagenome()
            elif choice == "2":
                self.get_alphagenome_results()
            elif choice == "0":
                break
            else:
                print("Invalid choice. Please select a number from 0-2.")


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
        print("5. AlphaGenome")
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
                        self.alphagenome_submenu()
                    case "0":
                        print("\nThank you for using ChatSAV!")
                        break
                    case _:
                        print("Invalid choice. Please select a number from 0-5.")
            
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