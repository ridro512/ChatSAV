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

import io
import os
import sys
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
from createFiles import (
    add_splicing_results,
    sanitise_name
)

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

    VALID_GTEX_ENDPOINTS = {
        "medianGeneExpression": "Median Gene Expression",
        "geneExpression": "Gene Expression (with sample data)",
        "medianExonExpression": "Median Exon Expression", 
        "medianJunctionExpression": "Median Junction Expression",
        "topExpressedGenes": "Top Expressed Genes (by tissue)",
        "getExons": "Get Exons"

    }
    ALPHAGENOME_SEQUENCE_LENGTHS = ["2KB", "16KB", "100KB", "500KB", "1MB"]
    
    def __init__(self):
        self.variant_coord: Optional[str] = None
        self.hg: Optional[int] = 38
        self.distance: int = 50
        self.mask: int = 0
        self.store: bool = False
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
        self.gtex_endpoints: List[str] = ["medianGeneExpression"]
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
        
        while True:
            storeFiles = input(
                "(Optional) Please enter if you would like to store your output in the current directory (1) or only view them in this session (0) (default: 0)\n"
            ).strip()

            try:
                store = int(storeFiles)
                if store in [0,1]:
                    self.store = bool(store)
                    if store == 1:
                        print(f"✓ Output set to be stored in current directory")
                    if store == 0:
                        print(f"✓ Output will not be stored")
                    break
                else:
                    print(f"{store} is not a valid option, output will not be stored 0\n")
                    break
            except ValueError:
                print(f"{storeFiles} is not a valid option, output will not be stored 0\n")
                break
            

        self.wait_for_user()

    def get_splice_results(self) -> None:
        """Get SpliceAI results for the current variant."""
        print("\n--- Getting SpliceAI Results ---")

        if self.store:
            buffer = io.StringIO()
            sys_stdout = sys.stdout
            sys.stdout = buffer
        
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

        finally:
            if self.store:
                sys.stdout = sys_stdout
                combined_output = buffer.getvalue()
                buffer.close()

                add_splicing_results(self.variant_coord, combined_output)
                filepath = os.path.join(sanitise_name(self.variant_coord), "splicing_results.txt")
                print(f"\n✓ Splicing results stored in folder: {self.variant_coord}/splicing_results.txt")
                with open(filepath, 'r') as f:
                    print(f.read())
    
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

    def select_gtex_endpoints(self) -> None:
        """Select which GTEx API endpoints to use (multiple selection)."""
        print("\n--- GTEx API Endpoints Selection ---")
        print("\nAvailable GTEx endpoints:")
        
        # Display available endpoints with numbers
        endpoint_list = list(self.VALID_GTEX_ENDPOINTS.keys())
        for i, (endpoint_key, endpoint_desc) in enumerate(self.VALID_GTEX_ENDPOINTS.items(), 1):
            print(f"{i}. {endpoint_desc}")
        
        endpoints_input = input(
            "\nPlease enter the endpoint numbers you wish to use, separated by commas.\n"
            "Press Enter without typing to default to Median Gene Expression:\n"
        ).strip()
        
        if not endpoints_input:
            # User pressed enter without input - set default
            self.gtex_endpoints = ["medianGeneExpression"]
            print("✓ Defaulting to Median Gene Expression")
        else:
            # Parse and validate endpoint selections
            try:
                endpoint_numbers = [int(num.strip()) for num in endpoints_input.split(',')]
                
                # Validate endpoint numbers
                invalid_numbers = []
                valid_endpoints = []
                
                for num in endpoint_numbers:
                    if 1 <= num <= len(endpoint_list):
                        endpoint_key = endpoint_list[num - 1]
                        if endpoint_key not in valid_endpoints:  # Avoid duplicates
                            valid_endpoints.append(endpoint_key)
                    else:
                        invalid_numbers.append(num)
                
                if invalid_numbers:
                    print(f"\nWarning: The following endpoint numbers are invalid: {invalid_numbers}")
                    print(f"Valid numbers are 1-{len(endpoint_list)}")
                
                if valid_endpoints:
                    self.gtex_endpoints = valid_endpoints
                    endpoint_names = [self.VALID_GTEX_ENDPOINTS[ep] for ep in valid_endpoints]
                    print(f"\n✓ GTEx endpoints set: {', '.join(endpoint_names)}")
                else:
                    print("\nError: No valid endpoint numbers provided. Keeping current selection.")
                    
            except ValueError:
                print("\nError: Please enter valid numbers separated by commas. Keeping current selection.")

        self.wait_to_continue()
    
    def get_matching_gencode_id_from_spliceai(self, spliceai_id: str, tissues: List[str]) -> str:
        """
        Loop downward from SpliceAI output to find an older GTEx compatible version.
        """
        base_id, version_str = spliceai_id.split('.')
        version = int(version_str)

        for v in range(version, 0, -1):
            test_id = f"{base_id}.{v}"
            result = medianGeneExpression(test_id, tissues)

            if result.get("expression_data"):
                return test_id

        return "Failed (Placeholder)"

    def get_gtex_results(self) -> None:
        """Get GTEx expression results for the current variant and tissue using selected endpoints."""
        print("\n--- Getting GTEx Results ---")
        
        # Check if we need tissues for the selected endpoints
        needs_tissues = any(ep != "getExons" for ep in self.gtex_endpoints)
        if needs_tissues and not self.tissues:
            print("Error: Please select tissue types first (Option 1 in GTEx submenu)")
            return
        
        if not self.variant_coord:
            print("Error: Please input variant coordinates first (Option 1)")
            return
        
        # Check if we need SpliceAI results for gene-based endpoints
        needs_genes = any(ep != "topExpressedGenes" for ep in self.gtex_endpoints)
        if needs_genes and not self.spliceai_results:
            print("Error: Please get SpliceAI results first (Option 2) to identify the gene")
            return
        
        print(f"Analyzing expression in {', '.join(self.tissues) if self.tissues else 'N/A'} for variant: {self.variant_coord}")
        print(f"Using GTEx endpoints: {', '.join([self.VALID_GTEX_ENDPOINTS[ep] for ep in self.gtex_endpoints])}")
        
        try:
            self.gtex_results = {}
            
            # Process each selected endpoint
            for endpoint in self.gtex_endpoints:
                print(f"\n{'='*20} Processing {self.VALID_GTEX_ENDPOINTS[endpoint]} {'='*20}")
                
                endpoint_results = {}
                
                # Handle topExpressedGenes differently (tissue-based, not gene-based)
                if endpoint == "topExpressedGenes":
                    for tissue in self.tissues:
                        print(f"\nQuerying top expressed genes in {tissue}")
                        result = topExpressedGenes(tissue)
                        
                        if 'error' not in result:
                            endpoint_results[tissue] = result
                            print(f"✓ Found {len(result.get('top_genes', []))} top expressed genes")
                            # Show top 5
                            for i, gene in enumerate(result.get('top_genes', [])[:5]):
                                print(f"  {i+1}. {gene['gene_symbol']} - {gene['median_expression']} TPM")
                        else:
                            print(f"Error: {result['error']}")
                
                else:
                    # Extract GENCODE IDs from SpliceAI results for gene-based endpoints
                    gencode_ids = set()
                    if endpoint != "getExons" or self.spliceai_results:
                        for transcript in self.spliceai_results.get('splice_scores', []):
                            gene_id = transcript.get('ensembl_id')
                            if gene_id:
                                gencode_ids.add(gene_id)
                    
                    if not gencode_ids and endpoint != "getExons":
                        print(f"Error: No GENCODE IDs found in SpliceAI results for {endpoint}")
                        continue
                    
                    # Get GTEx results for each gene
                    for gencode_id in gencode_ids:
                        resolved_id = self.get_matching_gencode_id_from_spliceai(gencode_id, self.tissues)
                        print(f"\nQuerying GTEx {endpoint} for resolved gene: {resolved_id}")
                        # print(f"\nQuerying GTEx {endpoint} for gene: {gencode_id}")
                        
                        # Call appropriate endpoint
                        if endpoint == "medianGeneExpression":
                            result = medianGeneExpression(resolved_id, self.tissues)
                        elif endpoint == "geneExpression":
                            result = geneExpression(resolved_id, self.tissues)
                        elif endpoint == "medianExonExpression":
                            result = medianExonExpression(resolved_id, self.tissues)
                        elif endpoint == "medianJunctionExpression":
                            result = medianJunctionExpression(resolved_id, self.tissues)
                        elif endpoint == "getExons":
                            result = getExons(resolved_id)
                        
                        if 'error' not in result:
                            endpoint_results[resolved_id] = result
                            
                            # Display results based on endpoint
                            print(f"\n✓ GTEx Results for {resolved_id}:")
                            print("=" * 50)
                            
                            if endpoint == "medianGeneExpression":
                                for tissue, data in result.get('expression_data', {}).items():
                                    print(f"\nTissue: {tissue}")
                                    print(f"  Median TPM: {data.get('median_tpm', 'N/A')}")
                                    print(f"  Expressed: {data.get('expressed', 'N/A')}")
                                    print(f"  Expression Level: {data.get('expression_level', 'N/A')}")
                            
                            elif endpoint == "geneExpression":
                                for tissue, data in result.get('expression_data', {}).items():
                                    print(f"\nTissue: {tissue}")
                                    print(f"  Sample count: {len(data.get('data', []))}")
                                    print(f"  Gene symbol: {data.get('gene_symbol', 'N/A')}")
                            
                            elif endpoint == "medianExonExpression":
                                for tissue, data in result.get('exon_data', {}).items():
                                    print(f"\nTissue: {tissue}")
                                    print(f"  Exon count: {data.get('exon_count', 0)}")
                                    # Show first few exons
                                    for exon in data.get('exons', [])[:3]:
                                        print(f"  - {exon['exon_id']}: {exon['median_expression']}")
                            
                            elif endpoint == "medianJunctionExpression":
                                for tissue, data in result.get('junction_data', {}).items():
                                    print(f"\nTissue: {tissue}")
                                    print(f"  Junction count: {data.get('junction_count', 0)}")
                            
                            elif endpoint == "getExons":
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
                
                # Store results for this endpoint
                if endpoint_results:
                    self.gtex_results[endpoint] = endpoint_results
            
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
            
            if self.gtex_endpoints:
                endpoint_names = [self.VALID_GTEX_ENDPOINTS[ep] for ep in self.gtex_endpoints]
                print(f"Current endpoints: {', '.join(endpoint_names)}")
            else:
                print("Current endpoints: None selected")
            
            print("\nOptions:")
            print("1. Select tissue types")
            print("2. Select GTEx endpoints")
            print("3. Get GTEx results")
            print("0. Return to main menu")
            print("-" * 60)
            
            choice = input("Your choice: ").strip()
            
            if choice == "1":
                self.select_tissue_type()
            elif choice == "2":
                self.select_gtex_endpoints()  # Updated method name
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
            print("Error: Please get GTEx results first (Option 3)")
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
                model="gpt-4o"
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


    def input_context(self) -> None:
        "Provide context to the LLM to tweak query/output"
        print("\n--- Context for LLM ---")

        context = input("\nPlease enter context you would like to provide to the LLM\n")
        if context:
            self.context = context
            print(f"✓ Context")
        else:
            print("Error: No context detected")
        


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
        print("4. AlphaGenome")
        print("5. LLM")
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
                        self.alphagenome_submenu()
                    case "5":
                        self.llm_submenu()
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