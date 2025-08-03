#!/usr/bin/env python3
"""
callGtex.py - File for calling GTEx expression data with correct API endpoints

Input: GENCODE ID and list of tissues
Output: GTEx expression data from various endpoints
"""
import requests
from typing import List, Dict, Any, Optional

def medianGeneExpression(gencode_id: str, tissues: List[str]) -> Dict[str, Any]:
    """
    Query the GTEx API for median gene expression across tissues.
    Uses: /api/v2/expression/medianGeneExpression

    Args:
        gencode_id (str): Versioned GENCODE ID (e.g. 'ENSG00000123456.15').
        tissues (List[str]): GTEx tissue site IDs (e.g. ['Whole_Blood', 'Brain_Cortex']).

    Returns:
        Dict containing median TPM values by tissue
    """
    results = {
        "gene_id": gencode_id,
        "tissues": tissues,
        "expression_data": {},
        "api_endpoint": "medianGeneExpression"
    }
    
    url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
    
    # GTEx API accepts arrays in the gencodeId parameter
    params = {
        "gencodeId": [gencode_id],
        "tissueSiteDetailId": tissues,
        "datasetId": "gtex_v10"
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        
        # Parse paginated response
        data = result.get("data", [])
        
        for item in data:
            tissue = item.get("tissueSiteDetailId")
            median_tpm = item.get("median")
            
            results["expression_data"][tissue] = {
                "median_tpm": median_tpm,
                "expressed": median_tpm > 1.0 if median_tpm else False,
                "expression_level": categorize_expression(median_tpm),
                "gencode_id": item.get("gencodeId"),
                "gene_symbol": item.get("geneSymbol"),
                "unit": item.get("unit")
            }

    except requests.exceptions.RequestException as e:
        print(f"Request failed for medianGeneExpression {gencode_id}: {e}")
        results["error"] = str(e)
    
    return results

def geneExpression(gencode_id: str, tissues: List[str] = None) -> Dict[str, Any]:
    """
    Query the GTEx API for gene expression across tissues.
    Uses: /api/v2/expression/geneExpression

    Args:
        gencode_id (str): Versioned GENCODE ID (e.g. 'ENSG00000123456.15').
        tissues (List[str], optional): Specific tissues to query.

    Returns:
        Dict containing expression data across tissues
    """
    results = {
        "gene_id": gencode_id,
        "expression_data": {},
        "api_endpoint": "geneExpression"
    }
    
    url = "https://gtexportal.org/api/v2/expression/geneExpression"
    params = {
        "gencodeId": [gencode_id],
        "datasetId": "gtex_v10"
    }
    
    # Add tissue filter if specified
    if tissues:
        params["tissueSiteDetailId"] = tissues

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        data = result.get("data", [])
        
        for tissue_data in data:
            tissue = tissue_data.get("tissueSiteDetailId")
            
            results["expression_data"][tissue] = {
                "gencode_id": tissue_data.get("gencodeId"),
                "gene_symbol": tissue_data.get("geneSymbol"),
                "unit": tissue_data.get("unit"),
                "data": tissue_data.get("data", []),  # Individual sample data
                "ontology_id": tissue_data.get("ontologyId"),
                "subset_group": tissue_data.get("subsetGroup")
            }

    except requests.exceptions.RequestException as e:
        print(f"Request failed for geneExpression {gencode_id}: {e}")
        results["error"] = str(e)
    
    return results

def medianExonExpression(gencode_id: str, tissues: List[str]) -> Dict[str, Any]:
    """
    Query the GTEx API for median exon expression data.
    Uses: /api/v2/expression/medianExonExpression

    Args:
        gencode_id (str): Versioned GENCODE ID
        tissues (List[str]): GTEx tissue site IDs

    Returns:
        Dict containing exon expression data
    """
    results = {
        "gene_id": gencode_id,
        "tissues": tissues,
        "exon_data": {},
        "api_endpoint": "medianExonExpression"
    }
    
    url = "https://gtexportal.org/api/v2/expression/medianExonExpression"
    params = {
        "gencodeId": [gencode_id],
        "tissueSiteDetailId": tissues,
        "datasetId": "gtex_v10"
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        data = result.get("data", [])
        
        # Group by tissue
        for item in data:
            tissue = item.get("tissueSiteDetailId")
            if tissue not in results["exon_data"]:
                results["exon_data"][tissue] = {"exons": [], "exon_count": 0}
            
            exon_info = {
                "exon_id": item.get("exonId"),
                "median_expression": item.get("median"),
                "gencode_id": item.get("gencodeId"),
                "gene_symbol": item.get("geneSymbol"),
                "unit": item.get("unit")
            }
            results["exon_data"][tissue]["exons"].append(exon_info)
        
        # Count exons per tissue
        for tissue in results["exon_data"]:
            results["exon_data"][tissue]["exon_count"] = len(results["exon_data"][tissue]["exons"])

    except requests.exceptions.RequestException as e:
        print(f"Request failed for medianExonExpression {gencode_id}: {e}")
        results["error"] = str(e)
    
    return results

def getExons(gencode_id: str) -> Dict[str, Any]:
    """
    Query the GTEx API for exon information for a gene.
    Uses: /api/v2/dataset/exon
    
    Args:
        gencode_id (str): Versioned GENCODE ID
        
    Returns:
        Dict containing exon structure information
    """
    results = {
        "gene_id": gencode_id,
        "exons": [],
        "api_endpoint": "exon"
    }
    
    url = "https://gtexportal.org/api/v2/dataset/exon"
    params = {
        "gencodeId": gencode_id,
        "datasetId": "gtex_v10"
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        data = result.get("data", [])
        
        for exon in data:
            exon_info = {
                "exon_id": exon.get("exonId"),
                "exon_number": exon.get("exonNumber"),
                "chromosome": exon.get("chromosome"),
                "start": exon.get("start"),
                "end": exon.get("end"),
                "strand": exon.get("strand"),
                "gencode_id": exon.get("gencodeId"),
                "transcript_id": exon.get("transcriptId"),
                "gene_symbol": exon.get("geneSymbol"),
                "length": exon.get("end", 0) - exon.get("start", 0) if exon.get("end") and exon.get("start") else None
            }
            results["exons"].append(exon_info)
        
        results["exon_count"] = len(results["exons"])
        
        # Sort exons by start position
        results["exons"].sort(key=lambda x: x.get("start", 0))
        
    except requests.exceptions.RequestException as e:
        print(f"Request failed for getExons {gencode_id}: {e}")
        results["error"] = str(e)
    
    return results

def medianJunctionExpression(gencode_id: str, tissues: List[str]) -> Dict[str, Any]:
    """
    Query the GTEx API for median junction expression data.
    Uses: /api/v2/expression/medianJunctionExpression

    Args:
        gencode_id (str): Versioned GENCODE ID
        tissues (List[str]): GTEx tissue site IDs

    Returns:
        Dict containing junction expression data
    """
    results = {
        "gene_id": gencode_id,
        "tissues": tissues,
        "junction_data": {},
        "api_endpoint": "medianJunctionExpression"
    }
    
    url = "https://gtexportal.org/api/v2/expression/medianJunctionExpression"
    params = {
        "gencodeId": [gencode_id],
        "tissueSiteDetailId": tissues,
        "datasetId": "gtex_v10"
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        data = result.get("data", [])
        
        # Group by tissue
        for item in data:
            tissue = item.get("tissueSiteDetailId")
            if tissue not in results["junction_data"]:
                results["junction_data"][tissue] = {"junctions": [], "junction_count": 0}
            
            junction_info = {
                "junction_id": item.get("junctionId"),
                "median_reads": item.get("median"),
                "gencode_id": item.get("gencodeId"),
                "gene_symbol": item.get("geneSymbol"),
                "unit": item.get("unit")
            }
            results["junction_data"][tissue]["junctions"].append(junction_info)
        
        # Count junctions per tissue
        for tissue in results["junction_data"]:
            results["junction_data"][tissue]["junction_count"] = len(results["junction_data"][tissue]["junctions"])

    except requests.exceptions.RequestException as e:
        print(f"Request failed for medianJunctionExpression {gencode_id}: {e}")
        results["error"] = str(e)
    
    return results

def topExpressedGenes(tissue: str, filter_mt_gene: bool = True) -> Dict[str, Any]:
    """
    Query the GTEx API for top expressed genes in a tissue.
    Uses: /api/v2/expression/topExpressedGene

    Args:
        tissue (str): GTEx tissue site ID
        filter_mt_gene (bool): Exclude mitochondrial genes

    Returns:
        Dict containing top expressed genes data
    """
    results = {
        "tissue": tissue,
        "top_genes": [],
        "api_endpoint": "topExpressedGene"
    }
    
    url = "https://gtexportal.org/api/v2/expression/topExpressedGene"
    params = {
        "tissueSiteDetailId": tissue,
        "filterMtGene": filter_mt_gene,
        "datasetId": "gtex_v10",
        "itemsPerPage": 100  # Get top 100 genes
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
        data = result.get("data", [])
        
        for item in data:
            gene_info = {
                "gencode_id": item.get("gencodeId"),
                "gene_symbol": item.get("geneSymbol"),
                "median_expression": item.get("median"),
                "unit": item.get("unit"),
                "tissue": item.get("tissueSiteDetailId"),
                "ontology_id": item.get("ontologyId")
            }
            results["top_genes"].append(gene_info)

    except requests.exceptions.RequestException as e:
        print(f"Request failed for topExpressedGenes {tissue}: {e}")
        results["error"] = str(e)
    
    return results

def categorize_expression(tpm_value: float) -> str:
    """Categorize expression level based on TPM value."""
    if tpm_value is None:
        return "no_data"
    elif tpm_value >= 10:
        return "high"
    elif tpm_value >= 1:
        return "moderate"
    elif tpm_value >= 0.1:
        return "low"
    else:
        return "very_low"

# Wrapper function for backward compatibility
def call_gtex(gencode_id: str, tissues: List[str]) -> Dict[str, Any]:
    """
    Main function for GTEx analysis - calls medianGeneExpression by default.
    
    Args:
        gencode_id (str): Versioned GENCODE ID
        tissues (List[str]): List of tissue types
    
    Returns:
        Dict containing expression results
    """
    return medianGeneExpression(gencode_id, tissues)

# Example usage
if __name__ == "__main__":
    # Test with the actual API endpoints
    test_gene = "ENSG00000167632.18"
    test_tissues = ["Whole_Blood", "Brain_Cortex"]
    
    print("Testing medianGeneExpression:")
    median_results = medianGeneExpression(test_gene, test_tissues)
    print(f"Results: {median_results}")
    
    print("\nTesting medianExonExpression:")
    exon_results = medianExonExpression(test_gene, test_tissues)
    print(f"Results: {exon_results}")
    
    print("\nTesting medianJunctionExpression:")
    junction_results = medianJunctionExpression(test_gene, test_tissues)
    print(f"Results: {junction_results}")
    
    print("\nTesting topExpressedGenes:")
    top_genes = topExpressedGenes("Whole_Blood")
    print(f"Top genes count: {len(top_genes.get('top_genes', []))}")