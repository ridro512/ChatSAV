#!/usr/bin/env python3
"""
callGtex.py - File for calling GTEx whole blood expression data

Input:
Output: GTEx whole blood expression status
"""
import requests
from typing import List, Dict, Any, Optional

def call_gtex(gencode_id, tissues):
    """
    Query the GTEx API to check median transcript expression in a specific tissue.

    Args:
        gencode_id (list[str]): Versioned GENCODE IDs (e.g. 'ENSG00000123456.15').
        tissue (list[str]): GTEx tissue site IDs (e.g. 'Whole_Blood').
        dataset_id (str): Dataset ID (default is 'gtex_v8').

    Returns:
        float or None: Median TPM value if found, otherwise None.
        splice junctions and exons
    """
    results_list = []
    for tissue in tissues:
        url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        params = {
            "gencodeId": [gencode_id],
            "tissueSiteDetailId": [tissue],
            "datasetId": "gtex_v8"
        }

        try:
            response = requests.get(url, params=params, timeout=200)
            response.raise_for_status()
            result = response.json()
            data = result.get("data", [])
            if not data:
                # print(f"No expression data found for {gencode_id} in {tissue}.")
                results_list.append(None)
                continue

            median_tpm = data[0].get("median")
            # print(f"Median expression (TPM) of {gencode_id} in {tissue}: {median_tpm}")
            results_list.append(median_tpm)

        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            results_list.append(None)
            continue
        return results_list

# example call
# if __name__ == "__main__":
#     call_gtex(["ENSG00000141510.16"], "Whole_Blood")
