#!/usr/bin/env python3
"""
callGtex.py - File for calling GTEx whole blood expression data

Input:
Output: GTEx whole blood expression status
"""
import requests

def call_gtex(gencode_id, tissue):
    """
    Query the GTEx API to check median transcript expression in a specific tissue.

    Args:
        gencode_id (str): Versioned GENCODE ID (e.g. 'ENSG00000123456.15').
        tissue (str): GTEx tissue site ID (e.g. 'Whole_Blood').
        dataset_id (str): Dataset ID (default is 'gtex_v8').

    Returns:
        float or None: Median TPM value if found, otherwise None.
    """
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
            return None

        median_tpm = data[0].get("median")
        # print(f"Median expression (TPM) of {gencode_id} in {tissue}: {median_tpm}")
        return median_tpm

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return None


# if __name__ == "__main__":
#     call_gtex("ENSG00000141510.16", "Whole_Blood")
