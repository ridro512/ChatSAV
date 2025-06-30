import requests
from typing import Tuple, Dict, Any, Optional

#Potential Errors: User inputs a build that doesnt have a corresponding ts
def call_pangolin(variant, hg, distance, mask):
    url = f"https://pangolin-38-xwkwwwxdwq-uc.a.run.app/pangolin/?hg={hg}&variant={variant}&distance={distance}&mask={mask}"

    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        results = []
        for score_set in data.get("scores", []):
            transcript_result = {
                "transcript_id": score_set.get("t_id"),
                "ensembl_id": score_set.get("g_id"),
                "g_name": score_set.get("g_name"),  # Added missing gene name
                "t_type": score_set.get("t_type"),  # Added missing transcript type
                "t_priority": score_set.get("t_priority"),  # Added missing transcript priority
                "t_strand": score_set.get("t_strand"),  # Added missing strand
                "t_refseq_ids": score_set.get("t_refseq_ids"),  # Added missing RefSeq IDs
                "DS_SG": score_set.get("DS_SG"),
                "DS_SL": score_set.get("DS_SL"),
                "DP_SG": score_set.get("DP_SG"),
                "DP_SL": score_set.get("DP_SL"),
                "SG_REF": score_set.get("SG_REF"),
                "SG_ALT": score_set.get("SG_ALT"),
                "SL_REF": score_set.get("SL_REF"),
                "SL_ALT": score_set.get("SL_ALT")
            }
            results.append(transcript_result)
        return {
            "variant": variant,
            "pangolin_scores": results
        }
    elif response.status_code == 503:
        return f"Variant {variant} does not map to any transcript in the GRCh{hg} reference genome"
    else:
        return {
            "variant": variant,
            "error": f"API call failed with status code {response.status_code}"
        }
