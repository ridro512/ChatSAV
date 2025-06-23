import requests

#Potential Errors: User inputs a build that doesnt have a corresponding ts
def call_splice(variant, hg):
    base_url = f"https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg={hg}&variant="
    url = base_url + variant

    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        results = []
        for score_set in data.get("scores", []):
            transcript_result = {
                "transcript_id": score_set.get("t_id"),
                "gene_id": score_set.get("gene_id"),
                "DS_AG": score_set.get("DS_AG"),
                "DS_AL": score_set.get("DS_AL"),
                "DS_DG": score_set.get("DS_DG"),
                "DS_DL": score_set.get("DS_DL"),
                "DP_AG": score_set.get("DP_AG"),
                "DP_AL": score_set.get("DP_AL"),
                "DP_DG": score_set.get("DP_DG"),
                "DP_DL": score_set.get("DP_DL")
            }
            results.append(transcript_result)
        return {
            "variant": variant,
            "splice_scores": results
        }
    elif response.status_code == 503:
        return f"Variant {variant} does not map to any transcript in the GRCh{hg} reference genome"
    else:
        return {
            "variant": variant,
            "error": f"API call failed with status code {response.status_code}"
        }
