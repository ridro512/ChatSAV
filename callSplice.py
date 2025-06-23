import requests

def callSplice(variant):
    base_url = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant="
    url = base_url + variant

    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        with open("variantSpliceResults.py", "w") as f:
            for score_set in data.get("scores", []):
                f.write(f"Transcript ID: {score_set.get('t_id')}\n")
                f.write(f"DS_AG: {score_set.get('DS_AG')}\n")
                f.write(f"DS_AL: {score_set.get('DS_AL')}\n")
                f.write(f"DS_DG: {score_set.get('DS_DG')}\n")
                f.write(f"DS_DL: {score_set.get('DS_DL')}\n")
                f.write(f"DP_AG: {score_set.get('DP_AG')}\n")
                f.write(f"DP_AL: {score_set.get('DP_AL')}\n")
                f.write(f"DP_DG: {score_set.get('DP_DG')}\n")
                f.write(f"DP_DL: {score_set.get('DP_DL')}\n")
                f.write("---\n")
    else:
        with open("variantSpliceResults.py", "w") as f:
            f.write(f"Error: {response.status_code}\n")
            f.write(response.text + "\n")