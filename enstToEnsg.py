import requests

def enst_to_ensg(enst_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{enst_id}?expand=0"
    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        return f"Invalid Request"
    decoded = response.json()
    return decoded.get("Parent")

# Testing
# if __name__ == "__main__":
#     test_enst = "ENST00000367770"  
#     gene_id = enst_to_ensg(test_enst)
#     if gene_id:
#         print(f"Transcript {test_enst} corresponds to gene {gene_id}")
#     else:
#         print(f"Gene ID for transcript {test_enst} not found")
