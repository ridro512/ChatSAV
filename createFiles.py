import os
import re

def sanitise_name(name):
    """Sanitise folder/file names by replacing unsafe characters with underscores."""
    return re.sub(r'[^a-zA-Z0-9_.-]', '_', name)

def create_variant_folder(variant):
    """Create the main variant folder if it doesn't exist."""
    variant = sanitise_name(variant)
    if not os.path.exists(variant):
        os.makedirs(variant)
    return variant 

def add_splicing_results(variant, stdout):
    """Add combined SpliceAI and Pangolin output to splicing_results.txt inside variant folder."""
    variant_folder = create_variant_folder(variant)
    file_path = os.path.join(variant_folder, "splicing_results.txt")
    with open(file_path, "w") as f:
        f.write(stdout)

def pretty_dat_thang(expr_data):
    if isinstance(expr_data, dict):
        lines = []
        for tissue, data in expr_data.items():
            lines.append(f"{tissue}:")
            if isinstance(data, dict):
                for key, value in data.items():
                    lines.append(f"  {key}: {value}")
            else:
                lines.append(f"  {str(data)}")
        return "\n".join(lines)
    else:
        return str(expr_data)


def add_gtex_results(variant_coord, gtex_results):
    parent = sanitise_name(variant_coord)
    gtex_folder = os.path.join(parent, "GTEx")
    os.makedirs(parent, exist_ok=True)
    os.makedirs(gtex_folder, exist_ok=True)

    for tissue, endpoints in gtex_results.items():
        tissue_folder = os.path.join(gtex_folder, tissue)
        os.makedirs(tissue_folder, exist_ok=True)

        for endpoint, content in endpoints.items():
            file_path = os.path.join(tissue_folder, f"{endpoint}.txt")

            if os.path.exists(file_path):
                with open(file_path, "r") as f:
                    existing_content = f.read()
            else:
                existing_content = ""

            gene_id = content.get("gene_id", "Unknown Gene")
            tissues = content.get("tissues", [])
            expression_data = content.get("expression_data", {})
            api_endpoint = content.get("api_endpoint", "Unknown Endpoint")

            pretty_expr = pretty_dat_thang(expression_data)

            new_entry = (
                f"gene_id: {gene_id}\n"
                f"tissues: {tissues}\n"
                f"expression_data:\n{pretty_expr}\n"
                f"api_endpoint: {api_endpoint}\n\n"
            )

            if new_entry.strip() not in existing_content:
                with open(file_path, "a") as f:
                    if existing_content and not existing_content.endswith("\n"):
                        f.write("\n") 
                    f.write(new_entry)

def add_alphagenome_results(variant_coord, sequence_length, alpha_stdout):
    """Store AlphaGenome results in its own directory"""
    variant_folder = create_variant_folder(variant_coord)
    ag_folder = os.path.join(variant_folder, "AlphaGenome")
    os.makedirs(ag_folder, exist_ok=True)

    file_path = os.path.join(ag_folder, f"{sequence_length}.txt")
    with open(file_path, "w") as f:
        f.write(alpha_stdout)