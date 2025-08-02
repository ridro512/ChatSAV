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
