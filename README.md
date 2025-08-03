![image](https://github.com/user-attachments/assets/5ec8bf02-9ff5-41e0-bda5-4f571991a106)

# ChatSAV
LLM-Powered Assistant for Evaluating Candidate Disease-Causing Splice-Altering Variants

Genetic variants that disrupt pre-mRNA splicing (splice-altering variants or SAVs) are a major but underdiagnosed cause of genetic disease. The current diagnostic workflow is slow, manual, and highly expert-driven, requiring analysts to gather scattered data to interpret and assess experimental feasibility for each candidate variant.

Our mission is to make this workflow more streamlined using OpenAI's LLM. We have integrated several API's to assist variant curators in more efficiently and accurately predicting potential SAV's. These tools include:

- SpliceAI/Pangolin: Splicing prediction scores
- GTEx: Tissue specific expression and isoform data
- Google's AlphaGenome: Predicts mRNA and protein sequence changes

Credits/Copyright: Dao M, Venkateswaran G, Rogers C, Zhao K
## Installation

1. **Install dependencies**

   You'll need the following Python packages:

   - [`openai`](https://pypi.org/project/openai/):
     ```bash
     pip install openai
     ```

   - [`alphagenome`](https://pypi.org/project/alphagenome/):
     ```bash
     pip install alphagenome
     ``` 
     > *Note: AlphaGenome requires Python ≥ 3.10*
     
2. **Set up API keys**

   Before running the program, make sure to set your API keys as environment variables:

   - **Windows:**
     ```powershell
     $env:OPENAI_API_KEY="your-openai-key"
     $env:ALPHA_GENOME_API_KEY="your-alphagenome-key"
     ```

   - **macOS / Linux:**
     ```bash
     export OPENAI_API_KEY="your-openai-key"
     export ALPHA_GENOME_API_KEY="your-alphagenome-key"
     ```

3. **Make the script executable**

   ```bash
   chmod +x chatSav.py

# Usage
Input: Coordinates of a variant (chr:pos:ref:alt format)

VARIANT FORMAT INFORMATION:
- Chromosome: chr
- Position: pos
- Reference allele: ref
- Alternate allele: alt

