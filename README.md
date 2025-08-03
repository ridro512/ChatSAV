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

## Usage

ChatSAV supports the following functionality:

### 1. Input Variant

- **Format**: `chr:pos:ref:alt` (e.g. `chr1:123456:A:G`)
- **Genome build**: Choose `hg38` or `hg19` (default: `hg38`)
- **Distance**: Window size for SpliceAI/Pangolin (default: `50`)
- **Score type**: Raw or masked (default: `Raw`)
- **Storage**: Option to save results locally

### 2. SpliceAI & Pangolin

- Generates prediction scores from both tools for the variant

### 3. GTEx Integration

- Supports selection of **multiple endpoints** and **multiple tissues** per query.
- Available endpoints include:
  - Median Gene Expression
  - Gene Expression with sample data
  - Median Exon Expression
  - Median Junction Expression
  - Top Expressed Genes (by tissue)
  - Exon retrieval by gene

- Supports **56 GTEx tissues**

### 4. AlphaGenome

- Predicts transcript and protein sequences across multiple sequence lengths:
  - **2KB** – Fastest; good for quick tests/debugging
  - **16KB** – Balanced; captures nearby regulatory elements
  - **100KB** – Default; comprehensive analysis
  - **500KB** – In-depth; long-range regulatory coverage
  - **1MB** – Most complete; for high-confidence research

> *AlphaGenome performs best with GRCh38/hg38. Accuracy may be reduced on GRCh37/hg19.*  
> *Start with 100KB for standard use; use 2KB for testing and 500KB+ for critical variants.*
