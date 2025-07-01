![image](https://github.com/user-attachments/assets/5ec8bf02-9ff5-41e0-bda5-4f571991a106)

# ChatSAV
LLM-Powered Assistant for Evaluating Candidate Disease-Causing Splice-Altering Variants

Genetic variants that disrupt pre-mRNA splicing (splice-altering variants or SAVs) are a frequent cause of genetic diseases, but most go undetected. The diagnostic process for SAVs begins with computational tools which score a patient’s millions of genetic variants to identify a manageable shortlist of candidate SAVs. Next, an analyst must manually evaluate each candidate to determine whether it warrants follow-up via an in vitro assay- a low throughput and costly process necessary for patient diagnosis. This investigation requires the analyst to assess:

The likelihood of disruption to splicing based on prediction scores
The likely impact at the mRNA and protein levels
Whether the predicted disruption matches the patient’s phenotype
Feasibility of in vitro testing methods (e.g. minigene or RT-PCR)
Tissue-specific expression of the gene and splice isoforms
Although much of this information is readily available, the analyst’s task is time-consuming and requires extensive domain expertise. Many researchers do not pursue candidate SAVs due to a lack of confidence navigating this complex analysis process, which we believe is highly amenable to assistance via a Large Language Model (LLM) agent. 

The proposed project is to build a software tool in which an LLM receives a query genetic variant, gathers relevant data using existing APIs and datasets, and provides recommendations. The system would guide users through SAV evaluation, helping reduce the expertise barrier and streamline the diagnostic process. This project would enable students to gain significant experience with a rapidly emerging and important technology: developing LLM agentic systems.
Credits/Copyright: Dao M, Venkateswaran G, Rogers C, Zhao K

# Installation
The required packages are OpenAI API, (pip install openai).
Ensure the file chatSav.py is executable and run using ./chatSav.py
# Usage
Input: Coordinates of a variant (chr:pos:ref:alt format)

VARIANT FORMAT INFORMATION:
- Chromosome: chr
- Position: pos
- Reference allele: ref
- Alternate allele: alt

