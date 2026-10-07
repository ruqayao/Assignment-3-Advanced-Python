# Assignment-3-Advanced-Python
# Gene Expression Analysis & Exploratory Data Analysis (EDA)
## Identifying Information
* **Programmer:** Qaliya Omar
* **Language:** Python 3 (Version 3.8+)
* **Date Submitted:** October 11, 2026
* **Purpose:** This project integrates multi-omics data files (`.xlsx`, `.csv`, `.tsv`) to calculate sample-level average gene expression, identify strongly differentially expressed genes ($\vert{}\text{Fold Change}\vert{} > 5$), and generate high-resolution statistical visualizations (histograms, bar charts, heatmaps, and hierarchical clustermaps).
## Required Input Files
The script expects the following three data files in the working directory:
1. **`Gene_Expression_Data.xlsx`**: Raw expression intensity values for genomic probes across Illumina GSM samples.
2. **`Gene_Information.csv`**: Gene annotations including Probe ID, Symbol, Entrez Gene ID, Chromosome, and Cytoband location.
3. **`Sample_Information.tsv`**: Tab-separated sample metadata linking GSM identifiers to sample groups (`tumor` or `normal`) and patient IDs.
## Software & Package Requirements
* **Python Environment:** Python 3.8 or higher
* **Required Libraries:**
  * `pandas` (for tabular data manipulation and file parsing)
  * `numpy` (for numerical operations)
  * `matplotlib` (for figure composition)
  * `seaborn` (for high-level statistical data visualization)
  * `openpyxl` (backend dependency for reading Excel files with pandas)
## Instructions for Execution
1. Open your terminal or command prompt.
2. Ensure analysis_script.py, Gene_Expression_Data.xlsx, Gene_Information.csv, and Sample_Information.tsv are all located in the same directory.
3. Run the Python script:
   ```bash
   python3 analysis_script.py

You can install all dependencies via pip:
```bash
pip install pandas numpy matplotlib seaborn openpyxl
