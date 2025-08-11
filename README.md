# Introduction 
Analysis code (in R) for integrated multi-omic analysis (transcriptome + proteome + phosphoproteome) of pediatric B-cell acute lymphoblastic leukemia (B-ALL) patient samples.

In this study, we used mass spectrometry to generate global proteomic and phosphoproteomic profiles of pediatric leukemia patient samples. We generated data at diagnosis (D) and remission (N) for five patients each with Ph-like and ETV6::RUNX1 B-ALL. Full methodological details can be found in the manuscript. We also incorporated existing bulk RNAseq data for these and other patient samples generated at Children's Mercy. RNAseq data available through [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002529.v2.p1). Finally, we integrated transcriptomic, proteomic, and phosphoproteomic data using Multiple Co-Inertia Analysis (MCIA) via the [omicade4 package](https://www.bioconductor.org/packages/release/bioc/html/omicade4.html).

# Inputs
Inputs for each dataset are as follows:
- Proteomics: protein-level output from ProteomeDiscover 3.0, inputs/proteo_input.csv
- Phosphoproteomics: peptide-level (phosphorylated) output from ProteomeDiscover 3.0, inputs/phospho_input.csv
- RNAseq: raw counts from kallisto, inputs/RNAseq_counts.csv; sample table, inputs/RNAseq_sampleTable.csv; batch-corrected, log2(normalized counts), inputs/RNAseq_counts_log2normalized_bc.csv

# Analysis files
We performed separate analysis for each individual 'omics dataset, as well as MCIA. Individual 'omics analyses can be performed independently, however MCIA.R requires input from the other three analyses.
1.	proteomics.R
2.	phosphoproteomics.R
3.	RNAseq.R
4.	MCIA.R

# Contact
For questions, please contact [Irina Pushel](mailto:ipushel@cmh.edu)