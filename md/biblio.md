# Methylome for diagnostic - Bibliography

**Aref-Eshghi, E. et al. Evaluation of DNA Methylation Episignatures for Diagnosis and Phenotype Correlations in 42 Mendelian Neurodevelopmental Disorders.** The American Journal of Human Genetics 106, 356–370 (2020).

- The study analyzed 42 genetic syndromes to identify DNA methylation episignatures.
34 robust, disease-specific episignatures were identified.
- The study included peripheral blood DNA samples from 787 subjects with confirmed diagnoses of one of the 42 genetic syndromes.
- DNA methylation analysis was performed using Illumina Infinium methylation arrays.
- A multiclass support vector machine (SVM) classifier was developed for the concurrent assessment of the 34 episignatures.
- The study demonstrated that episignatures could reliably classify syndromes and distinguish closely related genetic disorders.
- Episignatures helped resolve ambiguous clinical cases and identified previously undiagnosed cases in a large cohort of subjects with developmental delays and congenital anomalies.
- Episignatures are stable and can be detected in peripheral blood, making them practical for clinical diagnostics
- The study’s approach more than doubles the number of genetic syndromes with known episignatures, enhancing the diagnostic capabilities for a wide range of neurodevelopmental disorders.

**Husson, T. et al. Episignatures in practice: independent evaluation of published episignatures for the molecular diagnostics of ten neurodevelopmental disorders.** Eur J Hum Genet 32, 190–199 (2024).

- It aims to validate the predictive accuracy of episignatures published in previous studies, which is crucial for their use in clinical diagnostics.

#### method 
- They used a k-nearest-neighbour (kNN) classifier within a leave-one-out scheme to provide unbiased estimates of specificity and sensitivity for each episignature.
- Sensitivity and specificity estimates were derived from a leave-one-out cross-validation approach.
- DNA was analyzed using Illumina’s Infinium EPIC array, and data was processed with the Meffil R package.
- Bisulfite conversion 
- Briefly, probes which failed methylation detection (detection p value > 0.01) in more than 5% of samples were removed. Samples with >1% of failed probes or an outlier methylation distribution (methylation/ unmethylation ≥ 3 s.d. from mean) were flagged.
- Case-control gap: To assess whether the episignature strength was reproduced on our dataset, we computed the proportion of CpG positions whose absolute average difference between cases and controls met the 5 and 10% thresholds that are typically required at discovery.
  
#### specificity and sensitivity 
- The study achieved 100% specificity for most episignatures, indicating that they could perfectly distinguish between cases and controls
- Sensitivity varied significantly across different episignatures. Some, like ATRX, DNMT3A, KMT2D, and NSD1, displayed 100% sensitivity.
- Others, such as CREBBP-RSTS and one CHD8 signature, had sensitivities below 40%, showing instability in performance.

#### discussion 
- The variability in sensitivity highlights the need for cautious use of episignatures in diagnostics.
- The study recommends larger validation sample sizes and broader testing to better characterize each episignature’s validity and reliability.
- Visual inspections of principal components and heatmaps were used to assess the reliability of case/control discrimination, supporting the quantitative findings.
- Episignatures were also tested for their ability to classify VUS. Some VUS showed methylation profiles consistent with known episignatures, but there were cases of discordance and intermediate profiles.

**Afflerbach, A.-K. et al. Nanopore sequencing from formalin-fixed paraffin-embedded specimens for copy-number profiling and methylation-based CNS tumor classification.** Acta Neuropathol 147, 74 (2024).

#### objectif 
nanopore sequencing for copy-number profiling and methylation-based classification of CNS tumors using DNA extracted from formalin-fixed paraffin-embedded (FFPE) specimens

#### background
- Nanopore sequencing offers advantages such as rapid and scalable sequencing, direct measurement of methylated cytosines, and generation of copy-number profiles but has been limited to high-quality DNA from fresh or cryopreserved samples.
- microarray allow to detect cnv and methylation at the same time

#### methodology 
- FFPE samples from 40 CNS tumors
- The study utilized two classifiers for methylation analysis: nanoDx (random forest classifier) and Sturgeon (neural network classifier).

#### findings
-  The average sequencing run produced 205,000 reads and 201 Mb of data, resulting in an average genome coverage of 0.06x.
-  NanoDx correctly classified 63% of samples. Sturgeon achieved a higher accuracy, correctly classifying 93% of samples.

#### limitations 
- Focal copy-number alterations were not reliably detectable with nanopore sequencing data.
- Both classifiers are restricted to the methylation classes available in their training sets, potentially missing some tumor types.


**Berdasco, M. & Esteller, M. Clinical epigenetics: seizing opportunities for translation.** Nat Rev Genet 20, 109–127 (2019).
- Epigenetic biomarkers offer several advantages over traditional genetic biomarkers, including their ability to reflect environmental influences and their stability in clinical samples.
- The most successful epigenetic biomarker in clinical practice is the detection of CpG methylation at specific loci.
- need for cost-effective detection methods and the complexity of interpreting dynamic epigenetic changes.
- robust methods exist 
- Epigenetic regulation is crucial for brain development and function, with abnormalities linked to neurological diseases like Alzheimer's, Parkinson's, and psychiatric disorders.
  

- epigenetics markers are advantageous for disease diagnosis, prognosis, and treatment monitoring due to their stability and ability to incorporate environmental information.

**Yu, X. et al. Cancer epigenetics: from laboratory studies and clinical trials to precision medicine.** Cell Death Discov. 10, 28 (2024).
- Epigenetic alterations play a crucial role in all stages of tumor progression, including tumorigenesis, promotion, progression, and recurrence.
- These alterations can be reversed by epigenetic drugs, making them a target for cancer therapy
  

#### mechanisms
- DNA Methylation: Involves the addition of a methyl group to cytosine residues by DNA methyltransferases (DNMTs), with DNMT1 maintaining methylation patterns and DNMT3A/3B responsible for de novo methylation.
- Histone Modifications: Includes methylation, acetylation, phosphorylation, and ubiquitylation of histones, regulated by enzymes such as lysine methyltransferases (KMTs) and lysine demethylases (KDMs).
- RNA Epigenetics: Modifications like N6-Methyladenosine (m6A) in RNA influence mRNA stability, splicing, and translation, involving writers (e.g., METTL3/4), erasers (e.g., FTO, ALKBH5), and readers (e.g., YTHDCs).

#### Epigenetic Therapeutics:
- DNMT Inhibitors (DNMTis): Azacitidine and decitabine are FDA-approved for treating myelodysplastic syndromes (MDS) and acute myeloid leukemia (AML).
- Histone Deacetylase Inhibitors (HDACis): Vorinostat, romidepsin, belinostat, and panobinostat are approved for treating T-cell lymphoma and multiple myeloma.
- Bromodomain and Extraterminal Domain (BET) Inhibitors: Target proteins like BRD4 to disrupt transcriptional regulation in cancer cells, with inhibitors such as JQ1 showing promise in various cancers.
- KMT and KDM Inhibitors: Target histone methylation and demethylation, with drugs like ORY-1001 (LSD1 inhibitor) and GSK126 (EZH2 inhibitor) showing potential in clinical trials.

#### Combination Therapy:
- Combining epigenetic drugs with other therapies (e.g., DNMTi with HDACi) can overcome resistance and improve efficacy.
Examples include azacitidine with vorinostat in MDS and combining HDAC inhibitors with immunomodulatory drugs like bortezomib and dexamethasone in multiple myeloma.
Precision Medicine:


Multi-omics: Integration of genomics, transcriptomics, and epigenomics provides comprehensive insights for personalized treatment.


Artificial Intelligence (AI): Enhances precision medicine by automating data analysis, improving early cancer diagnosis, and personalizing clinical care.


#### Future Directions:

Continued advancements in sequencing technologies and AI will further integrate epigenomic indicators into clinical applications, advancing precision diagnostics and therapeutics.
