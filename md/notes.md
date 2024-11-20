# notes 
- [notes](#notes)
  - [EPIC microarrays](#epic-microarrays)
    - [Principles](#principles)
      - [type I design](#type-i-design)
      - [type II](#type-ii)
    - [output files](#output-files)
      - [`.idat files`](#idat-files)
      - [`.sdf` (sample definition files)](#sdf-sample-definition-files)
      - [Illumina manifest file `.bpm`](#illumina-manifest-file-bpm)
    - [Why Batch Effects Are a Concern in Methylation Arrays](#why-batch-effects-are-a-concern-in-methylation-arrays)
    - [QC](#qc)
      - [Beadarray](#beadarray)
        - [Staining green/red](#staining-greenred)
        - [Extension](#extension)
        - [Hybridization](#hybridization)
        - [Target removal](#target-removal)
        - [Bisulfite conversion I and II](#bisulfite-conversion-i-and-ii)
        - [Specificity I and II](#specificity-i-and-ii)
        - [Sample dependent](#sample-dependent)
      - [Normalization, R](#normalization-r)
      - [density plot](#density-plot)
        - [snp probes beta distribution](#snp-probes-beta-distribution)
      - [minfi](#minfi)
    - [Bisulfite sequencing](#bisulfite-sequencing)
      - [Bisulfite sequencing disadvantages](#bisulfite-sequencing-disadvantages)
  - [Tumor](#tumor)
    - [General Brain Tumor Terminology](#general-brain-tumor-terminology)
      - [Glioma](#glioma)
        - [Astrocytoma](#astrocytoma)
        - [Glioblastoma (GBM)](#glioblastoma-gbm)
        - [Oligodendroglioma](#oligodendroglioma)
        - [Ependymoma](#ependymoma)
      - [Medulloblastoma](#medulloblastoma)
      - [Meningioma](#meningioma)
    - [Tumor Grading and Staging](#tumor-grading-and-staging)
      - [WHO Grade](#who-grade)
      - [IDH Mutation](#idh-mutation)
      - [MGMT Promoter Methylation](#mgmt-promoter-methylation)
      - [1p/19q Codeletion](#1p19q-codeletion)
      - [Ki-67](#ki-67)
    - [Imaging and Histology](#imaging-and-histology)
      - [Enhancement](#enhancement)
      - [Necrosis](#necrosis)
      - [Infiltration](#infiltration)
      - [Peritumoral Edema](#peritumoral-edema)
      - [Mitoses](#mitoses)
    - [Molecular Markers and Genetic Terms](#molecular-markers-and-genetic-terms)
      - [EGFR Amplification](#egfr-amplification)
      - [TERT Promoter Mutation](#tert-promoter-mutation)
      - [ATRX Mutation](#atrx-mutation)
      - [H3K27M Mutation](#h3k27m-mutation)
      - [PDGFRA Amplification](#pdgfra-amplification)
  - [Analysis](#analysis)
    - [Dimentional reduction techniques](#dimentional-reduction-techniques)




## EPIC microarrays
- genome-wide methylation analysis
- ~930K unique CpGs 
- whole blood, FFPE
- bisulfite converted DNA
- Analyzed with GenomeStudio, illumina software

### Principles
1. Bisulfite treatment
2. Hybridation to probe (green or red, wether the position is methylated or not)
3. Fluorescence detection and measurement 
4. Calculate the beta value

#### type I design 
CpGs are measured using a single color, with two different probes in the same color channel providing the methylated and the unmethylated measurements. 

#### type II 
CpGs measured using a Type II design are measured using a single probe, and two different colors provide the methylated and the unmethylated measurements. Practically, this implies that on this array there is not a one-to-one correspondence between probes and CpG positions.   
Only tolerate 3 cpgs in the 50 bp probe.
Type I tolerate more, but assume they are all in the same state.
 
The reasoning behind using two different probe designs is simple, but illuminating for potential developers. Type II probes use only one probe per methylation locus and hence allow more loci on the array, at a fixed array size. However, owing to the chemistry used by the type II probe design, type II probes can only tolerate up to three CpGs within the 50-bp probe. The type I design tolerates more CpGs within the 50-bp probe, but assumes that all methylation loci in the probed sequence are in the same state, i.e. the probe measuring Meth assumes all CpGs within the probed sequence to be methylated, and the probe measuring Unmeth assumes all CpGs in the probed sequence are unmethylated. The ability to tolerate more CpGs in the probe allows type I probes to be used in regions of high CpG density, such as CpG islands (CGI).   

### output files
#### `.idat files`
Raw intensity data for each probes on the array. The Green channel and Red channel files correspond to the unmethylated and methylated signal intensities, respectively.   
It stores:
- Signal intensity 
- Detection p values
- Control probe data 

#### `.sdf` (sample definition files) 
Contain sample metadata, linking each sample to its IDAT files and providing essential information for data organization and downstream analysis. 

#### Illumina manifest file `.bpm`
The manifest is a tab-delimited text (*.txt) file that specifies the names and locations of targeted reference regions. The main section of the manifest file is the Regions section and contains the following data columns:
- Column
- Description
- Name
- Unique user-specified name for the target
- Chromosome
- Chromosome location (eg, chr10, chr5, etc.)
- Start
- 1-based index for the start position of the target
- Stop
- 1-based index for the stop position of the target
- Upstream Probe Length
- The length of the upstream probe.
- Downstream Probe Length
- The length of the downstream probe. 



### [Why Batch Effects Are a Concern in Methylation Arrays](https://doi.org/10.1186/s13148-022-01277-9)
Technical Variation Across Arrays (humidity, temperature, reagent...) and time-dependent changes increase noise and reduced statistical power.
It is more challenging to control batch effect when arrays are processed one by one.

Batch effect factors include:
- laboratory environment
- operating procedures
- preanalytic sample variables
- human factors such as fastidiousness and experience

Factors of particular concern to high throughput genomic technologies include: 
- sample quality
- preservation and shipment
- the choice of methods or procedures, including nucleic acid isolation techniques and wash/clean-up conditions
- ambient conditions, including room temperature and ozone levels
- differences in scanner hardware
- lot-to-lot variance 

In the instance of microarrays in particular, batch effects often arise due to intra-batch differences in: 
- fluor labelling efficiency
- dye bias (for two colour arrays)
- dye photobleaching
- pipetting accuracy
- buffer salt concentration
- hybridisation temperature and time
- array scanner variability
- artefacts such as air bubbles in the hybridisation solution.
- latter samples are more subject to photodegradation or the damaging effects of ozone. It is known that Cy5 dye is more prone to photobleaching and ozone degradation than Cy3 dye. 

Known susceptibility factors include 
- probe GC-content
- DNA secondary structure
- Gstacking and melting temperature
- the efficiency of DNA bisulphite conversion also
needs to be considered.

What are the batch variables I could have access to? 

A subset of probes have low correlation across technical and biological replicates.

There are normalsiation methods within samples and within a batch.

PCA allows to detect batch effect, and there are specific batch removal effect software. I wonder what it'll do.

### QC 
#### Beadarray 
As it is used in the lab, it generates an excel file that summarizes controls

##### Staining green/red
The array includes staining control probes that are specifically designed to bind fluorescent dyes used in the assay.
Bad staining value could be from:
- Reagent Quality
- Improper Reagent Preparation
- Array Quality
- Machine Calibration

Independant of hybridization and extension

##### Extension

##### Hybridization
Uses synthetic targets that are present in the hybridization buffer (RA1) at 3 different concentrations.

#####  Target removal
Tests the efficiency of the stripping step after the extension reaction. It ensures any unbound or non-specifically bound DNA is removed after the hybridization step. It's a cleaning step.

##### Bisulfite conversion I and II
Type I probes are more targeted and precise, while Type II probes provide broader, general coverage of bisulfite conversion across the sample. Poor bisulfite conversion will lead to false methylation signal.

##### Specificity I and II 
Refers to the ability of the probes to bind accurately and exclusively to their intended target DNA sequence, without cross-hybridizing to non-target sequences.


##### Sample dependent

#### Normalization, R
Normalization combination of noob + BMIQ/Quantile seems to be the norm. 
See [A systematic evaluation of normalization methods and probe replicability using infinium EPIC methylation data](https://doi.org/10.1186/s13148-023-01459-z)
- SeSAMe 2 (with Noob and pOOBAH masking) was the best normalization approach, outperforming Noob + BMIQ in most metrics, in this paper.
- Quantile-based methods, including BMIQ, showed the poorest performance overall, diverging from earlier studies on 450K arrays.

#### density plot 
These plots are a visual representation of the distribution of beta values (methylation levels) across all CpG sites for each sample. Beta values for methylation typically exhibit a bimodal distribution, One peak near 0 for unmethylated CpG sites. Another peak near 1 for fully methylated CpG sites. Density plots allow you to check if your data aligns with this expected pattern. Deviations (e.g., a flattened or irregular distribution) may indicate technical issues.

So, we should see a bimodal distribution in these plots. After filtering and normalization, the picks should align better and become smoother. Any distortion would be a sign of a technical problem, or a biological issue.

##### snp probes beta distribution 
It should be trimodal or multimodal distribution : near ), near !, and near 0.5 for heterozyguous individuals. SNP probes might show outliers caused by unexpected allele frequencies, population stratification, or errors in the SNP annotation.




#### minfi
[Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays](https://pmc.ncbi.nlm.nih.gov/articles/PMC4016708/)


### Bisulfite sequencing 
#### Bisulfite sequencing disadvantages
- DNA degradation 
- conversion bias
- cost, and quite complex 
- Can't make the difference between unmethylated cytosines, and thymines. I think the only way to determine the sequence is to align the reads and see where it matches, taking into account the conversion.
- PCR bias 
  - high CG regions are harder to denaturate 
  - primers competitivity 
  - first come bias
  - exponential amplification bias
- only CpGs

## Tumor 

### General Brain Tumor Terminology
#### Glioma
A type of tumor that arises from glial cells in the brain. This category includes several subtypes based on cell origin, such as astrocytomas, oligodendrogliomas, and ependymomas.

##### Astrocytoma 
A type of glioma derived from astrocytes, star-shaped glial cells. Astrocytomas can range from low-grade (slow-growing) to high-grade (aggressive).

##### Glioblastoma (GBM)
The most aggressive and common type of primary brain tumor in adults. It’s a grade IV astrocytoma and has a poor prognosis.

##### Oligodendroglioma
A glioma derived from oligodendrocytes, which are cells involved in producing myelin in the central nervous system. These tumors often have distinct genetic markers (e.g., 1p/19q codeletion) and a better prognosis than astrocytomas.

##### Ependymoma
A glioma that arises from ependymal cells, which line the ventricles of the brain and the central canal of the spinal cord. They can occur in both children and adults.

#### Medulloblastoma
A type of brain tumor that typically arises in the cerebellum and is most common in children. It’s a type of embryonal tumor and has specific molecular subtypes.

#### Meningioma
A tumor arising from the meninges, the layers of tissue covering the brain and spinal cord. Meningiomas are usually benign but can cause symptoms based on their location and size.

### Tumor Grading and Staging
#### WHO Grade
The World Health Organization grading system used to classify the malignancy of brain tumors. Ranges from Grade I (least aggressive) to Grade IV (most aggressive).

#### IDH Mutation
A mutation in the isocitrate dehydrogenase (IDH) gene, commonly found in lower-grade gliomas. IDH-mutant tumors generally have a better prognosis than IDH-wildtype tumors.

#### MGMT Promoter Methylation
A genetic alteration that makes certain tumors, especially glioblastomas, more responsive to the chemotherapy drug temozolomide.

#### 1p/19q Codeletion
A combined deletion of chromosome arms 1p and 19q, commonly found in oligodendrogliomas. It’s associated with better prognosis and response to therapy.

#### Ki-67
A marker of cell proliferation. The Ki-67 index is often used to indicate how rapidly tumor cells are dividing; a high Ki-67 index suggests a more aggressive tumor.

### Imaging and Histology
#### Enhancement
A term referring to areas of increased contrast uptake on MRI, indicating a breakdown of the blood-brain barrier, often associated with malignancy.

#### Necrosis
Dead tissue within the tumor, commonly seen in high-grade tumors like glioblastoma. It’s a sign of rapid tumor growth outpacing its blood supply.

#### Infiltration
The spread of tumor cells into surrounding brain tissue. Infiltrative tumors are harder to completely remove surgically.

#### Peritumoral Edema
Swelling around a tumor due to fluid accumulation. Often seen in aggressive tumors and can increase intracranial pressure.

#### Mitoses
Cell divisions observed in tumor tissue under a microscope. A high mitotic index indicates rapid growth and aggressiveness.

### Molecular Markers and Genetic Terms
#### EGFR Amplification
A common genetic alteration in glioblastoma where the epidermal growth factor receptor (EGFR) gene is overexpressed. It’s associated with a more aggressive tumor phenotype.

#### TERT Promoter Mutation
A mutation in the promoter region of the telomerase reverse transcriptase (TERT) gene, often found in high-grade gliomas and associated with poor prognosis.

#### ATRX Mutation
A mutation in the ATRX gene that’s often seen in IDH-mutant, 1p/19q non-codeleted astrocytomas. Loss of ATRX expression is a diagnostic marker for certain gliomas.

#### H3K27M Mutation
A mutation in histone H3 (often H3.3 or H3.1), characteristic of diffuse midline gliomas, including brainstem gliomas in children. It’s associated with a very poor prognosis.

#### PDGFRA Amplification
Amplification of the platelet-derived growth factor receptor alpha (PDGFRA) gene, observed in some gliomas and implicated in tumorigenesis.

## Analysis 
### Dimentional reduction techniques

| Aspect                     | UMAP                                      | t-SNE                                     | PCA                                     |
|----------------------------|-------------------------------------------|-------------------------------------------|-----------------------------------------|
| **Structure Preservation** | Preserves both local and global structure fairly well | Preserves mainly local structure          | Preserves global structure well but less focused on clusters |
| **Cluster Separation**     | Often separates clusters more cleanly     | Effective at clustering but may distort distances between clusters | Does not necessarily form clear clusters |
| **Speed and Scalability**  | Faster than t-SNE for large datasets      | Slower for very large datasets            | Fastest, since it’s a linear technique  |
| **Interpretability of Distances** | Distances between clusters are more interpretable | Distances between clusters are less reliable | Distances are reliable but only along main components |
| **Parameter Tuning**       | Flexible, with key parameters like `n_neighbors` and `min_dist` for tuning | Less flexible, typically tuned with `perplexity` | No tuning parameters required |

For methylation data, UMAP seems better.
