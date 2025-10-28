<h1 align="center"><em>Tribolium castaneum</em> transcriptomics pipeline developed for the manuscript "<em>Gene expression evolution is predicted by stronger indirect selection at more pleiotropic genes</em>"</h1>
<p align="center">
<img src="https://github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression/assets/25666459/c67546e7-e749-430c-8500-6670355d8478">
<br/>
</p>

<p align="justify">
This repository contains the complete transcriptomics pipeline, as well as post-analysis scripts, developed for the manuscript <strong><em>Gene expression evolution is predicted by stronger indirect selection at more pleiotropic genes</em></strong> (see <a href="https://doi.org/10.1093/evlett/qraf039">Evolution Letters, 2025</a>). The pipeline addresses various facets of the transcriptomics data obtained from the <em>Tribolium castaneum</em> laboratory adaptive experiment published in <a href="https://doi.org/10.1371/journal.pgen.1008768">Koch & Guillaume (2020a)</a>, <a href="https://doi.org/10.1111/mec.15607">Koch & Guillaume (2020b)</a> and <a href="https://doi.org/10.1111/evo.14119">Koch et al. (2020)</a>.
</p>

<p align="justify">
The pipeline was first designed to be deployed on the <a href="https://www.csc.fi/en/">CSC computing farm</a>. <strong>For this reason, parts of the code have been anonymized</strong> (<em>e.g.</em> to connect to the CSC servers).
</p>

<p align="justify">
(<strong>💡 Tip:</strong> if mermaid flowcharts do not display, please refresh the webpage).
</p>

# Table of contents
- [Authors](#authors)
- [Citation](#citation)
- [Usage](#usage)
- [Overview](#overview)
- [Dependencies](#dependencies)
- [Description of the pipeline](#pipeline_description)
  - [Scripts](#scripts)
    - [Reorganizing BAM files](#scripts_1)
    - [Variants detection](#scripts_2)
    - [Selecting populations](#scripts_3)
    - [Genotype imputation](#scripts_4)
    - [Allelic frequency changes (AFCs)](#scripts_5)
    - [Preparing read counts](#scripts_6)
    - [Preparing phenotypes](#scripts_7)
    - [eQTLs](#scripts_8)
    - [LD decay](#scripts_9)
  - [Analyses](#analyses)
  - [Data](#data)
- [Copyright](#copyright)
- [License](#license)

# Authors <a name="authors"></a>

- Charles Rocabert (https://orcid.org/0000-0002-9224-5819)
- Eva L. Koch (https://orcid.org/0000-0001-8366-4897)
- Frédéric Guillaume (https://orcid.org/0000-0003-0874-0081)

# Citation <a name="citation"></a>

Please cite our 2025 publication: <a href="https://doi.org/10.1093/evlett/qraf039">[Evolution Letters, 2025](https://doi.org/10.1093/evlett/qraf039)</a>.

```txt
@article{10.1093/evlett/qraf039,
    author = {Koch, Eva L and Rocabert, Charles and Beeravolu Reddy, Champak and Guillaume, Frédéric},
    title = {Gene expression evolution is predicted by stronger indirect selection at more pleiotropic genes},
    journal = {Evolution Letters},
    pages = {qraf039},
    year = {2025},
    month = {10},
    abstract = {Changes in gene expression levels are central to adaptation, yet predicting and understanding their evolution remains challenging. Here, we used transcriptome-wide variation in the red flour beetle Tribolium castaneum to identify genes under selection for expression changes during adaptation to heat and drought stress and to uncover the mechanisms driving these changes. We found that estimates of genetic selection on expression levels were predictive of their evolutionary changes after 20 generations across seven independent selection lines. Evolution was largely caused by indirect selection acting on genetically correlated genes rather than by direct selection on individual genes. Consequently, central genes in co-expression networks experienced stronger selection and larger expression changes. Our genomic analysis revealed that selection on expression levels is associated with parallel allele frequency changes in the respective genes, especially in pleiotropic genes and those carrying expression quantitative trait loci, with stronger genetic selection corresponding to greater parallelism. Contrary to previous evidence of constrained evolution at more connected genes, adaptation was driven by selection acting disproportionately on genes central to co-expression gene networks. Overall, our results demonstrated that selection measured at the transcriptome level not only predicts future gene expression evolution but also provides mechanistic insight into the genetic architecture of adaptation.},
    issn = {2056-3744},
    doi = {10.1093/evlett/qraf039},
    url = {https://doi.org/10.1093/evlett/qraf039},
    eprint = {https://academic.oup.com/evlett/advance-article-pdf/doi/10.1093/evlett/qraf039/64957241/qraf039.pdf},
}
```

# Usage

To re-run the pipeline, please follow the steps below:

**(1)** Download [the latest version of this repository](https://github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression/archive/refs/heads/main.zip) on your local computer,

**(2)** Download the `data.tar.gz` file from [https://doi.org/10.5281/zenodo.15298510](https://doi.org/10.5281/zenodo.15298510), uncompress it and include it to the repository as following:

      ├── scripts
      ├── analyses
      ├── data
      └── README.md

# Overview

![final_pipeline](https://github.com/user-attachments/assets/0c7933f0-b0de-4ac4-bef4-f3415fc46f95)

# Dependencies <a name="dependencies"></a>

### • Software:
- PLINK 1.9,
- PLINK 2.0,
- GATK-4.2.3.0,
- samtools 1.14
- picard-toolkit
- bcftools-1.14
- snpEff-5.1
- Beagle-5.4
- subread-2.0.3
- GEMMA-0.98.5

### • Programming languages:
- Python3+
- R
- shell

### • R-packages
- tidyverse
- cowplot
- edgeR
- limma
- ggfortify
- tibble
- WGCNA

### • Python libraries
- swiftclient
- paramiko

# Description of the pipeline <a name="pipeline_description"></a>

      ├── scripts
      ├── analyses
      ├── data
      └── README.md

The pipeline is splitted in three main folders:
- `scripts`: This folder contains all the specific -omics tasks to build the data that will be used for higher-level analyses. This includes _e.g._ variants detection, eQTLs analysis, calculation of AFCs, etc.
- `analyses`: This folder contains higher-level analyses.
- `data`: This folder contains all the data produced by the pipeline. This folder is **not included** in the repository as it is too large to be handled by Github. It must be dowloaded separately.

## Scripts <a name="scripts"></a>

      └── scripts
           ├── 1_BAM_files_reorganization
           ├── 2_variant_call
           ├── 3_select_population
           ├── 4_genotype_imputation
           ├── 5_AFCs
           ├── 6_read_counts
           ├── 7_phenotypes
           ├── 8_eQTLs
           └── 9_LD_decay

Transcriptomics tasks are separated into folders and numbered for clarity.
For each task, the scripts are also numbered in the order of their execution, and are split between local (`local` folder) and HPC (`hpc` folder).

Sometimes, a shell script is also available to run all **local scripts** in the right order (see below). HPC scripts are designed for the CSC computing farm. The user must update the code before deployment on a computer farm.

### 📂 Reorganizing BAM files <a name="scripts_1"></a>

      └── scripts
           └── 1_BAM_files_reorganization
                └── local
                     ├── 1_CreateBamMap.py
                     ├── 2_SplitBamMap.R
                     ├── 3_CheckBamReadgroups.py
                     └── 4_RelocateBamFiles.py

The objective here is to parse the original BAM files source, re-organize them, make some adjustements and extract information.

Associated data folder(s): `./data/tribolium_bam`.

#### ⚙️ `1_CreateBamMap.py (local)`:
> This script parses the original BAM files folder and extract the information in the form of a table `./data/tribolium_bam/bam_map.csv`.

#### ⚙️ `2_SplitBamMap.R (local)`:
> This script splits the file `bam_map.csv` into many sample subsets, depending on the reference genome, the environment, the line or the generation. All files are stored in `./data/tribolium_bam` folder.
**Sample files will be used all along the pipeline**.

#### ⚙️ `3_CheckBamReadgroups.py (local)`:
> This script parses every BAM files to check the absence of the "read group" entry ("RG" label).

#### ⚙️ `4_RelocateBamFiles.py (local)`:
> This script relocates BAM files from the original hard-drive for further analysis.
> **Ultimately, BAM files are transfered to a distant server with an independent script**.

```mermaid
flowchart LR
subgraph local
direction LR
A[("Source<br/>BAMs")] --> B("1_CreateBamMap.py<br/>(local)")
B --> C[("BAM<br/>map")]
C --> D("2_SplitBamMap.R<br/>(local)")
D --> E[("Samples")]
E --> F("3_CheckBamReadgroups.py<br/>(local)")
A[("Source<br/>BAMs")] --> F
E --> G("4_RelocateBamFiles.py<br/>(local)")
A[("Source<br/>BAMs")] --> G
G --> H[("Ready for<br/>analysis<br/>BAMs")]
end
```

### 📂 Variants detection <a name="scripts_2"></a>

      └── scripts
           └── 2_variant_call
                └── hpc
                     ├── 1_BamPreProcessing.py
                     ├── 2_HaplotypeCaller.py
                     └── 3_JointCaller.py
                └── local
                     └── 4_SelectFilterAnnotateVariants.py
                └── variant_detection_pipeline.sh

This task performs a variant call on the entire transcriptomic data (version Tcas3.30 of <em>T. castaneum</em> genome). It uses various software and ultimately produces a filtered VCF file containing raw SNPs.
This raw SNPs VCF file is stored in `./data/tribolium_vcf/Tribolium_castaneum_ALL_Tcas3.30.vcf.gz`.

Associated data folder(s): `./data/tribolium_vcf`.

#### ⚙️ `1_BamPreProcessing.py` (HPC):
> This script pre-processes BAM files by (in this order):
> - Importing the unedited BAM file from the distant server,
> - Importing the reference genome and generate indices,
> - Adding the read group,
> - Uncompressing the BAM file to SAM,
> - Editing the SAM file to recalibrate MAPQ values (255 --> 60),
> - Compressing edited SAM file to BAM,
> - Copying a version of the BAM file to mark duplicates,
> - Generating BAI index file,
> - Handling splicing events,
> - Exporting edited BAM files to the distant server.

#### ⚙️ `2_HaplotypeCaller.py` (HPC):
> This script uses marked duplicates BAM files to run the complete pipeline for per-sample variant call by (in this order):
> - Importing the BAM file from the distant server,
> - Importing the reference genome and generate indexes,
> - Generating BAI index file,
> - Running GATK HaplotypeCaller. **This task can take several hours**,
> - Exporting GVCF files to the distant server.

#### ⚙️ `3_JointCaller.py` (HPC):
> This script runs the complete pipeline for the join call by (in this order):
> - Generating GenomicsDB database by:
>   - Importing the list of samples,
>   - Importing all GVCFs from the distant server,
>   - Importing the reference genome from the distant server and compute index files,
>   - Generating the sample map,
>   - Generating the interval list,
>   - Running GATK GenomicsDBImport,
>   - Exporting the database to the scratch.
> - Performing the joint-call cohort by:
>   - Importing the consolidated database from the scratch if needed,
>   - Running GATK GenotypeGVCFs,
>   - Exporting the joint-call file (VCF) to a local server.

#### ⚙️ `4_SelectFilterAnnotateVariants.py` (local):
> This script selects bi-allelic SNP variants from the original VCF file by (in this order):
> 
> - Selecting bi-allelic SNP variants,
> - Tagging `DP=0` genotypes as missing (`./.`),
> - Filtering out low quality SNPs,
> - Annotating SNPs,
> - Adding SNP unique identifiers.

**➡️ Local script to run with the shell script `variant_detection_pipeline.sh`**.

```mermaid
flowchart TB
subgraph HPC
direction LR
A[("Unedited<br/>BAMs")] --> B("1_BamPreProcessing.py<br/>(HPC)")
B --> C[("Edited<br/>BAMs")]
C --> D("2_HaplotypeCaller.py<br/>(HPC)")
D --> E[("GVCFs")]
E --> F("3_JointCaller.py<br/>(HPC)")
F --> G[("VCF")]
end
subgraph local
direction LR
H("variant_detection_pipeline.sh<br/>(local)") --> I[("bi-allelic SNPs<br/>VCF")]
end
HPC --> |"Download<br/>VCF"| local
```

### 📂 Selecting populations <a name="scripts_3"></a>

      └── scripts
           └── 3_select_population
                └── local
                     ├── 1_SelectPopulation.py
                     ├── 2_FilterGenotypes.py
                     └── 3_VariantsToTable.py
                ├── VCF_ALL_Beagle_pipeline.sh
                ├── VCF_CT_G1_eQTL_imputed_pipeline.sh
                ├── VCF_HD_G1_eQTL_imputed_pipeline.sh
                ├── VCF_imputed_genotypes_line_separation_pipeline.sh
                └── VCF_CT_HD_G1_LepMAP3_pipeline.sh

This task is a general function which splits population-level files (VCF or read counts) into sub-population files on user request.
**This function is used in further tasks involving sub-population analyses**.

Associated data folder(s): `./data/tribolium_snp`, `./data/tribolium_counts`.

#### ⚙️ `1_SelectPopulation.py` (local):
> This script selects a sub-population of samples based on a sample list (see section [3.1.1](#scripts_1")) by selecting the subset of samples in a VCF file or a read counts file.

#### ⚙️ `2_FilterGenotypes.py` (local):
> This script filters SNPs based on genotype coverage by (in this order):
> - Filtering genotypes based on F_MISSING,
> - Removing non-variant SNPs,
> - Removing LOWQUAL SNPs.

#### ⚙️ `3_VariantsToTable.py` (local):
> This script extracts VCF statistics into a text table.

**➡️ 5 pipelines are available for the selection of VCF sub-populations:**
- `VCF_ALL_Beagle_pipeline.sh`: Select ALL genotypes for Beagle imputation pipeline,
- `VCF_CT_G1_eQTL_imputed_pipeline.sh`: Select CT-G1 imputed genotypes for the eQTLs analysis,
- `VCF_HD_G1_eQTL_imputed_pipeline.sh`: Select HD-G1 imputed genotypes for the eQTLs analysis,
- `VCF_imputed_genotypes_line_separation_pipeline.sh`: Separate imputed genotypes into lines and generations for the AFC pipeline,
- `VCF_CT_HD_G1_LepMAP3_pipeline.sh`: Select CT/HD-G1 genotypes for the lep-MAP3 pipeline.

```mermaid
flowchart LR
subgraph local
direction LR
A[("VCF")] --> B("1_SelectPopulation.py<br/>(local)")
B --> C("2_FilterGenotypes.py<br/>(local)")
C --> D("3_VariantsToTable.py<br/>(local)")
D --> E[("Sub-population<br/>VCF")]
end
```

**OR**

```mermaid
flowchart LR
subgraph local
direction LR
A[("Read<br/>counts")] --> B("1_SelectPopulation.py<br/>(local)")
B --> C[("Sub-population<br/>read<br/>counts")]
end
```

### 📂 Genotype imputation <a name="scripts_4"></a>

      └── scripts
           └── 4_genotype_imputation
                └── local
                     ├── 1_ExtractNoMissingMarkers.py
                     ├── 2_ImputationTests.py
                     └── 3_ImputeGenotypes.py
                └── VCF_ALL_imputation_pipeline.sh

This task imputes missing genotypes on the global VCF file (including all samples) using Beagle software.
It is also possible to run a test to evaluate the quality of genotype imputations.

Associated data folder(s): `./data/tribolium_snp`.

#### ⚙️ `1_ExtractNoMissingMarkers.py` (local):
> This script extracts from the raw VCF file all the markers with a 100% call rate (_i.e_ without missing genotypes).
> The dataset is saved in binary format in `./data/tribolium_snp/imputation_tests`.
> 
#### ⚙️ `2_ImputationTests.py` (local):
> This script tests imputation capabilities of Beagle with toy datasets.
> Test results are saved in `./data/tribolium_snp/imputation_tests`.

#### ⚙️ `3_ImputeGenotypes.py` (local):
> Imputes genotypes of the VCF file with Beagle.

**➡️ Local script to run with the shell script `VCF_ALL_imputation_pipeline.sh`**.

```mermaid
flowchart LR
subgraph "local (benchmark)"
direction LR
B("1_ExtractNoMissingMarkers.py<br/>(local)") --> C("2_ImputationTests.py<br/>(local)")
C --> D[("Imputation<br/>success rate")]
end
subgraph "local (main pipeline)"
direction LR
E("2_ImputeGenotypes.py<br/>(local)") --> F[("Imputed<br/>VCF")]
end
A[("VCF")] --> B
A --> E
```

### 📂 Allelic frequency changes (AFCs) <a name="scripts_5"></a>

      └── scripts
           └── 5_AFCs
                └── local
                     ├── 1_MergeAFs.R
                     ├── 2_ComputeAFCs.R
                     ├── 3_SplitAFCs.R
                     └── 4_PlotAFCDistributions.R

This task computes AFCs per line and environment and splits the resulting dataset by line for further analyses.

Associated data folder(s): `./data/tribolium_afc`.

#### ⚙️ `1_MergeAFs.R` (local):
> This script merges and saves allelic frequencies.

#### ⚙️ `2_ComputeAFCs.R` (local):
> This script computes AFCs and saves the result.

#### ⚙️ `3_SplitAFCs.R` (local):
> This script Split AFCs by line and environment and saves the result.

#### ⚙️ `4_PlotAFCDistributions.R` (local):
> Optional script to display AFCs distibutions.

```mermaid
flowchart LR
subgraph "local"
direction LR
A[("SNPs<br/>table")] --> B("1_MergeAFs.R<br/>(local)")
B --> C("2_ComputeAFCs.R<br/>(local)")
C --> D("3_SplitAFCs.R<br/>(local)")
D --> E[("AFCs<br/>per line<br/>and environment")]
end
```

### 📂 Preparing read counts <a name="scripts_6"></a>

      └── scripts
           └── 6_read_counts
                └── hpc
                     ├── 1_FeatureCounts.py
                     └── 2_MergeFeatureCounts.py
                └── local
                     ├── 3_TransformReadCounts.R
                     ├── 4_DetectLowExpressedTranscripts.R
                     └── 5_StandardizeReadCounts.R
                └── read_counts_preparation_pipeline.sh

This task handles BAM files (without marked duplicates) to perform a reads count. The final output is a text file containing read counts for every called samples.

Associated data folder(s): `./data/tribolium_counts`.

#### ⚙️ `1_FeatureCounts.py` (HPC):
> This script calculates feature counts for every individual samples by (in this order):
> - Importing `subread` package and compile it,
> - Importing reference genome annotation,
> - Importing the list of samples,
> - For each BAM file, if the read counts file does not exist:
>   - Running `subread FeatureCounts`,
>   - Exporting resulting files to the distant server.

#### ⚙️ `2_MergeFeatureCounts.py` (HPC):
> This script merges feature counts from every individual samples by (in this order):
> - Importing the list of samples and all read counts,
> - Mergeing read counts in a single file,
> - Exporting the resulting file to the distant server.

#### ⚙️ `3_TransformReadCounts.R` (local):
> This script transforms a gene expression dataset by filtering it, calculating the TMM normalization, and removing run, batch and line effects.

#### ⚙️ `4_DetectLowExpressedTranscripts.R` (local):
> This script detects low expressed transcripts and save the list of expressed ones.

#### ⚙️ `5_StandardizeReadCounts.R` (local):
> This script standardizes a gene expression dataset by quantile normalization.

**⚙➡️ Local scripts to run with the shell script `read_counts_preparation_pipeline.sh`**.

```mermaid
flowchart TB
subgraph HPC
direction LR
A[("BAMs")] --> B("1_FeatureCounts.py<br/>(HPC)")
B --> C[("Per-sample<br/>read counts")]
C --> D("2_MergeFeatureCounts.py<br/>(HPC)")
D --> E[("Global<br/>read counts")]
end
subgraph local
direction LR
F("3_TransformReadCounts.R<br/>(local)") --> H[("Transformed<br/>read counts")]
G("4_DetectLowExpressedTranscripts.R<br/>(local)") --> I[("Expressed<br/>reads")]
H --> J("5_StandardizeReadCounts.R<br/>(local)")
I --> J
J --> K[("Standardized<br/>read counts")]
end
HPC --> |"Download<br/>read counts"| local
```

### 📂 Preparing phenotypes <a name="scripts_7"></a>

      └── scripts
           └── 7_phenotypes
                └── local
                     ├── 1_CopyExpressionFiles.sh
                     └── 2_ComputeRelativeFitness.R

This task uses ready-to-use gene expression data to compute various phenotypes, including relative fitness.
**Calculated phenotypes will be used in further analyses, _e.g._ eQTLs**.

Associated data folder(s): `./data/tribolium_phenotypes`.

#### ⚙️ `1_CopyExpressionFiles.sh` (local):
> This script copies ready-to-use expression files to the folder `./data/tribolium_phenotypes/`.

#### ⚙️ `2_ComputeRelativeFitness.R` (local):
> This script computes relative fitnesses per line for CT and HD individuals (G1).

```mermaid
flowchart LR
subgraph "local"
direction LR
A[("read<br/>counts")] --> B("1_CopyExpressionFiles.sh<br/>(local)")
B --> C[("Normalized<br/>gene expression")]
G[("Fitness<br/>measurements")] --> H("4_ComputeRelativeFitness.R<br/>(local)")
H --> I[("Relative<br/>fitnesses")]
end
```

### 📂 eQTLs <a name="scripts_8"></a>

      └── scripts
           └── 8_eQTLs
                └── hpc
                     ├── CheckGemmaFiles.py
                     ├── DeleteGemmaFiles.py
                     ├── 3_Gemma.py
                     └── 4_CollectSignificantEQTLs.R
                └── local
                     ├── EditFam.R
                     ├── ExtractGenePos.py
                     ├── MergeEQTLsDatasets.R
                     ├── ComputeCorrelationToFitness.R
                     ├── 1_eQTLs_data_preparation_pipeline.sh
                     ├── 2_UploadGemmaFiles.py
                     ├── 5_DownloadGemmaFiles.py
                     └── 6_merge_eQTLs_pipeline.sh
                └── phenotype_preparation_pipeline.sh

This task runs GWAAs on various phenotypes using GEMMA software to detect significant eQTLs at genome-scale.

Associated data folder(s): `./data/tribolium_eqtl`.

#### ⚙️ `1_eQTLs_data_preparation_pipeline.sh` (local):
> This script generates input files necessary to the GWAA with GEMMA software:
> - `.bed` and `.bim` files using `plink2`,
> - `.fam` and `.pheno` files using the script `EditFam.R`.
> Files are generated at once for all the phenotypes (expression and fitness). 

#### ⚙️ `2_UploadGemmaFiles.py` (local):
> This script exports all the input files (`.bed`, `.bim`, `.fam` and `.pheno`) to the distant storage server.

#### ⚙️ `CheckGemmaFiles.py` (HPC maintenance script):
> This maintenance script checks if GEMMA output files are missing after a HPC job.

#### ⚙️ `DeleteGemmaFiles.py` (HPC maintenance script):
> This maintenance script delete GEMMA output files to prevent duplicates before a new run.

#### ⚙️ `3_Gemma.py` (HPC):
> This script runs GWAAs on a given number of phenotypes by (in this order):
> - Importing GEMMA software and datasets from the distant server,
> - Calculating the kinship matrix,
> - Running GEMMA on the right set of phenotypes,
> - Converting the output to RDS format,
> - Exporting resulting files to the distant server.

#### ⚙️ `4_CollectSignificantEQTLs.R` (HPC independent script):
> This script collect significant eQTLs in a dedicated file. Two steps are applied to filter significant eQTLs:
> Collect all significant eQTL associations given a p-value threshold.
> - Genomic correction is applied and the FDR is calculated,
> - eQTLs with a p-value < 0.05 are selected.

#### ⚙️ `5_DownloadGemmaFiles.py` (local):
> This script downloads significant eQTL files from the distant server, based on a list of phenotypes.

#### ⚙️ `6_merge_eQTLs_pipeline.sh (local):
> This pipeline collect and merge significant eQTLs with additional information such as gene position (`ExtractGenePos.py` script) and gene and phenotype annotations (`MergeEQTLsDatasets.R` script).
> All phenotypes (expression and fitness) are treated at once.
> This script also calculates the correlation to fitness of every phenotypes (`ComputeCorrelationToFitness.R` script).

```mermaid
flowchart TB
subgraph sg1["local (1)"]
direction LR
A[("Selected<br/>phenotypes")] --> B("1_eQTLs_data_preparation_pipeline.sh<br/>(local)")
C[("Reference<br/>genome")] --> B
B --> D[(".bed")]
B --> E[(".bim")]
B --> F[(".fam")]
B --> G[(".pheno")]
end
subgraph sg2["HPC"]
direction LR
H("3_Gemma.py<br/>(HPC)") --> I("4_CollectSignificantEQTLs.R<br/>(HPC)")
I --> J[("Significant<br/>eQTLs")]
end
subgraph sg3["local (2)"]
direction LR
K("6_merge_eQTLs_pipeline.sh<br/>(local)") --> L[("Annotated<br/>significant<br/>eQTLs")]
end
sg1 --> |"Upload<br/>input files<br/>(2_UploadGemmaFiles.py)"| sg2
sg2 --> |"Download<br/>output files<br/>(5_DownloadGemmaFiles.py)"| sg3
```

### 📂 LD decay <a name="scripts_9"></a>

      └── scripts
           └── 9_LD_decay
                └── local
                     └── 1_calculate_LD.sh

This task runs PLINK1.9 to calculate pairwise linkage desequilibrium correlations per chromosome.

Associated data folder(s): `./data/tribolium_ld`.

#### ⚙️ `1_calculate_LD.sh` (local):
> This script calculates pairwise LD per chromosome, based on the raw SNPs VCF file containing all individuals.
> One LD output is generated per chromosome.

```mermaid
flowchart TB
subgraph sg1["local (1)"]
direction LR
A[("VCF")] --> B("1_calculate_LD.sh<br/>(local)")
B --> C[(".ld")]
end
```

## Analyses <a name="analyses"></a>

Koch et al. analysis is located in the folder `indirect_selection_analysis`, and is organized as following:

      └── indirect_selection_analysis
           ├── scripts: contains the analysis scripts
                └── 1_BuildDatasets.R: build datasets based on the pipeline outputs
                ├── 2_ShowStatistics.R: display some statistics about the data
                └── 3_GenerateFigures.R: generate manuscript figures
           ├── data: contains all the datasets generated during the analysis
           └── plots: contains the generated figures for the manuscript

## Data <a name="data"></a>

      └── data
           ├── experiment_data: contains various gathered experimental data
           ├── tribolium_afc: contains the result of AFCs calculations
           ├── tribolium_bam: contains BAM files and sample information
           ├── tribolium_counts: contains read counts datasets
           ├── tribolium_eqtl: contains significant eQTLs for all phenotypes
           ├── tribolium_filters: contains various VCF filters
           ├── tribolium_genome: contains annotated genomes
           ├── tribolium_ld: contains linkage desequilibrium calculations per chromosome
           ├── tribolium_pedigree: contains pedigrees related to the family stucture of the population
           ├── tribolium_phenotypes: contains calculated phenotypes
           ├── tribolium_snp: contains VCF files containing bi-allelic variant SNPs
           └── tribolium_vcf: contains the raw VCF file obtained from the variant call

# Copyright <a name="copyright"></a>

Copyright © 2021-2025 Charles Rocabert, Eva L. Koch, Frédéric Guillaume. All rights reserved.

# License <a name="license"></a>
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

