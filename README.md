[README.md](https://github.com/user-attachments/files/33023372/README.md)
# LEC-proliferation

R scripts for the zebrafish endothelial single-cell RNA-sequencing analyses and interactive data browser associated with **Lineage-enriched cell proliferation code underlies spatio-temporal regulation of cell cycle dynamics in lymphatic vessels** by Gloger and colleagues.

The analyses compare lymphatic endothelial cells (LECs), venous endothelial cells (VECs) and their proliferative states at 36 hours post-fertilisation (hpf), 48 hpf, 5 days post-fertilisation (dpf) and 7 dpf.

The repository contains interactive analysis scripts that start from prepared Seurat objects. Input data, local paths and the analysis section must be selected before execution. Raw-read alignment, initial quality control and integration are described in the manuscript Methods; they are not implemented as an automated FASTQ-to-results pipeline in these scripts.

## Data and interactive browser

- **Zebrafish sequencing data:** [ArrayExpress E-MTAB-16797](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-16797).
- **Interactive exploration:** [Lymphatic proliferation zebrafish browser](https://lymphatic-proliferation-zebrafish.serve.scilifelab.se/app/lymphatic-proliferation-zebrafish).
- **Related datasets reanalysed in the study:** [GSE201916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201916) (Petkova et al., 2023), [GSE216970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216970) (Chen et al., 2024), [GSE169039](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169039) (Chavkin et al., 2022) and [GSE296391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296391) (Pi et al., 2025).

The manuscript Source Data and supplementary tables contain the reported quantitative results. The public browser allows exploration without installing the local analysis environment.

## Repository files

| File | Purpose |
| --- | --- |
| `Full script prolif paper MG.R` | Main downstream analysis: cluster inspection, marker expression, cluster differential expression, GO biological-process enrichment, proliferation comparisons, volcano plots, heatmaps and dot plots. |
| `GO term analysis.R` | Additional interactive GO analysis examples, including gene-identifier conversion, over-representation analysis and ranked gene-set enrichment. |
| `functions_for_example.R` | Helper functions for differential-expression filtering, gene-list overlap, marker plotting and topGO analyses. |
| `ShinyAppKKLab.R` | ShinyCell examples for generating single-dataset and multi-dataset applications and adjusting display settings. |
| `ShinyAppKKLab script updated.R` | Multi-dataset ShinyCell generation and configuration examples for the four developmental stages. |
| `ShinyAppKKLab script updated short.R` | Package setup, input loading and launch of an already generated `shinyAppMulti` application. |
| `server BUshinyMG.R` | Saved Shiny server code that expects the corresponding generated application data files. |
| `ui BUshinyMG.R` | Saved Shiny user-interface code that expects the corresponding generated application configuration files. |
| `CITATION.cff` | Machine-readable software citation and author information. |

The Shiny scripts are alternative examples, rather than consecutive steps that all need to be run.

## Required inputs

The main analysis and Shiny generation scripts load these prepared Seurat objects:

| Filename | Developmental stage |
| --- | --- |
| `data_36h.Rds` | 36 hpf |
| `data_48h.Rds` | 48 hpf |
| `data_5days.Rds` | 5 dpf |
| `data_7days.Rds` | 7 dpf |

The initial inspection section additionally loads `seurat_integration_SCT-2.Rds` and the stage-specific objects `dataI.SCT.seurat_36h.Rds`, `dataI.SCT.seurat_48h.Rds`, `dataI.SCT.seurat_5days.Rds` and `dataI.SCT.seurat_7days.Rds`.

These binary input objects are not included in the repository. The sequencing accession and public browser above do not by themselves supply every intermediate object required by the scripts. Obtain the matching processed objects before running the relevant sections. If their download location is not provided with the data record, contact the study authors or open a repository issue to request it.

The scripts use existing cluster identities, dimensional reductions and metadata. Relevant fields include `seurat_clusters`, `Phase`, `CellCycle`, `Part` and `Timepoint`, depending on the object and analysis. Proliferation comparisons use `CellCycle` labels `G1` and `prol`; phase plots use `Phase` labels `G1`, `S` and `G2M`. Differential expression uses the `RNA` assay; the Shiny generation examples use the `SCT` assay and stored UMAP coordinates. Cluster numbers are specific to the study objects and should not be transferred to a newly clustered dataset without checking cell identities.

## Software and dependencies

The following versions are documented in the manuscript Methods. They describe the study analysis environment and are not a complete dependency lockfile or a claim that every script has been tested under both R versions.

| Software or package | Version reported in the study |
| --- | --- |
| R | 4.3.1 and 4.4.0 |
| RStudio, optional IDE | 2023.06.1 and 2024.04.1+748 |
| Seurat | 5.1.0 |
| SeuratObject | 5.0.2 |
| sctransform | 0.4.1 |
| dplyr | 1.1.4 |
| ggplot2 | 3.5.1 |
| enrichplot | 1.24.4 |
| viridis | 0.6.5 |
| pheatmap | 1.0.12 |
| ggVennDiagram | 1.5.6 |
| fgsea | 1.30.0 |
| clusterProfiler | 4.12.6 |
| org.Dr.eg.db | 3.19.1 |
| ShinyCell | 2.1.0 |
| openxlsx | 4.2.8.1 |
| DoMultiBarHeatmap | 0.1.0 |

Some entries above support upstream processing or other study analyses. The repository scripts also load or request the following packages, for which exact study versions are not recorded here:

- **Analysis and plotting:** `Matrix`, `gridExtra`, `scales`, `rlang`, `tidyverse`, `topGO`, `knitr` and `kableExtra`.
- **Shiny generation and application:** `shiny`, `shinyhelper`, `data.table`, `DT`, `magrittr`, `hdf5r`, `reticulate`, `ggdendro`, `ggrepel`, `glue`, `readr`, `RColorBrewer`, `R.utils` and `devtools`.

The `grid` package is supplied with R. Install CRAN packages through `install.packages()` and Bioconductor packages, including `clusterProfiler`, `enrichplot`, `org.Dr.eg.db`, `topGO` and `fgsea`, using a compatible Bioconductor environment. ShinyCell installation instructions are provided in [SGDDNB/ShinyCell](https://github.com/SGDDNB/ShinyCell); DoMultiBarHeatmap is distributed through [elliefewings/DoMultiBarHeatmap](https://github.com/elliefewings/DoMultiBarHeatmap).

Installing current package versions does not reproduce the recorded study environment automatically. Check package versions and compatibility with the saved Seurat objects, particularly in helper functions that access assay slots directly. Save `sessionInfo()` with any new analysis output.

## Running the analyses

1. Download or clone this repository and obtain the required processed data objects.
2. Open the relevant script in RStudio or another R editor. Set `datadir` and, for the initial inspection section, `datadir2` to the directories containing your inputs.
3. Replace the hard-coded input and output paths throughout the chosen sections with your local paths. Check `outdir`, `scan()`, `source()` and file-export calls as well as the initial data directory. The main script's call to a local `functions.R` should point to the supplied `functions_for_example.R` for the helper functions used here.
4. Load the packages and data, then run the sections needed for the selected comparison. Some sections use objects created earlier in the interactive session; the scripts should not be sourced in full before paths, inputs and dependencies have been configured.
5. Use the manuscript Methods, figure legends and Source Data to identify the relevant clusters, comparisons and thresholds. The scripts also contain exploratory examples and alternative settings.

The main analysis is organised as follows:

| Stage | Analysis and outputs |
| --- | --- |
| Cluster inspection | UMAPs, marker-expression feature plots, violin plots and cluster-level marker displays for Figs. 3 and S3. |
| Cluster differential expression | RNA-assay Wilcoxon marker tests and cluster gene lists. The `get_DEG_lists_all_clusters()` function uses `min.pct = 0.1`, `logfc.threshold = 0.25`, positive markers and an adjusted P-value filter of `< 0.01`. Outputs include `DEG_all_clusters_all_timepoints.xlsx`. |
| GO biological-process analysis | Gene-symbol conversion to Entrez IDs, removal of duplicate identifiers and `enrichGO()` analysis with Benjamini-Hochberg adjustment. The main script uses `pvalueCutoff = 0.5` for exported results; this output filter is not a significance threshold of adjusted P < 0.05. |
| Proliferation comparisons | RNA-assay Wilcoxon comparisons of `G1` versus `prol` in selected LEC populations, with tabular exports and volcano plots for Figs. 4 and S4. In these comparisons, negative log2 fold change indicates higher expression in proliferating cells. Volcano highlighting uses absolute log2 fold change > 0.5 and unadjusted P < 0.1. |
| Expression visualisation | Comparisons of LEC/VEC identity and cell-cycle states, with heatmaps and dot plots for selected gene sets. |

`GO term analysis.R` contains additional GO examples using pre-existing differential-expression tables or gene lists. Match its input filenames to the lists generated for the selected comparison; those example filenames are not automatically identical to the main script's export names. Ranked enrichment additionally requires a correctly named and ordered gene-level statistic vector.

## Building or running the Shiny application

For a new local application, use the generation sections in `ShinyAppKKLab.R` or `ShinyAppKKLab script updated.R` after loading the four stage-specific objects:

1. Create the four ShinyCell configurations using `createConfig()`.
2. Apply the desired metadata labels, colours and display settings before generating the final application files. Check the example metadata selections against each object's actual fields.
3. Generate the data files using `makeShinyFiles()` with prefixes `sc1`, `sc2`, `sc3` and `sc4`, corresponding to 36 hpf, 48 hpf, 5 dpf and 7 dpf.
4. Generate the application code using `makeShinyCodesMulti()` in the same output directory. Set the title and footnote to the appropriate study information before generation.
5. Launch the complete generated directory:

```r
shiny::runApp("shinyAppMulti")
```

The scripts contain launch calls before the generation examples; use those early calls only when the application directory already exists. The short script launches an existing application and does not generate its data files. The saved `server BUshinyMG.R` and `ui BUshinyMG.R` files also depend on generated configuration, metadata and expression files, which are not bundled in this repository. Renaming or copying those two scripts alone is insufficient to create a working application.

## Citation

If you use this code, please cite:

Gloger, M., Johansson, A., Smialowska, A. & Koltowska, K. **LEC-proliferation** [Computer software]. https://github.com/LymphaticsLab/LEC-proliferation.

Machine-readable citation information is provided in [`CITATION.cff`](CITATION.cff). When reporting a new analysis, also record the repository commit used so that its code state can be identified.
