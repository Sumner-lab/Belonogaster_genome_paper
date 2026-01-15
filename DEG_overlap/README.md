# Orthogroup Analysis: Queen vs Worker Gene Expression

This R script analyzes differential gene expression patterns across three wasp species using orthogroup information to identify conserved regulatory patterns.

## Required Input Files

### 1. Gene Lists (6 files total)
Name your files following this pattern:
- `BJ_QUEEN_UP.txt` - Belonogaster juncea queen-upregulated genes
- `BJ_WORKER_UP.txt` - Belonogaster juncea worker-upregulated genes
- `LF_QUEEN_UP.txt` - Liostenogaster flavolineata queen-upregulated genes
- `LF_WORKER_UP.txt` - Liostenogaster flavolineata worker-upregulated genes
- `PD_QUEEN_UP.txt` - Polistes dominula queen-upregulated genes
- `PD_WORKER_UP.txt` - Polistes dominula worker-upregulated genes

Each file should contain one gene ID per line (no header):
```
g6387
g11318
g6205
g9230
...
```

### 2. Orthogroups File
Name: `Orthogroups.tsv`

Tab-separated file with columns:
```
Orthogroup	Belonogaster_juncea	Liostenogaster_flavolineata	Polistes_dominula
OG0000440	g11086.t1, g11100.t1	g11203.t1, g489.t1	g272.t1
OG0000441	g12345.t1	g23456.t1, g23457.t2	g34567.t1
...
```

## Installation

Install required R packages:

```r
install.packages(c("tidyverse", "ggplot2", "RColorBrewer", "VennDiagram", "futile.logger"))

# For ComplexHeatmap (from Bioconductor)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## Usage

1. Place all input files in the same directory as the script
2. Open R or RStudio
3. Set your working directory:
   ```r
   setwd("path/to/your/files")
   ```
4. Run the script:
   ```r
   source("orthogroup_analysis.R")
   ```

## Output Files

### Figures
1. **figure1_intersections.pdf/png** - Bar chart showing orthogroup intersections between species for each caste (with significance stars: *, **, ***)
2. **figure2_heatmap.pdf** - Heatmap showing presence/absence of orthogroups across all conditions
3. **figure3_patterns.pdf/png** - Bar chart summarizing conservation patterns
4. **figure4_within_species.pdf/png** - Stacked bar chart comparing queen vs worker overlap within each species
5. **figure5_enrichment.pdf/png** - Bar chart showing fold enrichment of overlaps (observed/expected) with significance indicators
6. **figure6_volcano.pdf/png** - Scatter plot showing fold enrichment vs statistical significance for all comparisons
7. **figure7_venn_queen.pdf/png** - Venn diagram showing overlap of queen-upregulated orthogroups across the three species
8. **figure7_venn_worker.pdf/png** - Venn diagram showing overlap of worker-upregulated orthogroups across the three species

### Data Files
1. **orthogroup_patterns.csv** - Complete matrix with pattern classifications for all orthogroups
2. **overlap_significance_tests.csv** - Statistical tests (hypergeometric/Fisher's exact) for all pairwise overlaps including observed counts, expected counts, fold enrichment, and P-values
3. **conserved_queen_orthogroups.txt** - Orthogroups upregulated in queens across all 3 species
4. **conserved_worker_orthogroups.txt** - Orthogroups upregulated in workers across all 3 species
5. **mixed_regulation_orthogroups.txt** - Orthogroups showing both queen and worker upregulation
6. **intersection_counts.csv** - Detailed intersection counts for plotting

## What the Script Does

1. **Loads Data** - Reads all gene lists and the orthogroups file
2. **Maps Genes to Orthogroups** - Identifies which orthogroups contain differentially expressed genes
3. **Analyzes Patterns** - Categorizes orthogroups by conservation pattern:
   - Conserved across all 3 species (queen or worker)
   - Present in 2 species
   - Mixed regulation (queen in some species, worker in others)
   - Single species
4. **Statistical Testing** - Performs hypergeometric tests on all pairwise overlaps to determine if they're greater than expected by chance:
   - Within-species tests (Queen vs Worker for each species)
   - Between-species tests (same caste across species pairs)
   - Calculates fold enrichment (observed/expected)
   - Reports P-values with significance codes (*, **, ***)
5. **Creates Visualizations** - Generates 8 different figures showing:
   - Species-level intersections (with significance annotations)
   - Overall heatmap of presence/absence
   - Conservation pattern summary
   - Within-species queen/worker comparison
   - Fold enrichment analysis across all comparisons
   - Volcano plot (enrichment vs significance)
   - **Venn diagram of queen-upregulated orthogroups** (3-way overlap)
   - **Venn diagram of worker-upregulated orthogroups** (3-way overlap)

## Customization

### Adjusting Species Names
If your files use different species abbreviations, modify lines 24-25:
```r
species <- c("BJ", "LF", "PD")
species_names <- c("B. juncea", "L. flavolineata", "P. dominula")
```

### Adjusting Column Names in Orthogroups File
If your orthogroups file has different column names, modify line 49:
```r
orthogroup_data <- orthogroups %>%
  rowwise() %>%
  mutate(
    BJ_genes = list(extract_gene_ids(Your_Column_Name_1)),
    LF_genes = list(extract_gene_ids(Your_Column_Name_2)),
    PD_genes = list(extract_gene_ids(Your_Column_Name_3))
  )
```

### Filtering Heatmap
By default, the heatmap shows orthogroups present in at least 2 conditions (line 297):
```r
filtered_matrix <- matrix_data[rowSums(matrix_data) >= 2, ]
```
Change the `>= 2` to a different threshold if desired.

## Interpreting Results

### Statistical Testing
The script performs **hypergeometric tests** (similar to Fisher's exact test) to assess whether observed overlaps between gene sets are statistically significant. This tests the null hypothesis that overlaps occur by chance given:
- The size of each gene set
- The total number of orthogroups in the analysis

**Key metrics reported:**
- **Observed overlap**: Actual number of shared orthogroups
- **Expected overlap**: Number expected by chance = (size1 × size2) / universe_size
- **Fold enrichment**: Observed / Expected (values > 1 indicate more overlap than expected)
- **P-value**: Probability of observing this overlap (or greater) by chance
  - P < 0.05 (*): Significant
  - P < 0.01 (**): Highly significant  
  - P < 0.001 (***): Very highly significant

**Example interpretation:**
If "BJ vs LF QUEEN" shows:
- Observed: 50 orthogroups
- Expected: 20.3
- Fold enrichment: 2.46×
- P-value: 1.2e-08 (***)

This means the two species share 2.46× more queen-upregulated orthogroups than expected by chance, which is extremely statistically significant. This suggests conserved queen-biased gene expression between these species.

### Conserved Orthogroups
- **Conserved queen-up**: Genes consistently upregulated in queens across all species (candidate caste-determining genes)
- **Conserved worker-up**: Genes consistently upregulated in workers across all species
- **Mixed regulation**: Genes with opposite regulation patterns between species (interesting for understanding species-specific adaptations)

### Intersection Analysis
- Shows how many orthogroups are shared between different species combinations
- Helps identify species-specific vs. conserved regulatory programs

### Pattern Distribution
- Reveals the overall evolutionary conservation of caste-biased gene expression
- High numbers in "all 3 species" categories suggest strong selective constraint
- High numbers in "single species" categories suggest rapid evolution or species-specific adaptations

## Troubleshooting

**Error: "cannot open file 'Orthogroups.tsv'"**
- Check that the file is in your working directory
- Verify the filename matches exactly (case-sensitive)

**Error: "package 'ComplexHeatmap' not found"**
- Install from Bioconductor: `BiocManager::install("ComplexHeatmap")`

**Heatmap is too large/small**
- Adjust PDF dimensions in line 313: `pdf("figure2_heatmap.pdf", width = 8, height = 10)`

**Empty gene lists**
- Verify that gene list files have no header row
- Check that gene IDs match the format in the orthogroups file (without .t1, .t2 suffixes)

## Citation

If you use this script in your research, please cite the appropriate R packages:
- ggplot2: Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
- ComplexHeatmap: Gu Z, et al. (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
