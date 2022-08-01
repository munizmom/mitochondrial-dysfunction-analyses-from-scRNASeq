# mitochondrial-dysfunction-analyses-from-scRNASeq

This collection can be used to study in deep mitochondrila dysfunction in scRNASeq data by analysing not only mitochondrial encoded genes but also nuclear encoded genes that are known to be associated to mitochondrial dysfunction based on Habermann lab database [MitoXplorer](http://mitoxplorer.ibdm.univ-mrs.fr) that is publicly available to download for mouse and human data from the Interactome area.

## ANALYTICAL STEPS:

1. Download the database
2. Identify those Differential expressed genes per cluster (DEGs) that are known to be involved on each of the 38 mitochondrial dysfunction categories
3. Calculate the significance of each MitoXplorer category (noted as MX) as a custom defined pathway using [pathFindR](https://cran.r-project.org/web/packages/pathfindR/index.html).
4. Perform different visualization plots: 
- heatmap per MX category 
- circle plots showing the number of genes found altered linked to each MX category.
- Dotplot highlighting the most deregulated genes and having less intrasample variance for prioritization analyses.
