# PD_RNAseq
Transcriptome changes in human post-mortem PD striatum and blood

## runLimma.R
This function was used for all limma, GSEA and GSVA analyses in the manuscript, with the exception of MS analysis and splines analysis. This function generates limma results for all contrasts in a constrasts matrix object and then GSEA for all contrasts. Next, the GSVA package is used to generate GSVA scores for all GO pathways, and then these scores are used for limma analysis using the same contrasts matrix.  In addition, plotting objects are generated that merge the expression (lcpm for genes) or enrichment (GSVA for pathways) that are ready for ggplotting to visualize results.

## Classifier.Rmd
The Classifier.Rmd contains the code used to generate a ML classifier to distinguish PD and CTRL brain samples (Figure 1g).  The classifier was trained using 70% of the caudate samples in the manuscript, and tested on the remaining samples and on putamen samples for the same patients.

## Splines.Rmd
Splines.Rmd contains the code for Figure 5, in which splines analysis was used to find pathways that are associated with disease duration in PD brain samples.

### If you use this code, please cite the manuscript (include citation when it's available).

### Acknowledgements
