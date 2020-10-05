# CRISPR-Screen-Analysis
Scripts and Functions to analyze CRISPR screen data after alignment to the guide library. This analysis pipeline will normalize the data and create several plots for QC and data visualization. The analysis pipeline requires the plotting functions in the screen graphing file. There are additional graphing functions that are not called by the pipeline and can be used to make customized scatterplots and rankplots with labels.

```
usage: crispranalysisv1.0(libSeqFile = "tkov3_index.txt", normalizationFile = "Day0Vehicle.results", 
                          nontargeting = c("luciferase", "EGFP", "LacZ"))

libSeqFile  = File containing the guide library sequences and gene information.
normalizationFile = Screen results file that all the other screen results will be compared to.
nontargeting  = Vector containing the gene symbols for any non-targeting control sgRNA's contained 
                in the library. If there are none, then input NULL.
```
The input libSeqFile must be a 2 column tab-delimited file formatted as sgRNA UID then sgRNA sequence:
```
chr19:58864777-58864796_A1BG_+	CAAGAGAAAGACCACGAGCA
chr19:58864319-58864338_A1BG_+	GCTCAGCTGGGTCCATCCTG
chr10:52619622-52619641_A1CF_+	TGCGCTGGACCAGTGCGCGG
chr10:52575966-52575985_A1CF_-	AGTTATGTTAGGTATACCCG
chr12:8975265-8975284_A2ML1_-	ATAGGGCCAACATTCCTAGA
etc...
```

The input ".results" CRISPR data files needs to be in a space-delimited file with number of reads aligned then sgRNA information with the following format, no headers:
```
2151 chr1:156212037-156212056_BGLAP_+
1234 chr15:34634170-34634189_NOP10_-
1201 chr7:72966536-72966555_BCL7B_-
1182 chr20:1629947-1629966_SIRPG_-
1180 chrX:48457776-48457795_WDR13_-
1167 chr3:129034703-129034722_H1FX_-
1152 chr19:8468376-8468395_RAB11B_+
etc...
```
The control essential genes file is uploaded in this repository as "CoreEssentialGenes.txt" and a control non-essential genes file is uploaded as "NonEssentialGenes.txt". Both files have the following format but any tab-delimited gene file with these names can be used:
```
Gene	HGNC_ID
AARS	HGNC:20
ABCE1	HGNC:69
ABCF1	HGNC:70
ACTB	HGNC:132
etc..
```
