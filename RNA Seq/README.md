# 1-Process-HTSeq
* Code for RNA sequence quality control and trimming, as well as extraction of HTSeq count data (abundance of transcripts)
* Combine `HTSeq` files (*.htseq.out) into one CSV for proccessing

# 2-Differential_gene_expression
* Uses `DESeq` to compare groups, years and interaction for differentially expressed genes
* DESeq uses Wald Tests to produce files like `*_vs_*_p0.05.txt`
    * Only genes significantly differentially expressed are saved in the file (`padj < 0.05`)
    * Each row is a transcript ID (`nbisL1-transcript-*`)
    * Each row has a `log2FoldChange` (how strong is the effect of differential expression) as well as `lfcSE`, `stat`, `pvalue`, `padj`
* Also gives a count of differentially expressed genes (via likelihood ratio test - LRT) to produce the DEG bargraph by group, temperature and their interaction 

↓  
**One file per comparison**  
↓

# 4-Functional_annotation
* groups genes to give a list of significantly enriched terms
* For each comparison done in `Differential_gene_expression`...
* Handles genes differentially expressed "up" (`stat > 0`) separetely from "down" (`stat < 0`) 
* Deals with each Ontology (BP MF CC) separately
* `topGO` calculates the "Weighted Fisher" for each GO term and adjusted P value and only significantly enriched terms are used in the following step

↓  
**One file per comparison, per direction, per ontology** (6 files per comparison)  
↓

# 5-Goslim
* For each comparison done in `Differential_gene_expression` and `Functional_annotation`...
* Uses `goSlim` to group GO terms into GO slims
* Handles genes "up" separetely from "down"
* Deals with each Ontology (BP MF CC) separately 
* Outputs GO slims with counts of GO Terms grouped in each GO slim
* We augment this with the list of GO Terms using `goWide`

↓  
**One file per comparison, per direction, per ontology** (6 files per comparison)

---
---

# 3-Reference_to_GO
* Processes files from the reference genome to produce a CSV mapping Gene IDs (`nbisL1-transcript-*`) to GO Terms for use by `topGO`
