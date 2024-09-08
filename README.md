# 1-Process-HTSeq
* Combine `HTSeq` files (*.htseq.out) into one CSV for proccessing

# 2-Differential_gene_expression
* Uses `DESeq` to compare groups, years and interaction for differentially expressed genes
* Produces files like `*_vs_*_p0.05.txt`
    * Only genes significantly differentially expressed are saved in the file (`padj < 0.05`)
    * Each row is a Gene ID (`nbisL1-transcript-*`)
    * Each row has a `log2FoldChange` (how strong is the effect) as well as `lfcSE`, `stat`, `pvalue`, `padj`
* ? Are these Gene IDs Genes?

↓  
**One file per comparison**  
↓

# 4-Functional_annotation
* For each comparison done in `Differential_gene_expression`...
* Handles genes differentially expressed "up" (`stat > 0`) separetely from "down"
    * ? Check that's right
* Deals with each Ontology (BP MF CC) separetely 
* `topGO` calculates the "Weighted Fisher" for each GO term and adjusted P value

↓  
**One file per comparison, per direction, per ontology** (6 files per comparison)  
↓

# 5-Goslim
* For each comparison done in `Differential_gene_expression` and `Functional_annotation`...
* Uses `goSlim` group GO terms into Slims
* Handles genes "up" separetely from "down"
* Deals with each Ontology (BP MF CC) separetely 
* Outputs GO Slims with counts or GO Terms in each one
* We augment this with the list of GO Terms using `goFat`

↓  
**One file per comparison, per direction, per ontology** (6 files per comparison)

---
---

# 3-Reference_to_GO
* Processes files from the feference genome to produce a CSV mapping Gene IDs (`nbisL1-transcript-*`) to GO Terms for use by `topGO`