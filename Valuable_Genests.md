# Valuable Gene set

- Transcription Factor
  - [AnimalTFDB (v3.0)](http://bioinfo.life.hust.edu.cn/AnimalTFDB/)
    - AnimalTFDB is a comprehensive database including classification and annotation of genome-wide transcription factors (TFs), and transcription cofactors in 97 animal genomes. The TFs are further classified into 73 families based on their DNA-binding domain (DBD) and cofactors are classified into 83 families and 6 categories.
  
- Surface Molecule
  - [Cell Surface Protein Atlas](http://wlab.ethz.ch/cspa/)
    - A Mass Spectrometric-Derived Cell Surface Protein Atlas
  - [cluster of differentiation molecules](https://www.genenames.org/data/genegroup/#!/group/471)
    - Gene group: CD molecules from HGNC
  
- Cell type markers/Feature genes

  - [PanglaoDB](https://panglaodb.se)  (metadata on [Github](https://github.com/oscar-franzen/PanglaoDB))  IT IS SO GREAT!

    - [Cell type gene expression markers](https://panglaodb.se/markers.html)

      > This is a list of gene expression markers are used to define cell types. Green rows indicate canonical markers (classical markers used to define the cell type).
      >
      > Summary: 8286 associations (178 cell types, 4679 gene symbols, 29 tissues); Last updated: 27/03/2020 10:44:00 CET

    - [Ubiquitousness index](https://panglaodb.se/ui.html)

      > Ubiquitousness Index (0-1). 0 indicates the gene is not expressed in any cell cluster and 1 (maximum) indicates that the gene is expressed in all cell clusters. Can be useful for spotting genes that are likely to be housekeeping genes.

    - [Publications in the scientific literature](https://panglaodb.se/papers.html) (Yes, it's not markers 😂)

      >
      > This list contains 3030 papers covering various topics of single cell sequencing. It is compiled by searching for keywords in abstracts and titles using data from ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/ and ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/. It is automatically updated daily. Rows in blue indicate papers in "high impact" journals.

  - Arterial

    - Dll4, Igfbp3, Unc5b, Gja4, Hey1, Mecom, Efnb2, Epas1, Vegfc, Cxcr4, Bmx, Gja5

      > Note arterial genes might be of spatiotemporal heterogeneity

  - Venous

    - Nr2f2, Nrp2, Aplnr

  - Hematopoietic

  - Endothelial

  - Mescenchymal

  - Cell cycle
    - G1/S
    - G2/M

- Cell-cell interaction

  - [Ligand & Receptor](https://doi.org/10.1038/ncomms8866) 

    - [2015-NC] A draft network of ligand–receptor-mediated multicellular signalling in human

      > It also gave 6 consensus classes for all genes.

  - [CellPhoneDB](https://www.cellphonedb.org/)

    - CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions.

- Angiocrine/Secreted

  - [2015-NC] A draft network of ligand–receptor-mediated multicellular signalling in human
  - Subcategory
    - Cytokine, chemokine & growth factor
      - see MSigDB
    - Signaling (or Ligand)
      - see [Ligand & Receptor](https://doi.org/10.1038/ncomms8866) 
    - Extracellular Matrix (ECM)
      - see GO term

- [Metabolism]( https://www.genome.jp/dbget-bin/www_bget?pathway+hsa01100)

  - genes included in all metabolic pathways
  
- [MSigDB genesets](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp)

  > The 25724 gene sets in the Molecular Signatures Database (MSigDB) are divided into 8 major collections, and several sub-collections.
  - [**H**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=H) (hallmark gene sets, 50 gene sets)
  - [**C1**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C1) (positional gene sets, 299 gene sets)
  - [**C2**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C2) (curated gene sets, 5529 gene sets)
    - [**CGP**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CGP) (chemical and genetic perturbations, 3297 gene sets)
    - [**CP**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP) (Canonical pathways, 2232 gene sets)
      - [**CP:BIOCARTA**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:BIOCARTA) (BioCarta gene sets, 289 gene sets)
      - [**CP:KEGG**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:KEGG) (KEGG gene sets, 186 gene sets)
      - [**CP:PID**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:PID) (PID gene sets, 196 gene sets)
      - [**CP:REACTOME**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:REACTOME) (Reactome gene sets, 1532 gene sets)
  - [**C3**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C3) (regulatory target gene sets, 3735 gene sets)
    - [**MIR**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=MIR) (microRNA targets, 2598 gene sets)
      - [**MIR:MIR_Legacy**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=MIR:MIR_Legacy) (Legacy microRNA targets, 221 gene sets)
      - [**MIR:MIRDB**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=MIR:MIRDB) (MIRDB microRNA targets, 2377 gene sets)
    - [**TFT**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=TFT) (All transcription factor targets, 1137 gene sets)
      - [**TFT:GTRD**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=TFT:GTRD) (GTRD transcription factor targets, 526 gene sets)
      - [**TFT:TFT_Legacy**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=TFT:TFT_Legacy) (Legacy transcription factor targets, 611 gene sets)
  - [**C4**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C4) (computational gene sets, 858 gene sets)
    - [**CGN**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CGN) (cancer gene neighborhoods, 427 gene sets)
    - [**CM**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CM) (cancer modules, 431 gene sets)
  - [**C5**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C5) (GO gene sets, 10192 gene sets)
    - [**BP**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=BP) (GO biological process, 7530 gene sets)
    - [**CC**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CC) (GO cellular component, 999 gene sets)
    - [**MF**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=MF) (GO molecular function, 1663 gene sets)
  - [**C6**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C6) (oncogenic signatures, 189 gene sets)
  - [**C7**](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C7) (immunologic signatures, 4872 gene sets)