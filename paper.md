---
title: 'scRNAClust: A Clustering Pipeline For mammalian brain single-cell RNA-sequencing data'
tags:
  - Python
  - R
  - Snakemake
  - RNA
authors:
  - name: Emiliano Penaloza^[first author]
    affiliation: "1" 
  - name:  Jiayue Tian
    affiliation: "1"
  - name: Bazillah Zargar
    affiliation: "1"
  - name: Fatemeh Gholizadeh
    affiliation: "1" 
  - name: Camila P. E. de Souza
    affiliation: "1" 
affiliations:
 - name: University of Western Ontario
   index: 1
date: 22 December 2021
bibliography: paper.bib
link-citations: true

references:
  - type: article-journal
    id: McKenzieWangHauberg2018
    author:
    - family: McKenzie
      given: A.T
    - family: Wang
      given: M
    - family: Hauberg
      given: M.E
    - family: et al
    issued:
      date-parts:
      - - 2018
    title: 'Brain Cell Type Specific Gene Expression and Co-expression Network Architectures'
    title-short: 
    container-title: Nature
    volume: 
    issue: 
    page: 
    DOI: 10.1038/s41598-018-27293-5
    URL: 
    language: en-GB
  - type: article-journal
    id: LouCallaway2008
    author:
    - family: Lou
      given: L
    - family: Callaway
      given: E.M
    issued:
      date-parts:
      - - 2008
        - 
        - 
    title: 'Genetic dissection of neural circuits'
    title-short: 
    container-title: Nature
    volume: 
    issue: 
    page: 
    DOI: 10.1016/j.neuron.2008.01.002
    URL: 
    language: en-GB
  - type: article-journal
    id: CiortanDefrance2021
    author:
    - family: Ciortan
      given: M
    - family: Defrance
      given: M
    issued:
      date-parts:
      - - 2021
        - 
        - 
    title: 'Contrastive self-supervised clustering of scRNA-seq data'
    title-short: 
    container-title: NCBI
    volume: 
    issue: 
    page: 
    DOI: 10.1186/s12859-021-04210-8
    URL: 
    language: en-GB
  - type: article-journal
    id: JiJi2016
    author:
    - family: Ji
      given: Zhicheng
    - family: Ji
      given: Hongkai
    issued:
      date-parts:
      - - 2016
        - 
        - 
    title: 'TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis'
    title-short: 
    container-title: Nucleic Acids Research
    volume: 44
    issue: 13
    page: 
    DOI: 10.1093/nar/gkw430
    URL: 
    language: en-GB
  - type: article-journal
    id: Gini
    author:
    - family: Lan 
      given: Jiang
    - family: Chen
      given: Huidont
    - family: Pinello
      given: Luca
    - family: Yuan
      given: geo-Cheng
    issued:
      date-parts:
      - - 2016
        - 
        - 
    title: 'Giniclust: detecting rare cell types from single-cell gene expression data with gini index'
    title-short: 
    container-title: Genome biology
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
  - type: article-journal
    id: sc3
    author:
    - family: Kiselev 
      given: Vladimir Yu
    - family: Kirshner
      given: Kristina
    - family: et al
      given:
    issued:
      date-parts:
      - - 2016
        - 
        - 
    title: 'Sc3: consensus clustering of single-cell rna-seq data'
    title-short: 
    container-title: Nature methods
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
  - type: article-journal
    id: backspin
    author:
    - family: Zeisel 
      given: Amit
    - family: Munoz-Manchado
      given: Ana B
    - family: Codeluppi
      given: Simone
    - family: Lonnerberg
      given: Peter
    - family: La Manno
      given: Gioele
    - family: Jureus
      given: Anna
    - family: Marques
      given: Sueli
    - family: Munguba
      given: Hermany
    - family: He
      given: Liqun
    - family: Betsholtz
      given: Christer
    - family: et al
      given: 
    issued:
      date-parts:
      - - 2015
        - 
        - 
    title: 'Cell types in the mouse cortex and hippocampus revealed by single-cell rna-seq'
    title-short: 
    container-title: Science
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
  - type: article-journal
    id: seurat
    author:
    - family: Satija 
      given: Rahul
    - family: Farrell
      given: Jeffrey A
    - family: Gennert
      given: David
    - family: Schier
      given: Alexander F
    - family: Regev
      given: Aviv
    issued:
      date-parts:
      - - 2015
        - 
        - 
    title: 'Spatial reconstruction of single-cell gene expression data'
    title-short: 
    container-title: Nature biotechnology
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
  - type: article-journal
    id: vmeasure
    author:
    - family: Rosenberg
      given: Andrew
    - family: Hirshberg
      given: Julia
    - family: Schier
      given: Alexander F
    issued:
      date-parts:
      - - 2007
        - 
        - 
    title: 'V-measure: A conditional entropy-based external cluster evaluation measure'
    title-short: 
    container-title: 
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
  - type: article-journal
    id: ARI
    author:
    - family: Hubert
      given: Lawrence
    - family: Arabie
      given: Phipps
    issued:
      date-parts:
      - - 1985
        - 
        - 
    title: 'Comparing parititions'
    title-short: 
    container-title: Journal of classification
    volume: 
    issue: 
    page: 
    DOI: 
    URL: 
    language: en-GB
---


# Summary 

Unsupervised clustering is the primary way of identifying cell subpopulation groups in single-cell
RNA-sequencing (RNA-seq) data sets. The presented pipeline, scRNAClust, uses Python and R and provides a lightweight interface for preprocessing, clustering, and analyzing mammalian brain single-cell RNA-seq data. scRNAClust includes four state-of-the-art single-cell RNA-seq clustering algorithms, Giniclust(@Gini), Sc3(@sc3), Seurat(@seurat) and Backspin(@backspin). This easy-to-use pipeline intakes single-cell RNA-seq data through a preprocessing
step in which data are prepared into various formats digestible by each of the stated clustering techniques. The following clustering performance evaluation metrics are implemented: cluster purity,
V-measure (@vmeasure) and ARI (@ARI). Furthermore, visualization of evaluation criteria and a visual representation of cluster groups are produced
and easily retrieved for further analysis. This pipeline simply executed through Snakemake’s
command-line interface provides users with the necessary tools to produce powerful analysis in
a swift matter. Due to this being a popular field of research, the pipeline supports the integration
of new clustering algorithms. The software is easily accessible through its GitHub repository
and can be executed on any local or remote workspace.

# Statement of Need 

The mammalian brain is a complex system composed of specialized cells that vary in morphology and functionality. Distinct cell types in the brain play different and specialized roles
in electrical signaling, metabolic coupling, axonal unsheathing, regulation of blood flow, and
immune surveillance (@LouCallaway2008). Only five primary cell types have been identified: neurons, astrocytes, oligodendrocytes, microglia, and endothelial cells, which are believed
to be responsible for the outlined functions. However, depending on the source, roughly 20-50
cell types have been identified as separate entities (@McKenzieWangHauberg2018). Over the past
few years, a series of comprehensive RNA-seq experiments in different brain cell types have
been published in humans and mice to try and more concretely outline the cell landscape of
the brain. One of the main interests in single-cell RNA sequencing is clustering, which enables the
extraction of the underlying subpopulations of various cell groups. By considering the transcriptome of each cell, single-cell RNA-seq clustering can generate a high-resolution view of gene
expression in cell populations (@JiJi2016). Often the compiled datasets present challenges as
they are often zero-inflated due to specific lab protocols, as well as they have a higher quantity
of genes than the number of cells. Unfortunately, the described characteristics are not
optimal for typical clustering algorithms. As such, specialized single-cell RNA seq clustering algorithms have been engineered to compensate for these data characteristics (@CiortanDefrance2021). This introduces the need to evaluate such methods and compare them against each other. Many pipelines and workflows have been introduced to tackle this
specific issue, but there is little work done addressing the need to evaluate methods for clustering
brain cells based on their transcriptome profiles. Generating a more concrete picture of the brain’s landscape can be crucial in solving issues related to learning, memory, and other cognitive functions. The presented software
enables researchers to address the outlined problems and more efficiently conduct research


# References

