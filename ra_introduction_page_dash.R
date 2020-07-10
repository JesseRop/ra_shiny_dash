HTML("
     
<div>
<h1>Welcome</h1>
<p>
Here, you can explore our single-cell RNA sequencing data to investigate the phenotypic spectrum of synovial tissue macrophage during development and resolution of rheumatoid arthritis (RA). 
Analyses were performed using the Seurat pipeline and custom scripts which are available upon request. In this web application we use a random subset of 8000 cells for faster computation. Further details can be found in our paper summarized below which you should cite if you use our data.</p>
<p>

</p>

<div>
<h4><a href='https://doi.org/10.1038/s41591-020-0939-8' target='_blank'>Distinct synovial tissue macrophage subsets regulate inflammation and remission in
rheumatoid arthritis </a></h4>
<p>Stefano Alivernini, Lucy MacDonald, Aziza Elmesmari, Samuel Finlay, Barbara
Tolusso, Maria Rita Gigante, Luca Petricca, Clara Di Mario, Laura Bui, Simone Perniola,
Moustafa Attar, Marco Gessi, Anna Laura Fedele, Sabarinadh Chilaka, Domenico Somma,
Stephen Sansom, Andrew Filer, Charles McSharry, Neal L. Millar, Kristina Kirschner,
Alessandra Nerviani, Myles J. Lewis, Costantino Pitzalis, Andrew R. Clark, Gianfranco
Ferraccioli, Irina Udalova, Christopher D. Buckley, Elisa Gremese, Iain B. McInnes,
Thomas D. Otto and Mariola Kurowska-Stolarska
</p>
</div>


<div>

<h4>Abstract</h4>
<p>Immune-regulatory mechanisms of drug-free remission in rheumatoid arthritis (RA) are unknown.
We hypothesised that synovial tissue macrophages (STM), which persist in remission, contribute
to joint homeostasis. We used single-cell transcriptomics to profile 32000 STMs and identified
phenotypic changes in patients with early/active RA, treatment-refractory/active RA and RA in
sustained remission. Each clinical state was characterised by different frequencies of 9 discrete
phenotypic clusters within 4 distinct STM subpopulations with diverse homeostatic, regulatory and
inflammatory functions. This cellular atlas combined with deep-phenotypic, spatial and functional
analyses of synovial biopsy FACS-sorted STMs revealed two STM subpopulations
(MerTK<sup>pos</sup>/TREM2<sup>high</sup> and MerTK<sup>pos</sup>/LYVE1<sup>pos</sup>) with unique remission transcriptomic signatures
enriched in negative-regulators of inflammation. These STMs were potent producers of
inflammation-resolving lipid mediators and induced the repair response of synovial fibroblasts <i>in
vitro</i>. A low proportion of MerTK<sup>pos</sup> STMs in remission was predictive of flare after treatment
cessation. Therapeutic fostering of MerTK<sup>pos</sup> STM-subpopulations could therefore be a successful
strategy for RA treatment.
</p>
</div>

<div>
The results are displayed in 2 sections
<ul>
<li><a onclick = 'openTab(\"all_cluster_res\")' href='#'> Clustering</a> - Allows cluster visualization and exploration of top cluster markers.</li>
<li><a onclick = 'openTab(\"ge_vis_gene\")' href='#'> Differential expression (DE)</a> - Comparison of gene expression in macrophages from healthy, UPA, Naive RA, Resistant RA and Remission RA using
parameters specified in our <a href='https://doi.org/10.1038/s41591-020-0939-8' target='_blank'> paper </a> (i.e the first 12 Principal components and 0.5 as the clustering resolution). </li>
</ul>
</div>

<div style='border:thin solid black ; padding:30px ; margin:30px'>
<figure>  

<img src='ra_scrna_multiplot_final.png' alt='Results' style='width:100%;'>
<figcaption>scRNAseq defines heterogeneity within MerTK<sup>pos</sup> /CD206<sup>pos</sup> and
MerTK neg /CD206 neg STM populations. <b>(a)</b> UMAP of 9 STM clusters identified by scRNAseq
analysis. <b>(a-h)</b> show data from Healthy (n=4), UPA (n=4), naïve-active RA (n=5), treatment-
resistant RA (n=6) and RA in remission (n=6) in 5 independent experiments. <b>(b)</b> Heatmap of
the top 20 differentially expressed genes per cluster. Top cluster markers and the total
number of genes characterized each cluster are provided. <b>(c)</b> Violin plots represent log-
normalized expression values of STM clusters’ markers with median marked by black dots
while cluster identity by unique colour. <b>(d)</b> Relationship between clusters embedded in the
top 3 Diffusion-map Components. <b>(e)</b> Hierarchical clustering of STMs. <b>(f)</b> MerTK expression
in the 9 STM clusters. <b>(g)</b> Proposed classification of human STMs. <b>(h)</b> Split UMAP and dot
plots of relative changes in the STM clusters between groups. Significant differences (*p less than 0.05) between the given condition and at least two other conditions in Two-way ANOVA
with Tukey’s correction. Precise p-values in Supplementary Figure 3a which can be found in our <a href='https://doi.org/10.1038/s41591-020-0939-8' target='_blank'> paper </a>.</figcaption>
</figure>
</div>


</div>")