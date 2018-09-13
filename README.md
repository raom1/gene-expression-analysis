# Gene Expression Analysis
Analysis of gene expression changes over time and in different cell types

### Summary
Ongoing advances in <a href="https://en.wikipedia.org/wiki/DNA_sequencing#High-throughput_methods" target="_blank">High Throughput Sequencing</a> has greatly expanded the data generated in biology. All published data must be made available for
use by the scientific community. However, many times these data go unused because of difficulties comparing between projects
or inconsistencies with methodology. These obstables can be overcome but there is also a great opportunity in further
analyzing existing data in new ways to uncover additional insights. Here I have selected two pieces of previously published
data for analysis beyond the conclusions of the paper. The original papers can be found here: <a href="https://www.cell.com/developmental-cell/fulltext/S1534-5807(17)30873-0" target="_blank">Timecourse of mouse (and human) retinal development</a> and <a href="https://www.cell.com/fulltext/S0092-8674(15)00549-8" target="_blank">Single cell sequencing of mouse retinas</a>.

#### Technologies used
R<br>
Tidyverse<br>
ggplot2<br>
apcluster<br>
seurat<br>

### Timecourse of mouse retinal development
In the mouse, just as in other vertebrates, development of the retina is a hightly dynamic and tightly regulated process. It
takes places over an extended period of time and involves communication between many different parts of the growing body,
rapid changes in the organization of tissues, and the production of many different cell types. Perhaps intuitively, and
supported by a large body of research, these coordinated changes rely heavily on careful regulation of gene expression.
Investigation of the importance of genes in development have typically been studied on an individual basis. High Throughput
Sequencing has opened the door to study multiple genes simultaneously at a very detailed level. The <a href="http://faculty.washington.edu/tomreh/" target="_blank">Reh lab</a> published a paper comparing the development of humans
and mice, but the data can also be used to study mouse development in greater detail.

Initial exploratory data analysis revealed that there were groups of genes that had similar gene expression changes over time.
Based on this, genes were clustered to reveal groups that change in similar ways, possibly highlighting discrete gene sets
that may work cooperatively to accomplish a given developmental goal.
![gene_expression](/images/gene_expression.png)

#### Methodology
Data were collected from the <a href="https://www.ncbi.nlm.nih.gov/geo/" target="_blank">Gene Expression Omnibus</a> and
read into R. The data contain each gene as observations and different timepoints as variables. The values represent the number
of transcripts of a gene at a given timepoint. The data were cleaned to remove extra columns and taking the mean of
replicates. Gene subsets were initially identified by applying cutoffs at different timepoints to identify genes that follow a
particular pattern. In order to make the change between each timepoint more comparable, fold change was calculated between
each adjacent timepoint. These data were then taken as an input to calculate a similarity score. Euclidian distance was
manually calculated for each gene at each time point. This function took a long time to run so was instead performed on a
subset of genes (10,000). Agglomerative clustering was then performed to cluster genes and the fold change of clusters were
visualized to identify discrete gene sets and expression pattern changes.

#### Results
Overall the clustering performed well and was able to identify discrete gene sets that all appeared to behave similarly across
the observed timecourse. When multiple clusters are compared to each other and remaining genes it can be seen that the
clustering is able to identify diverse gene expression patterns that could highlight import differences in gene function
![gene_clustering](/images/gene_clustering_3.png)

### Single cell sequencing of mouse retinas
Single cell sequencing has been a great technological advancement that provides much more granular information. The analysis
pipeline is still being optimized but there are new insights that can be gained. The <a href="http://mccarrolllab.org/" target="_blank">McCarroll Lab</a> pioneered a technique called <a href="http://mccarrolllab.org/dropseq/" target="_blank">Drop-seq</a>. With this technique, the gene expression similarities and differences between each cell can be
identified and used to cluster and identify novel cell types. The mouse retina was used as a proof of concept for the Drop-seq
technique. I recreated the original main finding of the study and then dove deeper into different cell types to try and
identify different subtypes.

#### Methodology
The R package <a href="https://satijalab.org/seurat/" target="_blank">seurat</a> was used the the analysis. The original data
was collected from the <a href="https://www.ncbi.nlm.nih.gov/geo/" target="_blank">Gene Expression Omnibus</a> and processed
similarly to the methods published in the original paper. Once the processing and analysis was done individual clusters could
be identified.
![all_retinas_tsne](/images/all_retinas_tsne.png)

Multiple clusters were identified, more than the 7 major cell types, indicating subtypes were also able to be separated.
Certain genes are known markers of cell types in the retina and these were used to identify which clusters corresponded t
what cell type. Rlbp was used for Muller Glia, Vsx2 was used for Bipolar Cells, and Opn1sw was used for cones (both short and
long wave).
![all_retinas_tsne_rlbp](/images/all_retinas_tsne_rlbp.png)
![all_retinas_tsne_vsx2](/images/all_retinas_tsne_vsx2.png)
![all_retinas_tsne_opn1sw](/images/all_retinas_tsne_opn1sw.png)

Muller Glia have previously been thought to not have any subtypes. To check this, only cells that were associated with the
Muller Glia cluster were selected from the original data, and then reprocessed through the same pipeline. The tsne plot
revealed that multiple subtypes may actually exist in that population.
![mg_tsne](/images/mg_tsne.png)

## Conclusions
So much data is currently available and more is added to the literature each day. Novel experiments are important for
extending our understanding in new directions, but interesting insights can also be gained from existing data and these data
can also be combined to approach new and interesting questions.
