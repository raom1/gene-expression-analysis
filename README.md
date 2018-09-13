# gene-expression-analysis
Analysis of gene expression changes over time and in different cell types

### Summary
Ongoing advances in <a href="https://en.wikipedia.org/wiki/DNA_sequencing#High-throughput_methods" target="_blank">High Throughput Sequencing</a> has greatly expanded the data generated in biology. All published data must be made available for
use by the scientific community. However, many times these data go unused because of difficulties comparing between projects
or inconsistencies with methodology. These obstables can be overcome but there is also a great opportunity in further
analyzing existing data in new ways to uncover additional insights. Here I have selected two pieces of previously published
data for analysis beyond the conclusions of the paper. The original papers can be found here: <a href="https://www.cell.com/developmental-cell/fulltext/S1534-5807(17)30873-0" target="_blank">Timecourse of mouse (and human) retinal development</a> and <a href="https://www.cell.com/fulltext/S0092-8674(15)00549-8" target="_blank">Single cell sequencing of mouse retinas</a>

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
![gene_clustering](/images/gene_clustering_2.png)

### Single cell sequencing of mouse retinas
Single cell sequencing has been a great technological advancement that provides much more granular information
