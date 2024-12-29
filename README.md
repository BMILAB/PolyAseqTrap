# PolyAseqTrap
A universal tool for genome-wide identification and quantification of polyadenylation sites from different 3′ end sequencing data

About
====================
PolyAseqTrap is an R package designed to identify and quantify polyA sites from various 3′ sequencing datasets (e.g., DRS, PAT-seq, PASC-seq, 3′READS). It utilizes a polyA read prioritization strategy with detailed post-inspection to minimize false polyA site calls and accurately determine their precise locations. Notably, PolyAseqTrap incorporates a transferable, cross-species deep learning model to resolve the persistent challenge of internal priming. Furthermore, it includes a weighted density peak clustering method that considers the microheterogeneity of polyadenylation across species to define polyA site clusters (PACs). The package also provides extensive tools for annotation, validation, and visualization of polyA sites, delivering well-structured reports to facilitate seamless analysis.

* The PolyAseqTrap package consists of six main modules.

<img src="https://github.com/BMILAB/PolyAseqTrap/blob/master/img/schema.png" alt="schema" style="width:60%;"/>

**a**. Preprocessing and mapping of 3′ seq data across diverse techniques. 

**b**. Identification and correction of polyA sites. 

**c**. Removing internal priming sites by deep learning. 

**d**. Clustering polyA sites to address species-specific cleavage microheterogeneity. 

**e**. Annotation of polyA sites. 

**f**. Quantitative evaluation of polyA sites. PAS, polyA site; PAC, polyA site cluster; PAR, polyA read; IP, internal priming; UMI, unique molecular identifier.


Installing PolyAseqTrap
=============
Mandatory 
---------

* R (>=4.0.0). [4.3.1 ](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [plyr](https://cran.r-project.org/web/packages/plyr/index.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html), [BiocGenerics](https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html), [tibble](https://cran.r-project.org/web/packages/tibble/index.html), [outliers](https://cran.r-project.org/web/packages/outliers/index.html), [pbmcapply](https://cran.r-project.org/web/packages/pbmcapply/index.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [methods](https://cran.r-project.org/web/packages/R.methodsS3/index.html)

Suggested R Packages
---------
* [movAPA](https://github.com/BMILAB/movAPA), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [knitr](https://cran.r-project.org/web/packages/knitr/index.html), [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/PolyAseqTrap")
library(PolyAseqTrap)

browseURL("https://bmilab.github.io/PolyAseqTrap/doc/PolyAseqTrap_tutorial.html")

##or you can download ZIP, and then unzip
devtools::install_local("your_path_of_PolyAseqTrap-master.zip", build_vignettes = TRUE)
```

Application examples
=============
We evaluated PolyAseqTrap against existing 3' sequencing pipelines using data from 16 different 3' sequencing techniques across multiple species. This comprehensive evaluation demonstrates the effectiveness and robustness of PolyAseqTrap. In this guide, we use demo data from three species—human, mouse, and Arabidopsis to illustrate how PolyAseqTrap can be applied for unified and user-friendly polyA site identification and analysis across different types of 3' sequencing data.
The demo includes the following, please refer to the vignette ([PDF](https://github.com/BMILAB/PolyAseqTrap/blob/main/doc/PolyAseqTrap_tutorial.pdf), [HTML](https://bmilab.github.io/PolyAseqTrap/doc/PolyAseqTrap_tutorial.html)) for full details. 

**Note**: To ensure efficient distribution and maintain a lightweight structure for the PolyAseqTrap R package, the demo data previously stored in the *`inst/extdata`* directory has been relocated to the [refer branch](https://github.com/BMILAB/PolyAseqTrap/tree/refer) under the [demo_data](https://github.com/BMILAB/PolyAseqTrap/tree/refer/demo_data) directory. Additionally, the training data for **DeepIP**, which was previously stored in the *`scripts`* directory, has now been moved to the [DeepIP_train](https://github.com/BMILAB/PolyAseqTrap/tree/refer/DeepIP_train) directory in the [refer branch](https://github.com/BMILAB/PolyAseqTrap/tree/refer).
* **Preparations**
* **Identify PACs at varying confidence levels from BAM file**

```
library(PolyAseqTrap,  warn.conflicts = FALSE, quietly=TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <-  BSgenome.Hsapiens.UCSC.hg38

# load 3'UTR annotation for detecting V8 polyA site
threeUTR_path <- system.file("extdata",
                             "ThreeRegion_Homo_sapiens.Rdata",
                             package = "PolyAseqTrap")
threeUTRregion <- readRDS(threeUTR_path)

# get bam file
bam_T_file <- system.file("extdata",
                          "SRR299116_T_chr22_hg_sorted.bam",
                          package = "PolyAseqTrap")


# identify and quantify PACs
pa.hg.result <- FindPTA(bam=bam_T_file, 
        yieldSize=10^7,
        reverse=F,
        bsgenome=bsgenome,
        d=24,
        poly='T',
        adjust.chr=TRUE,
        threeUTRregion=threeUTRregion,
        cutoffCount = 5,
        ext3UTRlen =   1000 ,
        isDRS = FALSE,
        run.quantify=TRUE)
head(pa.hg.result)
```


* **Remove internal priming artifacts**
* **Mitigating Microheterogeneity in PACs**
* **Annotate PACs**
* **Summary report**
** The summary report provides an overview of PACs at different levels, including their genomics region and length distributions, signal distribution, and nucleotide frequency distribution. It helps users quickly understand the predicted results by displaying the classification of PACs into categories (C1, C2, C3) and further subdividing them into subclasses (V1 to V8) based on polyA tail length, base composition, and alignment results. Additionally, chi-square values are provided to assess the similarity between identified PACs and reference sequences.
  
```
library(PolyAseqTrap,  warn.conflicts = FALSE, quietly=TRUE)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
library(movAPA, warn.conflicts = FALSE, quietly=TRUE)
library(ggpubr)
bsgenome <-BSgenome.Hsapiens.UCSC.hg38
data(data.PACds)
faFiles=faFromPACds(data.PACds, bsgenome, what='updn', fapre='updn', 
                    up=-100, dn=100,byGrp='ftr')

## plot single nucleotide profiles for 3'UTR PACs
plotATCGforFAfile("updn.3UTR.fa", ofreq=FALSE, opdf=FALSE, 
                  refPos=101, mergePlots = TRUE)


## calculate a chi-square metric to assess the similarity between the single nucleotide profile of identified polyA sites and the reference profile
# load reference PACs profile
data(ref.data)

#This demonstrates functionality with 10 iterations and 200 random PACs. Larger parameters are recommended (e.g., iteration = 100, use.size = 5000)
chisq.result <- cal.chisq(fafile="updn.3UTR.fa",ref.data=ref.data,
                          iteration=10,use.size=200)
statistic <- as.numeric(unlist(chisq.result$statistic))  
summary(statistic )
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.219   8.878   8.957   8.997   9.261   9.589 

boxplot(statistic ,col="orange",
        ylab = "Chi-squared metric",
        border = "brown"
)

## plot distributions of PACs across different categories, subclasses, genomic features, and polyA signals
data.PACds@anno$coord_class <- factor(data.PACds@anno$coord_level,                            
          levels=c("V1","V2","V3","V4","V5","V6","V7","V8","Count"),
        labels= c("C1","C1","C2","C2","C1","C2","C2","C3","Count"))

p<- plot_summary(data=data.PACds@anno)
ggarrange(p$p1,p$p2,p$p3,p$p4,    labels = c("A", "B", "C","D"))

```


The summary report provides an overview of PACs at different levels, including their genomics region and length distributions, signal distribution, and nucleotide frequency distribution. It helps users quickly understand the predicted results by displaying the classification of PACs into categories (C1, C2, C3) and further subdividing them into subclasses (V1 to V8) based on polyA tail length, base composition, and alignment results. Additionally, chi-square values are provided to assess the similarity between identified PACs and reference sequences.

<img src="https://github.com/BMILAB/PolyAseqTrap/blob/main/img/features.png" alt="features" style="width:60%;"/>



Citation
=============
If you are using PolyAseqTrap, please cite [Wenbin Ye, Xin Cheng, and Xiaohui Wu, PolyAseqTrap: a universal tool for genome-wide identification and quantification of polyadenylation sites from different 3’ end sequencing data, 2024](https://github.com/APAexplorer/PolyAseqTrap)

