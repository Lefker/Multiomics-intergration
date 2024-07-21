These files contain the code for my master thesis analysis: Machine Learning in multiomics data.

Here is the abstract:

Prostate cancer, for men, has the highest chance of malignancy and is one of the leading causes of death. 
In most cases it tends to grow slowly, is small in size, has low risk and limited aggressiveness. Especially 
if confined to the prostate it is potentially curable. In modern times, there have been several attempts to 
define changes at the level of the molecule that control the onset and progression of cancer.
New technologies that have been introduced in the life sciences are generating data at an extremely high rate and volume. 
Acquiring valuable knowledge from each individual omics data category, such as transcriptomics, epigenomics, 
metabolomics and proteomics, is a challenge of bioinformatics. Using computers and in statistics as key tools, 
it is possible to process the data and draw conclusions with the ultimate goal of developing therapeutic 
and prognostic strategies in the treatment of serious diseases. By applying machine learning algorithms, 
the predictive model is evaluated to correctly classify case and control samples. The complexity of molecular 
processes in biological systems makes it difficult to make confident decisions from one type of molecular analysis 
therefore generalization of conclusions from bioinformatics analysis of a single type of omics data is in most cases impractical.
The next step is to integrate the different kinds of data into polymorphic datasets with the initial idea that 
it can provide a completer and more comprehensive picture. It is expected that there are several problems in this process.
The heterogeneity of multi-modal data of different origins can cause noise and the curse of dimensionality are some of them.
Many methods have been proposed to integrate multiomics data that largely address the challenges that arise. In this paper, 
we analysed freely harmonized transcriptomics (mRNA-Seq, miRNA-Seq) and epigenomics (methylation beta values) data from the 
NIH (National Cancer Institute) TCGA-PRAD program. Using statistical methods (normalization, Tsne, t-test, FDR) we extracted 
gene signatures for each individual molecular analysis species and evaluated the accuracy of sample classification 
using machine learning models. We then integrated the data into a polymorphic dataset and followed the same methodology to 
examine the potential improvement of the predictive models.

Here is the link (Greek language): https://apothesis.eap.gr/archive/item/183861

In order to successfully merge the three datasets, only the common samples were used in all four analyses, three for each 
individual set and one for the integrated set. In the rest of the paper, the data where RNA-Seq data are referred to will 
be written as Rseq, for the corresponding DNA Methylation data as Meth and finally for miRNA-Seq as miRseq. The Rseq, 
Meth and miRseq datasets contain 554, 553 and 551 samples, respectively, with definitions Primary solid Tumor (Tumor), 
Solid Tissue Normal (Normal) and one Metastatic, which was removed due to a small number of samples in the category. 
All samples are from patients with cancer but we do not have samples from both categories (matched samples) for all patients.
In each dataset the samples have unique barcodes due to the the different NCBI analysis they received, despite being from 
the same patient and therefore the patient identity rather than the barcode was preferred as the identifier for integration. 
In the common cancer samples a C was added to the patient identities while in the normal ones an N was added and replaced 
the barcodes in all sets. For this reason, 4 patients who gave two cancer samples were kept only one of them. 
The samples (tumor-normal) in the three sets were initially 501-52 for Rseq, 502-50 for Meth and 498-52 for miRseq, 
while their common is 493-35.

In the Meth data table there were missing values and the corresponding probes were removed. In the remaining sets, 
variables that had null values in all samples were removed. This reduced the Rseq dataset by 3300, Meth by 146000 
and miRseq by 370. Therefore the embedded dataset has a size of 528x398073.
