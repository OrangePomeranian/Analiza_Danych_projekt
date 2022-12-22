# Analiza_Danych_projekt

We received 4 files in .vcf format of Genetic data of healthy and diseased Holstein-Friesian cows. 

The files had more than 14.3 million records for diseased individuals and more than 13.8 mln for healthy individuals. 

The aim was to find SNP-type mutations that may have a biological basis for the development of the disease and to determine their relationship with selected parameters for each chromosome.

In this project were used libraries from Python such as Pandas, os Modules, Python multiprocessing, NumPy, SciPy, and seaborn and from R dplyr, microbenchmark, and parallel. SNPs were detected with the use of the chi2 test with Yamates correction. Additionally, we examined the Pearson correlation between our results and the length of the chromosome.

Due to the enormous data, we used parallelization in Python (makeCluster and clusterApply) and R. And in a result, we discovered that project in R appeared to be faster than in Python.
