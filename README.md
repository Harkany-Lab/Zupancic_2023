# Zupancic_2023

1. System requirements

- All software dependencies and operating systems (including version numbers) are available in the methods section of the manuscript and in the sessionInfo() output of the Rmarkdown and Quarto analysis notebooks at the bottom of each analysis.
- Archived version of this repo has been tested on x86_64-pc-linux-gnu (64-bit) platform running under Ubuntu 22.04.1 LTS with R versions 4.2.1 and 4.2.2, and Python version 3.8.8
- No specific hardware is required to run the code in this repo. However, we recommend using a computer with at least 64 GB of RAM and 8 CPU cores. Moreover, we recommend using a GPU with at least 8 GB of memory to clear raw expression data matrix from ambient RNA as we did for the analysis of lateral hypothalamus dataset (Mickelsen, et al., 2019).

2. Installation guide

- Installation is available via Rstudio framework infrastructure. The workflowr pipeline is in the project root.
- Typical install time on a "normal" desktop computer with "_vanilla_" environment should take several hours.

3. Demo

- We include links to the preprocessed data that one can use to reproduce figures related to scRNA-seq data from the manuscript.
- It is expected that output images can slightly differ due to random factor native to dimension reduction technics.
- Expected run time for demo on a "normal" desktop computer is within several hours.

4. Instructions for use

- To run the demo on preprocessed data you need to render Rmarkdown and Quarto analysis notebooks via workflowr and quarto packages in Rstudio environment.
