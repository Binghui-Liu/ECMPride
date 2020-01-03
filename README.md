# ***<u>ECMPride</u>***

## What is it?

ECMPride is a flexible and scalable tool developed for predicting extracellular matrix (ECM) proteins. ECMPride can directly perform ECM prediction by taking UniProt IDs in CSV (*.csv) file
format as input. The core of ECMPride was written in [R 3.6.1 language](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/) on the [RStudio 1.1.442](https://rstudio.com/products/rstudio/download/) under Windows System. The function in ECMPride are based on [R statistical environment](https://www.r-project.org/).

## The release version

### ECMPride version v 1.2.0 (1/3/2020)

- Update the positive dataset
  - We introduced more ECMs validated by the proteomic experiment of the healthy samples
  - We excluded some ECMs which appeared only in disease samples and don't have GO annotation of ECM.
  - The total number of ECMs is changed from 478 to 521.
- Update the prediction model.

### ECMPride version v 1.1.0 (12/4/2019)

- Predictive model optimization: 99 submodels were used to solve the problem of unbalanced data sets;
- Enable prediction results be reproduced accurately by setting random number seeds before training models by random forest;
- Upload [user manual](https://github.com/Binghui-Liu/ECMPride/blob/master/User%20Manual%20for%20ECMPride.pdf)  to make it easier and clearer for users to use ECMPride;
- Fixed other bugs.

### ECMPride version v 1.0.0 (3/11/2019)

- The first release version of ECMPride.

## Hardware requirements

- 2.0 GHz CPU minimum
- 2 GB RAM minimum

## Software requirements

- Supported operating system (OS) versions (32-bit or 64-bit) 
  Windows 7
  Windows 10

- [R 3.6.1](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/) or higher (for Windows) from R project

## Note

The required R packages and their installation commands are listed below:

1. `install.packages("randomForest")`  
2. `install.packages("plyr")`  
3. `install.packages("dplyr")`  
4. `install.packages("xlsx")`  
5. `install.packages("mRMRe")`  
6. `install.packages("caret")`  
7. `install.packages("parallel")`

You should install these R packages by R 3.6.1 (not R 3.5.3 or the older version)
ahead of time. In fact, ECMPride will install these R packages itself the first
time it runs, but this approach may face some unknown errors. Therefore, we
recommend users to install these R packages before running ECMPride.

## Installation

Download link: https://github.com/Binghui-Liu/ECMPride.git.

Please see  [User Manual for ECMPride.pdf](https://github.com/Binghui-Liu/ECMPride/blob/master/User%20Manual%20for%20ECMPride.pdf) for details.

## License

 Please see the file called [license.pdf](https://github.com/Binghui-Liu/ECMPride/blob/master/license.pdf).

## Contact

For any questions involving ECMPride, please contact Mr. Binghui Liu (Email: l_binghui@163.com)
