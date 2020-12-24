\name{REW-ISA}
\alias{REW-ISA}
\title{REW-ISA: RNA Expression Weighted Iterative Signature Algorithm}
\description{
REW-ISA is used to find potential local functional blocks (LFB) based on MeRIP-Seq data, where sites are hyper-methylated simultaneously across the specific conditions. REW-ISA adopts RNA expression levels of each site as weights to make sites of lower expression level less significant.
}
\usage{
### Given the range and step size of row threshold and column threshold, optimize the selection of thresholds in the above range.
REWISA_result <- REWISA(FPKM_IP, FPKM_INPUT, optimization=TRUE, repeat_num,
                        thr_row_interval, row_step, thr_col_interval, col_step)

### Run REW-ISA under the selected optimal row and column threshold combination
REWISA_bicluster <- REWISA(FPKM_IP, FPKM_INPUT, optimization=FALSE, optimal_LFB_num, 
                           optimal_thr_row, optimal_thr_col)
}
\arguments{
  \item{FPKM_IP}{IP sample data of m6A epi-transcriptome.}
  \item{FPKM_INPUT}{input sample data of m6A epi-transcriptome.}
  \item{MethylationLevel}{Methylation level matrix.}
  \item{ExpressionLevel}{Expression level matrix.}
  \item{optimization}{Logical. If it is TRUE, start looking for the best threshold values.}
  \item{optimal_thr_row}{The optimal row threshold found by grid search.}
  \item{optimal_thr_col}{The optimal col threshold found by grid search.}
  \item{optimal_LFB_num}{The optimal number of LFBs found by grid search}
  \item{repeat_num}{The number of times to run REW-ISA repeatedly under each pair of threshold parameter settings.}
  \item{thr_row_interval}{Range of row threshold.}
  \item{row_step}{The step size of the row threshold within its range.}
  \item{thr_col_interval}{Range of col threshold.}
  \item{col_step}{The step size of the col threshold within its range.}
}
\value{
  \item{ASwC}{In each repeated calculation, the Average Similarity within Clusters three-dimensional array calculated for each pair of threshold combinations.}
  \item{SDwC}{In each repeated calculation, the Standard Deviation within Clusters three-dimensional array calculated for each pair of threshold combinations.}
  \item{LFB_num}{In repeated experiments, a three-dimensional array of LFB numbers generated under each pair of threshold combinations.}
  \item{ASwC_mean}{The average ASwC value of each repeated calculation result in each pair of threshold combinations.}
  \item{SDwC_mean}{The average SDwC value of each repeated calculation result in each pair of threshold combinations.}
  \item{LFB_num_mode}{Under the combination of each pair of thresholds, the mode of the number of LFB is generated.}
  \item{Return value}{Function returns a list that stores optimized threshold combinations, the number of LFBs, or specific LFBs.}
}
\references{
To use REW-ISA, please cite the following reference:
Lin Zhang, Shutao Chen, Jingyi Zhu, Jia Meng and Hui Liu. REW-ISA: unveiling local functional blocks in epi-transcriptome profiling data via an RNA expression-weighted iterative signature algorithm. BMC bioinformatics, 2020, 21(1), 1-22.
}
\author{
Shutao Chen <shutao.chen@cumt.edu.cn>, Lin Zhang, Jingyi Zhu, Jia Meng and Hui Liu.
}
\note{
Give a set of FPKM_IP and FPKM_INPUT (or MethylationLevel and ExpressionLevel) to run REW-ISA
}
\examples{
### Load data
data(Methylation_level, package = "REWISA")
data(Expression_level, package = "REWISA")

data <- as.matrix(data)
weight <- as.matrix(weight)
### Using grid search to find the optimal row threshold, column threshold and the number of LFBs
REWISA_result <- REWISA(MethylationLevel = data, ExpressionLevel = weight, 
                        optimization = TRUE, repeat_num = 40, 
                        thr_row_interval = seq(1, 3, 0.1), row_step = 0.1, 
                        thr_col_interval = seq(0.1, 1.5, 0.05), col_step = 0.05)

### The final LFBs are determined according to the optimal row threshold, column threshold and the number of LFBs.
REWISA_bicluster <- REWISA(Ratio = data, weight = weight, optimization = FALSE, optimal_LFB_num = 3, 
                           optimal_thr_row = 1.2, optimal_thr_col = 0.35)
}

