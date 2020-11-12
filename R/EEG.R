#' @title  ElectroEncephaloGraphy(EEG) Data
#' @description This data arises from a large study to examine EEG
#'  correlates of genetic predisposition to alcoholism. It contains
#'  measurements from 64 electrodes placed on the scalp sampled at 
#'  256 Hz (3.9-msec epoch) for 1 second.
#' @details There were two groups of subjects: 77 alcoholic and 45
#'  control. In the original data, each subject was exposed to either 
#'  a single stimulus (S1) or to two stimuli (S1 and S2) which were 
#'  pictures of objects chosen from the 1980 Snodgrass and Vanderwart 
#'  picture set. Here, this dataset include the average of 120 trials 
#'  under S1. It is a list including two arrays,
#'  \itemize{
#'    \item alcholic Array of \eqn{256 \times 64 \times 77}.
#'    \item control  Array of \eqn{256 \times 64 \times 45}.}
#' @source \url{https://archive.ics.uci.edu/ml/machine-learning-databases/eeg-mld/}
"EEG"
#> [1] "EEG"
