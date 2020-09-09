REWISA <- function(FPKM_IP, FPKM_INPUT, Ratio, weight, optimal_thr_row, optimal_thr_col,
                   optimization, repeat_num, thr_row_interval, row_step,
                   thr_col_interval, col_step,
                   fixed_side, side, row_fix_temp, col_fix_temp){

  # optimization: Logical variables. If TRUE, turn on threshold optimization
  # repeat_num indicates the number of times to run REW-ISA repeatedly under each pair of threshold parameter settings.
  # thr_row_interval and thr_col_interval represent the selection range of row and column thresholds.
  # row_step and col_step indicate that the row and column threshold is within the step size of the selection.
  # For the setting of the line threshold, it is recommended to change it in the range of 1 to 3 in steps of 0.1.
  # For the setting of the column threshold, it is recommended to change it in steps of 0.05 within 0.1 to 1.5.
  # fixed_side: Usually set to FALSE. Set to TRUE if you want to manually fix the threshold on one side.
  # side: The value 1 indicates a fixed row threshold; the value 2 indicates a fixed column threshold.
  # row_fix_temp: Fixed row threshold value.
  # col_fix_temp: Fixed column threshold value.

  # Output of REW-ISA parameter optimization result:
  # ASwC: In each repeated calculation, the Average Similarity within Clusters three-dimensional array calculated for each pair of threshold combinations.
  # SDwC: In each repeated calculation, the Standard Deviation within Clusters three-dimensional array calculated for each pair of threshold combinations.
  # LFB_num: In repeated experiments, a three-dimensional array of LFB numbers generated under each pair of threshold combinations
  # ASwC_mean and SDwC_mean: The average value of each repeated calculation result in each pair of threshold combinations.
  # LFB_num_mode: Under the combination of each pair of thresholds, the mode of the number of LFB is generated.
  # find_TR: Optimized row threshold.
  # find_TC: Optimized col threshold.
  # LFB_number: The optimal number of LFB after optimization.

  # Import necessary packages
  library(biclust)
  library(isa2)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)

  if (missing(optimization)) {
    stop("You must specify whether to enable the parameter optimization process.
         Select TRUE or FALSE for optimization.")
  }
  if ((missing(thr_row_interval) && missing(optimal_thr_row) && missing(row_fix_temp))) {
    stop("Missing row threshold related parameters")
  }
  if ((missing(thr_col_interval) && missing(optimal_thr_col) && missing(col_fix_temp))) {
    stop("Missing col threshold related parameters")
  }
  if (optimization == TRUE && fixed_side == TRUE) {
    stop("Conflict between optimized parameters and fixed one-sided")
  }

  if ((missing(Ratio) || missing(weight))) {
    FPKM_IP <- FPKM_IP + 0.001
    FPKM_INPUT <- FPKM_INPUT + 0.001
    data_sum <- FPKM_IP + FPKM_INPUT - 0.002
    Ratio <- FPKM_IP / data_sum
    weight <- log2(data_sum + 1)
  }
  len_row <- nrow(Ratio)
  len_col <- ncol(Ratio)

  row_min <- as.numeric(apply(Ratio, 1, min))
  row_max <- as.numeric(apply(Ratio, 1, max))
  row_dis <- row_max - row_min
  row_min_matrix <- matrix(rep(row_min,len_col), ncol=len_col)
  row_dis_matrix <- matrix(rep(row_dis,len_col), ncol=len_col)
  test_1 <- (Ratio - row_min_matrix) / row_dis_matrix
  test_1 <- t(test_1)
  col_min <- apply(Ratio, 2, min)
  col_max <- apply(Ratio, 2, max)
  col_dis <- col_max - col_min
  col_min_matrix <- t(matrix(rep(col_min,len_row), ncol=len_row))
  col_dis_matrix <- t(matrix(rep(col_dis,len_row), ncol=len_row))
  test_2 <- (Ratio - col_min_matrix) / col_dis_matrix
  nm_test <- list()
  nm_test[[1]] <- test_1
  nm_test[[2]] <- test_2
  names(nm_test)[1] <- paste("Er")
  names(nm_test)[2] <- paste("Ec")
  nm_w <- list()
  nm_w[[1]] <- nm_test[[1]] * t(weight)
  nm_w[[2]] <- nm_test[[2]] * weight
  names(nm_w)[1] <- paste("Er")
  names(nm_w)[2] <- paste("Ec")
  data_weight <- Ratio * weight

  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }


  if (fixed_side == FALSE) {
    # Grid search threshold optimization
    if (optimization == TRUE) {
      if(repeat_num < 1){
        stop("The number of repetitions must be greater than 0")
      }
      if ((missing(thr_row_interval) || missing(thr_col_interval))) {
        stop("Incorrect interval, please check interval")
      }
      if ((missing(row_step) || missing(col_step))) {
        stop("Incorrect row_step or col_step, please check step!")
      }

      thr_row_initial <- thr_row_interval[1]
      thr_row_length <- length(thr_row_interval)
      thr_row_num <- thr_row_length
      thr_row_step <- row_step
      thr_row_final <- thr_row_interval[thr_row_length]

      thr_col_initial <- thr_col_interval[1]
      thr_col_length <- length(thr_col_interval)
      thr_col_num <- thr_col_length
      thr_col_step <- col_step
      thr_col_final <- thr_col_interval[thr_col_length]

      ASwC <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      SDwC <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      LFB_num <- array(dim = c(repeat_num, thr_col_num, thr_row_num))
      ASwC_mean <- array(dim = c(thr_row_num, thr_col_num))
      SDwC_mean <- array(dim = c(thr_row_num, thr_col_num))
      LFB_num_mode <- array(dim = c(thr_row_num, thr_col_num))

      for (row_index in 1:thr_row_num) {
        thr_row_temp <- thr_row_interval[row_index]
        cat("\nthr_row_position:", row_index)
        cat("\nthr_row:", thr_row_temp)

        for (cycle in 1:repeat_num) {
          cat("\nrepeat:", cycle)
          epoch <- 1
          for (thr in seq(thr_col_initial, thr_col_final, thr_col_step)) {
            cat("\nthr_col:", thr)
            ## Random seeds
            seeds <- generate.seeds(length=len_row, count=100)
            isares <- isa.iterate(nm_w, row.seeds = seeds,
                                  thr.row = thr_row_temp, thr.col = thr)
            ## Eliminate duplicated modules
            isares.unique <- isa.unique(nm_w, isares)
            ## Filter out not robust ones
            isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
            bc <- isa.biclust(isares2)

            isa_row <- bc@RowxNumber
            isa_col <- bc@NumberxCol
            isa_col_sum <- apply(isa_col, 1, sum)
            isa_drop_index <- which(isa_col_sum == 1)
            if(length(isa_drop_index) != 0) {
              isa_row_new <- as.matrix(isa_row[, -isa_drop_index])
              if(ncol(isa_row_new) == 1) {
                isa_col_new <- t(as.matrix(isa_col[-isa_drop_index, ]))
              } else {
                isa_col_new <- as.matrix(isa_col[-isa_drop_index, ])
              }
              isa_number <- nrow(isa_col_new)
              bc@RowxNumber <- as.matrix(isa_row_new)
              bc@NumberxCol <- as.matrix(isa_col_new)
              bc@Number <- isa_number
            }

            # Cluster analysis
            isa_row <- bc@RowxNumber
            isa_col <- bc@NumberxCol
            aswc_temp <- vector()
            mean <- vector()
            std <- vector()
            cluster_number <- ncol(isa_row)
            if(cluster_number < 1) {
              ASwC[cycle, epoch, row_index] <- 0
              SDwC[cycle, epoch, row_index] <- 0
            }
            if(cluster_number >= 1) {
              for (temp_num in 1:cluster_number) {
                temp_pcc_vec_col <- vector()
                temp_pcc_post_col <- 1
                r <- isa_row[, temp_num]
                c <- isa_col[temp_num, ]
                rr <- list()
                i <- 1
                for(var in 1:len_row) {
                  if(r[var]==TRUE) {
                    rr[i] <- var
                    i <- i+1
                  }
                }
                cc <- list()
                i <- 1
                for(var in 1:len_col) {
                  if(c[var]==TRUE) {
                    cc[i] <- var
                    i <- i+1
                  }
                }
                length_c <- length(cc)
                length_r <- length(rr)
                # Calculate ASwC
                if (length_c == 1) {
                  aswc_temp[temp_num] <- 0
                }
                else {
                  for (var_1 in cc) {
                    for (var_2 in cc) {
                      temp_vector_1 <- vector()
                      temp_vector_2 <- vector()
                      temp_post <- 1
                      for (var_3 in rr) {
                        temp_vector_1[temp_post] <- data_weight[var_3, var_1]
                        temp_vector_2[temp_post] <- data_weight[var_3, var_2]
                        temp_post <- temp_post + 1
                      }
                      temp_pcc_vec_col[temp_pcc_post_col] <- cor(temp_vector_1, temp_vector_2, method='pearson')
                      temp_pcc_post_col <- temp_pcc_post_col + 1
                    }
                  }
                  aswc_temp_sum <- (sum(abs(temp_pcc_vec_col)) / 2) - 0.5 * length_c
                  aswc_temp[temp_num] <- (2 / (length_c * (length_c - 1))) * aswc_temp_sum
                }

                # Calculate SDwC
                temp <- 0
                for (var_1 in rr) {
                  for (var_2 in cc) {
                    temp <- temp + data_weight[var_1, var_2]
                  }
                }
                mean_temp <- temp / (length_c * length_r)
                mean[temp_num] <- mean_temp
                count <- 0
                for (var_1 in rr) {
                  for (var_2 in cc) {
                    temp <- (data_weight[var_1, var_2] - mean_temp)^2
                    count <- count + temp
                  }
                }
                std_temp <- sqrt(count / (length_c * length_r))
                std[temp_num] <- std_temp
              }
              ASwC[cycle, epoch, row_index] <- mean(aswc_temp)
              SDwC[cycle, epoch, row_index] <- mean(std)
            }
            # Statistics cluster number
            LFB_num[cycle, epoch, row_index] = cluster_number
            epoch <- epoch + 1
          }
        }
      }
      ASwC[is.na(ASwC)] <- 0
      SDwC[is.na(SDwC)] <- 0
      # calculate the mean of ASwC, SDwC and the mode of LFB_num
      for (var_1 in 1:thr_row_num) {
        for (var_2 in 1:thr_col_num) {
          ASwC_mean[var_1, var_2] <- mean(ASwC[, var_2, var_1])
          SDwC_mean[var_1, var_2] <- mean(SDwC[, var_2, var_1])
          LFB_num_mode[var_1, var_2] <- getmode(LFB_num[, var_2, var_1])
        }
      }

      LFB_all <<- LFB_num_mode
      nrow_temp <- 3
      ncol_temp <- 5
      nrow_step <- (nrow_temp - 1) / 2
      ncol_step <- (ncol_temp - 1) / 2
      LFB_max <- which(LFB_num_mode == max(LFB_num_mode), arr.ind = T)
      row_max_index <- max(LFB_max[, 1]) + nrow_step
      col_max_index <- max(LFB_max[, 2]) + ncol_step
      LFB_num_all <- LFB_num_mode
      LFB_num_mode <- LFB_num_mode[1:row_max_index, 1:col_max_index]

      thr_row_num_second <- row_max_index
      thr_col_num_second <- col_max_index
      LFB_sliding <- matrix(1, nrow = nrow_temp, ncol = ncol_temp)
      sliding_score <- array(dim = c((thr_row_num_second - nrow_step), (thr_col_num_second - ncol_step)))
      for (var_1 in 1:(thr_row_num_second - nrow_step)) {
        for (var_2 in 1:(thr_col_num_second - ncol_step)) {
          if (var_1 <= nrow_step || var_2 <= ncol_step){
            sliding_row <- var_1 - nrow_step
            if (sliding_row < 1) {
              sliding_row <- 1
            }
            sliding_col <- var_2 - ncol_step
            if (sliding_col < 1) {
              sliding_col <- 1
            }
            LFB_select_min <- LFB_num_mode[sliding_row:(var_1 + nrow_step), sliding_col:(var_2 + ncol_step)]
          }

          if (var_1 > nrow_step && var_2 > ncol_step) {
            LFB_select_min <- LFB_num_mode[(var_1 - nrow_step) : (var_1 + nrow_step),
                                           (var_2 - ncol_step) : (var_2 + ncol_step)]
          }
          LFB_select_min <- apply(LFB_select_min, 2, as.numeric)


          LFB_min_filter <- table(as.matrix(as.data.frame(LFB_select_min)))
          LFB_min_filter <- as.data.frame(LFB_min_filter)
          LFB_min_filter[, 1] <- as.numeric(as.character(LFB_min_filter[, 1]))
          LFB_min_filter[, 2] <- as.numeric(LFB_min_filter[, 2])
          LFB_min_mode <- getmode(LFB_select_min)
          LFB_min_fre_index <- which(LFB_min_filter[, 1] == LFB_min_mode)
          LFB_min_fre <- LFB_min_filter[LFB_min_fre_index, 2]
          LFB_min_fre_all <- sum(LFB_min_filter[, 2])
          if (LFB_num_mode[var_1, var_2] == LFB_min_mode) {
            sliding_score[var_1, var_2] <- 1
          } else {
            sliding_score[var_1, var_2] <- 0
          }
        }
      }

      #####################现在的阈值寻找######################
      #####################找到最大的LFB##############
      LFB_weight <- sliding_score * LFB_num_mode[(1:(nrow(LFB_num_mode)-nrow_step)),
                                                 (1:(ncol(LFB_num_mode)-ncol_step))]
      LFB_weight_max <- max(LFB_weight)

      row_max_index <- row_max_index - nrow_step
      col_max_index <- col_max_index - ncol_step
      LFB_num_mode <- LFB_num_mode[1:row_max_index, 1:col_max_index]

      LFB_number_filter <- table(as.matrix(as.data.frame(LFB_weight)))
      LFB_number_filter <- as.data.frame(LFB_number_filter)
      delete0 <- which(LFB_number_filter[, 1] == 0)
      LFB_number_filter <- LFB_number_filter[-delete0, ]
      LFB_sum <- sum(LFB_number_filter[, 2])
      LFB_number_filter[, 1] <- as.numeric(as.character(LFB_number_filter[, 1]))
      LFB_number_filter[, 2] <- LFB_number_filter[, 2] / LFB_sum
      LFB_filter_find <- LFB_number_filter[, 1] * LFB_number_filter[, 2]
      LFB_filter_index <- which(LFB_filter_find == max(LFB_filter_find))
      LFB_filter_final <- as.numeric(LFB_number_filter[LFB_filter_index, 1])
      LFB_number <- LFB_filter_final
      LFB_weight[LFB_weight != LFB_filter_final] <- 0
      LFB_weight[LFB_weight == LFB_filter_final] <- 1

      ASwC_mean <- ASwC_mean[1:row_max_index, 1:col_max_index]
      SDwC_mean <- SDwC_mean[1:row_max_index, 1:col_max_index]
      SDwC_mean_norm <- (SDwC_mean - min(SDwC_mean)) / (max(SDwC_mean) - min(SDwC_mean))
      SDwC_mean_norm_2 <- 1 - SDwC_mean_norm
      ASwC_mean_norm <- (ASwC_mean - min(ASwC_mean)) / (max(ASwC_mean) - min(ASwC_mean))
      ASwC_mean_norm_2 <- 1 - ASwC_mean_norm

      SDwC_weight <- SDwC_mean_norm_2 * LFB_weight
      SDwC_useful_num <- length(which(SDwC_weight != 0))
      SDwC_temp <- SDwC_mean_norm_2 * LFB_weight
      SDwC_temp[SDwC_temp < (sum(SDwC_weight) / SDwC_useful_num)] <- 0
      SDwC_temp[SDwC_temp != 0] <- 1

      ASwC_temp <- ASwC_mean_norm_2 * SDwC_temp
      thr_row_find_index <- as.numeric(as.data.frame(which(ASwC_temp == max(ASwC_temp),
                                                           arr.ind = T))[1, 1])
      thr_col_find_index <- as.numeric(as.data.frame(which(ASwC_temp == max(ASwC_temp),
                                                           arr.ind = T))[1, 2])
      thr_row_find <- thr_row_interval[thr_row_find_index]
      thr_col_find <- thr_col_interval[thr_col_find_index]
      LFB_number <- LFB_num_mode[thr_row_find_index, thr_col_find_index]

      find_TR <<- thr_row_find
      find_TC <<- thr_col_find
      LFB_number <<- LFB_number
      ASwC <<- ASwC
      SDwC <<- SDwC
      LFB_num <<- LFB_num
      ASwC_mean <<- ASwC_mean
      SDwC_mean <<- SDwC_mean
      LFB_num_mode <<- LFB_num_mode

      cat("\n\nThe results obtained by REW-ISA are as follows:")
      cat("\noptimal_thr_row:", find_TR)
      cat("\noptimal_thr_col:", find_TC)
      cat("\noptimal_LFB_number:", LFB_number)
    }



    # Finding LFBs under optimized threshold
    if (optimization == FALSE) {
      if ((missing(optimal_thr_row) && missing(optimal_thr_col))) {
        stop("No optimal thresholds, please optimize the row and column thresholds first")
      }
      if ((!missing(optimal_thr_row) && !missing(optimal_thr_col))) {
        # Random seeds
        seeds <- generate.seeds(length = len_row, count = 100)
        isares <- isa.iterate(nm_w, row.seeds = seeds,
                              thr.row = optimal_thr_row, thr.col = optimal_thr_col)
        # Eliminate duplicated modules
        isares.unique <- isa.unique(nm_w, isares)
        # Filter out not robust ones
        isares2 <- isa.filter.robust(nm_w[[2]], nm_w, isares.unique)
        bc <- isa.biclust(isares2)
        LFB <<- bc
        return(bc)
      }
    }
  }
}
