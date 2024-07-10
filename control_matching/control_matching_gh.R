rm(list = ls())

## Load the data
library(readxl)
library(tidyverse)

alldat <- list()

## RMFT fitness times
alldat$"RMFT" <-
  read_excel(
    "..."
  )
#View(alldat$"RMFT")

## PT Rank fitness times
alldat$"PT" <-
  read_excel(
    "..."
  )
#View(alldat$"PT")

## Bottom Field fitness times
alldat$"BFT" <-
  read_excel(
    "..."
  )
#View(alldat$"BFT")

## data loaded into list, now create z scores for BMI & fitness times by each specific troop.
calculate_z_scores <-
  function(data_list,
           vars_to_standardize = c("BMI", "Time")) {
    ## Loop through each data frame in the list
    for (dataset_name in names(data_list)) {
      dataset <- data_list[[dataset_name]]
      
      ## Loop through each variable to standardize (e.g., "BMI" and "Time")
      for (var_name in vars_to_standardize) {
        ## Calculate z-scores within each troop group
        z_score_col_name <- paste(var_name, "z", sep = "_")
        dataset[[z_score_col_name]] <-
          ave(
            dataset[[var_name]],
            dataset$Troop,
            FUN = function(x)
              scale(x, scale = TRUE)
          )
        
        ## If the variable is "Time," convert the z-scores to numeric
        if (var_name == "Time") {
          dataset[[z_score_col_name]] <-
            as.numeric(dataset[[z_score_col_name]])
        }
      }
      
      ## Update the data frame in the list
      data_list[[dataset_name]] <- dataset
    }
    
    return(data_list)
  }

## Now apply to data frames in list
alldat <- calculate_z_scores(alldat)

## Check data frames
View(alldat$"RMFT")
summary(alldat$"RMFT")
View(alldat$"PT")
summary(alldat$"PT")
View(alldat$"BFT")
summary(alldat$"BFT")

## Run across all data sets and find the matching control samples
for (ndat in names (alldat)) {
  ## Scan through each troop in turn
  dat <- alldat[[ndat]]
  all_troops <- unique(dat[["Troop"]])
  case_ctrl <- list()
  for (t in all_troops) {
    ## Get pairwise Euclidean distance between all samples
    dat_sub <- as.data.frame(subset(dat, Troop == t))
    rownames(dat_sub) <- dat_sub[["UIN"]]
    n <- nrow(dat_sub)
    dat_sub_diff <-
      matrix(
        NA,
        nrow = n,
        ncol = n,
        dimnames = list(rownames(dat_sub), rownames(dat_sub))
      )
    for (i in 1:n) {
      for (j in 1:n) {
        if (is.na(dat_sub[i, "Time_z"]) || is.na(dat_sub[j, "Time_z"])) {
          dat_sub_diff[i, j] <-
            sqrt((dat_sub[i, "BMI_z"] - dat_sub[j, "BMI_z"]) ^ 2)
        } else if (is.na(dat_sub[i, "BMI_z"]) ||
                   is.na(dat_sub[j, "BMI_z"])) {
          dat_sub_diff[i, j] <-
            sqrt((dat_sub[i, "Time_z"] - dat_sub[j, "Time_z"]) ^ 2)
        } else {
          dat_sub_diff[i, j] <-
            sqrt((dat_sub[i, "BMI_z"] - dat_sub[j, "BMI_z"]) ^ 2 + as.numeric(dat_sub[i, "Time_z"] -
                                                                                dat_sub[j, "Time_z"]) ^ 2)
        }
      }
    }
    
    ## For each case, find the best matching control
    ncase <- subset(dat_sub, Group %in% c(1, 2, 3))[["UIN"]]
    nctrl <- subset(dat_sub, Group %in% c(0))[["UIN"]]
    dat_sub_diff_ctrl <- dat_sub_diff[nctrl, ]
    case_ctrl_tmp <-
      matrix(
        NA,
        nrow = length(ncase),
        ncol = 2 + length(ncase),
        dimnames = list(ncase, c(
          "Troop", "Case", paste0("Control", 1:length(ncase))
        ))
      )
    case_ctrl_tmp <- data.frame(case_ctrl_tmp)
    case_ctrl_tmp[["Troop"]] <- t
    case_ctrl_tmp[["Case"]] <- ncase
    rownames(case_ctrl_tmp) <- ncase
    for (c in ncase) {
      case_diff <- dat_sub_diff_ctrl[, c]
      case_diff_sort <- sort(case_diff)
      case_ctrl_tmp[c, 3:ncol(case_ctrl_tmp)] <-
        names(case_diff_sort[1:length(ncase)])
    }
    
    ## If every sample has a unique match, use this
    case_ctrl[[t]] <- case_ctrl_tmp[, c(1:3)]
    names(case_ctrl[[t]])[3] <- "Control"
    
    ## If some samples have matching controls, select the optimum pairing
    if (length(unique(case_ctrl_tmp[["Control1"]])) != length(ncase)) {
      ## Look at the unique matched values
      all_ctrl <- unique(case_ctrl_tmp[["Control1"]])
      for (c in all_ctrl) {
        ## Get the Euclidean distances for all pairs
        ctrl_conflict <- subset(case_ctrl_tmp, Control1 == c)
        ctrl_conflict_sub <- ctrl_conflict[,-c(1:2)]
        ctrl_conflict_val <- ctrl_conflict_sub
        for (i in 1:nrow(ctrl_conflict_sub)) {
          for (j in 1:ncol(ctrl_conflict_sub)) {
            ctrl_conflict_val[i, j] <-
              dat_sub_diff[rownames(ctrl_conflict_sub)[i], ctrl_conflict_sub[i, j]]
          }
        }
        
        ## If there is a conflict, correct this
        if (nrow(ctrl_conflict) == 2) {
          ## Check all pairings and select the one with the least total distance
          pairings <- list()
          for (i in 1:ncol(ctrl_conflict_sub)) {
            pairings[[ctrl_conflict_sub[1, i]]] <- list()
            for (j in 1:ncol(ctrl_conflict_sub)) {
              if (i == j)
                next
              pairings[[ctrl_conflict_sub[1, i]]][[ctrl_conflict_sub[2, j]]] <-
                as.numeric(ctrl_conflict_val[1, i]) + as.numeric(ctrl_conflict_val[2, j])
            }
          }
          pairings <- unlist(pairings)
          
          ## Find the best match (making sure not to match to one of the other used controls)
          donotuse <- NULL
          for (nope in case_ctrl[[t]][case_ctrl[[t]][['Control']] != c, 'Control'])
            donotuse <- c(donotuse, grep(nope, names(pairings)))
          if (length(donotuse) != 0)
            pairings <- pairings[-donotuse]
          best_match <-
            strsplit(names(sort(pairings))[1], "\\.")[[1]]
          
          ## Update the list
          case_ctrl[[t]][rownames(ctrl_conflict_sub)[1], "Control"] <-
            best_match[1]
          case_ctrl[[t]][rownames(ctrl_conflict_sub)[2], "Control"] <-
            best_match[2]
          
        } else if (nrow(ctrl_conflict) > 2) {
          cat("Warning for Troop", t, ":\n")
          print(ctrl_conflict_sub)
          print(ctrl_conflict_val)
        }
      }
      
    }
    
  }
  
  ## Combine together
  case_ctrl <- do.call(rbind.data.frame, case_ctrl)
  write.csv(
    case_ctrl,
    paste0(
      "...",
      ndat,
      ".csv"
    )
  )
  
}
