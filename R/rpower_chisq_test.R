# Load necessary libraries
library(reshape2)
library(crayon)
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)

# Core power calculation function
power_chisq_test <- function(n = NULL,
                           w = NULL,  # Effect size (w)
                           power = NULL,
                           alpha = 0.05,
                           drop = 0,
                           type = c("gof", "independence"),  # goodness-of-fit or independence
                           df = NULL,  # degrees of freedom
                           population = NULL) {
  
  # Match arguments
  type <- match.arg(type)
  
  # Check population size only when n is not NULL
  if (!is.null(population) && !is.null(n)) {
    if (n > population) {
      stop("Sample size cannot be larger than population size")
    }
  }
  
  # If df is NULL, stop and request it
  if (is.null(df)) {
    stop("Degrees of freedom (df) must be specified")
  }
  
  # Define the power calculation function
  p = quote({
    ncp <- n * w^2  # Non-centrality parameter
    crit <- qchisq(1 - alpha, df)  # Critical value
    pchisq(crit, df, ncp, lower.tail = FALSE)  # Power calculation
  })
  
  # Solve based on which parameter is NULL
  if (is.null(power)) {
    res = eval(p)
  } else if (is.null(n)) {
    # For n calculation, ensure we don't exceed population size if specified
    upper_bound <- if (!is.null(population)) population else 1e+09
    lower_bound <- df + 1  # Minimum sample size based on df
    
    tryCatch({
      res <- uniroot(function(n) {
        pow <- try(eval(p), silent = TRUE)
        if (inherits(pow, "try-error")) return(1)
        return(pow - power)
      }, c(lower_bound, upper_bound))$root
    }, error = function(e) {
      if (!is.null(population)) {
        stop("Could not achieve desired power with given population size constraints")
      } else {
        stop("Could not find valid sample size. Try adjusting effect size or power.")
      }
    })
  } else if (is.null(w)) {
    res <- uniroot(function(w) eval(p) - power, c(0.01, 2))$root
  } else if (is.null(alpha)) {
    res <- uniroot(function(alpha) eval(p) - power, c(1e-10, 1 - 1e-10))$root
  }
  
  return(res)
}

# Main function
rpower_chisq_test <- function(n = NULL,
                            w = NULL,
                            power = NULL,
                            alpha = 0.05,
                            drop = 0,
                            population = NULL,
                            type = c("gof", "independence"),
                            df = NULL,
                            plot = TRUE) {
  
  # Enhanced error checking
  if (!is.null(population)) {
    if (!is.numeric(population))
      stop("Population size must be numeric")
    if (population <= 0)
      stop("Population size must be positive")
    if (population %% 1 != 0)
      stop("Population size must be a whole number")
    if (population < df + 1)
      stop(paste("Population size must be at least", df + 1, "for this analysis"))
  }
  
  if (!is.null(alpha) && (!is.numeric(alpha) || any(0 > alpha | alpha > 1)))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  
  if (!is.null(power) && (!is.numeric(power) || any(0 > power | power > 1)))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  
  if (!is.null(drop) && (!is.numeric(drop) || any(0 > drop | drop >= 1)))
    stop(sQuote("dropout rate"), " must be numeric in [0, 1)")
  
  if (is.null(df))
    stop("Degrees of freedom (df) must be specified")
  
  type <- match.arg(type)
  
  # Define lists for storing results
  p_list <- list()
  n_list <- list()
  w_list <- list()
  a_list <- list()
  
  # When n is NULL, solve for sample size
  # When n is NULL, solve for sample size
  if (is.null(n)) {
    for (i in seq_along(power)) {
      n_list2 <- list()
      for (j in seq_along(alpha)) {
        n_list3 <- list()
        for (k in seq_along(w)) {
          n_list3[[k]] <- power_chisq_test(
            alpha = alpha[j],
            n = NULL,
            w = w[k],
            power = power[i],
            type = type,
            df = df,
            population = population
          )
        }
        names(n_list3) <- w
        n_list2[[j]] <- n_list3
      }
      names(n_list2) <- alpha
      n_list[[i]] <- n_list2
    }
    names(n_list) <- power
    
    # Create result data frame
    res <- data.frame(
      "Target_Power" = rep(as.numeric(names(n_list)), each = length(n_list[[1]][[1]])),
      "Effect_Size" = rep(as.numeric(names(n_list[[1]][[1]])), times = length(n_list)),
      "Alpha" = as.numeric(names(n_list[[1]])),
      "n" = ceiling(as.numeric(unlist(n_list))),
      stringsAsFactors = FALSE
    )
    
    # Adjust n for dropout rate if drop > 0
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    }
    
    # Calculate actual power using the appropriate n value
    res$Actual_Power <- mapply(function(n_value, es, alpha) {
      n_to_use <- if(drop > 0) ceiling(n_value/(1-drop)) else n_value
      as.numeric(power_chisq_test(
        n = n_to_use,
        w = es,
        power = NULL,
        alpha = alpha,
        type = type,
        df = df,
        population = population
      ))
    }, n_value = res$n, es = res$Effect_Size, alpha = res$Alpha)
  }
  # When power is NULL, calculate power given n
  else if (is.null(power) && !is.null(n)) {
    p_list <- list()
    for (i in seq_along(n)) {
      p_list2 <- list()
      for (j in seq_along(alpha)) {
        p_list3 <- list()
        for (k in seq_along(w)) {
          p_list3[[k]] <- power_chisq_test(
            alpha = alpha[j],
            n = n[i],
            w = w[k],
            power = NULL,
            type = type,
            df = df,
            population = population
          )
        }
        names(p_list3) <- w
        p_list2[[j]] <- p_list3
      }
      names(p_list2) <- alpha
      p_list[[i]] <- p_list2
    }
    names(p_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(p_list)), each = length(p_list[[1]][[1]])),
      "Effect_Size" = rep(as.numeric(names(p_list[[1]][[1]])), times = length(p_list)),
      "Alpha" = as.numeric(names(p_list[[1]])),
      "Actual_Power" = as.numeric(unlist(p_list)),
      stringsAsFactors = FALSE
    )
    
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    } else {
      res$`n'` <- res$n
    }
  }
  
  # When w is NULL, calculate effect size
  else if (is.null(w)) {
    w_list <- list()
    for (i in seq_along(n)) {
      w_list2 <- list()
      for (j in seq_along(power)) {
        w_list3 <- list()
        for (k in seq_along(alpha)) {
          w_list3[[k]] <- power_chisq_test(
            n = n[i],
            w = NULL,
            power = power[j],
            alpha = alpha[k],
            type = type,
            df = df,
            population = population
          )
        }
        names(w_list3) <- alpha
        w_list2[[j]] <- w_list3
      }
      names(w_list2) <- power
      w_list[[i]] <- w_list2
    }
    names(w_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(w_list)), each = length(w_list[[1]][[1]])),
      "Target_Power" = rep(as.numeric(names(w_list[[1]])), times = length(w_list)),
      "Alpha" = as.numeric(names(w_list[[1]][[1]])),
      "Effect_Size" = as.numeric(unlist(w_list)),
      stringsAsFactors = FALSE
    )
    
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    } else {
      res$`n'` <- res$n
    }
    
    res$Actual_Power <- res$Target_Power
  }
  
  # When alpha is NULL, calculate alpha given n, w, and power
  else if (is.null(alpha)) {
    # Code to solve for alpha and create res
    a_list <- list()
    for (i in seq_along(n)) {
      a_list2 <- list()
      for (j in seq_along(power)) {
        a_list3 <- list()
        for (k in seq_along(w)) {
          a_list3[[k]] <- power_chisq_test(
            n = n[i],
            w = w[k],
            power = power[j],
            alpha = NULL,
            type = type,
            df = df,
            population = population
          )
        }
        names(a_list3) <- w
        a_list2[[j]] <- a_list3
      }
      names(a_list2) <- power
      a_list[[i]] <- a_list2
    }
    names(a_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(a_list)), each = length(a_list[[1]][[1]])),
      "Target_Power" = rep(as.numeric(names(a_list[[1]])), times = length(a_list)),
      "Effect_Size" = rep(as.numeric(names(a_list[[1]][[1]])), times = length(a_list) * length(a_list[[1]])),
      "Alpha" = as.numeric(unlist(a_list)),
      stringsAsFactors = FALSE
    )
    
    # Adjust n for dropout rate if drop > 0
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    } else {
      res$`n'` <- res$n
    }
    
    # Calculate actual power 
    res$Actual_Power <- res$Target_Power
  }
  
  # Plotting code
  if (plot) {
    if (is.null(power) && !is.null(n)) {
      title = "Power vs Effect Size by Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = `Effect_Size`, y = `Actual_Power`, col = as.factor(Alpha))) +
        geom_point(aes(x = `Effect_Size`, y = `Actual_Power`, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Effect Size (w)",
             y = "Power",
             color = "Alpha") +
        theme_minimal()
    } else if (is.null(w)) {
      title = "Effect Size vs n by Power and Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = n, y = `Effect_Size`, col = as.factor(Alpha))) +
        geom_point(aes(x = n, y = `Effect_Size`, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Sample Size (n)",
             y = "Effect Size (w)",
             color = "Alpha") +
        facet_grid(~ `Target_Power`) +
        theme_minimal()
    } else {
      title = ifelse(length(unique(res$`Effect_Size`)) > 1,
                     "Power vs n by Alpha and Effect Size",
                     "Power vs n by Alpha")
      
      if (length(unique(res$`Effect_Size`)) > 1) {
        res$`Effect_Size` = round(as.numeric(as.character(res$`Effect_Size`)), 3)
        
        p <- ggplot(data = res) +
          geom_line(aes(x = n, y = `Actual_Power`, col = as.factor(Alpha))) +
          geom_point(aes(x = n, y = `Actual_Power`, col = as.factor(Alpha))) +
          labs(title = title,
               x = "Sample Size (n)",
               y = "Power",
               color = "Alpha") +
          facet_grid(~ `Effect_Size`) +
          theme_minimal()
      } else {
        p <- ggplot(data = res) +
          geom_line(aes(x = n, y = `Actual_Power`, col = as.factor(Alpha))) +
          geom_point(aes(x = n, y = `Actual_Power`, col = as.factor(Alpha))) +
          labs(title = title,
               x = "Sample Size (n)",
               y = "Power",
               color = "Alpha") +
          theme_minimal()
      }
    }
    
    # Add population size reference line if specified
    if (!is.null(population)) {
      p <- p + geom_vline(xintercept = population,
                          linetype = "dashed",
                          color = "red",
                          alpha = 0.5) +
        annotate("text",
                 x = population,
                 y = 0,
                 label = "Population Size",
                 vjust = -0.5,
                 angle = 90,
                 color = "red",
                 alpha = 0.5)
    }
    
    p <- p + theme(legend.position = "bottom") +
      scale_y_continuous(limits = c(0, NA))
    
    print(p)
  }
  
  
  # Printing the details
  if (is.null(n)) {
    cat("Solve For ...........................", bold("Sample Size"), "\n")
  } else if (!is.null(alpha) & !is.null(power) & is.null(w)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else if (is.null(alpha)) {
    cat("Solve For ...........................", bold("Alpha"), "\n")
  } else if (is.null(power) && !is.null(n)) {
    cat("Solve For ...........................", bold("Power"), "\n")
  } else if (is.null(w)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else {
    cat("Solve For ...........................", bold("Power"), "\n")
  }
  
  # Print test type
  cat("Test ................................", 
      bold(ifelse(type == "gof", "Chi-square Goodness-of-Fit Test", 
                  "Chi-square Test of Independence")), "\n")
  
  # Print all parameters
  if (drop > 0) {
    cat("Dropout Rate ........................", bold(paste0(drop * 100, "%")), "\n")
  }
  
  if (!is.null(power)) {
    cat("Power ...............................", bold(power), "\n")
  }
  
  if (!is.null(alpha)) {
    cat("Alpha ...............................", bold(alpha), "\n")
  }
  
  if (!is.null(n)) {
    cat("n (Sample Size) .....................", bold(n), "\n")
  }
  
  if (!is.null(w)) {
    cat("Effect size (w) .....................", bold(w), "\n")
  }
  
  if (!is.null(population)) {
    cat("Population Size .....................", bold(population), "\n")
  }
  
  cat("Degrees of Freedom ..................", bold(df), "\n")
  
  cat(".............................................................\n")
  cat("\n")
  
  # Final result messages
  if (!is.null(power) & is.null(n)) {
    if (nrow(res) == 1) {
      if (drop > 0) {
        cat(bold("The required total sample size, n' =", res$`n'`), "\n")
        cat(bold("This includes an adjustment for a", paste0(drop * 100, "%"), "dropout rate."), "\n")
        cat(bold("The sample size before dropout adjustment, n =", res$n), "\n")
      } else {
        cat(bold("The required total sample size, n =", res$n), "\n")
      }
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (!is.null(alpha) & !is.null(power) & is.null(w)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, w =", res$`Effect_Size`), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (is.null(alpha)) {
    if (nrow(res) == 1) {
      cat(bold("The alpha level, alpha =", res$Alpha), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (is.null(w)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, w =", res$`Effect_Size`), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (is.null(power) && !is.null(n)) {
    if (nrow(res) == 1) {
      cat(bold("The actual power of the test, p =", res$`Actual_Power`), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  }
  
  # Remove columns with all NA values and n' when drop = 0 before returning
  res <- res[, colSums(is.na(res)) < nrow(res)]
  if (drop == 0 && "n'" %in% names(res)) {
    res <- res[, !names(res) %in% "n'"]
  }
  
  return(res)
  
  # Remove columns with all NA values before returning
  res <- res[, colSums(is.na(res)) < nrow(res)]
  
  return(res)
}