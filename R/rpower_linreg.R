# Load necessary libraries
library(reshape2)
library(crayon)
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)

# Modified core power calculation function with corrected calculations
power_linreg <- function(n = NULL,
                           predictors = NULL,
                           f2 = NULL,
                           power = NULL,
                           alpha = 0.05,
                           drop = 0,
                           population = NULL) {
  
  if (is.null(predictors) || predictors < 1) {
    stop("Number of predictors must be specified and at least 1")
  }
  
  if (!is.null(population) && !is.null(n)) {
    if (n > population) {
      stop("Sample size cannot be larger than population size")
    }
  }
  
  # Degrees of freedom
  u <- predictors  # numerator df
  if (!is.null(n)) v <- n - predictors - 1  # denominator df
  
  # Modified power calculation with correct F distribution parameters
  p = quote({
    v <- n - predictors - 1
    # Calculate non-centrality parameter
    lambda <- n * f2
    # Critical F value
    f_crit <- qf(1 - alpha, u, v)
    # Power calculation
    pf(f_crit, u, v, lambda, lower.tail = FALSE)
  })
  
  if (is.null(power)) {
    res = eval(p)
  } else if (is.null(n)) {
    upper_bound <- if (!is.null(population)) population else 1e+06
    lower_bound <- predictors + 2
    
    tryCatch({
      res <- uniroot(
        function(n) {
          pow <- try(eval(p), silent = TRUE)
          if (inherits(pow, "try-error")) return(1)
          return(pow - power)
        },
        c(lower_bound, upper_bound),
        tol = 1e-10,  # Further increased precision
        maxiter = 1000  # Increased maximum iterations
      )$root
    }, error = function(e) {
      if (!is.null(population)) {
        stop("Could not achieve desired power with given population size constraints")
      } else {
        stop("Could not find valid sample size. Try adjusting effect size or power.")
      }
    })
  } else if (is.null(f2)) {
    res <- uniroot(function(f2) eval(p) - power, c(0.001, 2))$root
  } else if (is.null(alpha)) {
    res <- uniroot(function(alpha) eval(p) - power, c(1e-10, 1 - 1e-10))$root
  }
  
  return(res)
}

# Main function
rpower_linreg <- function(n = NULL,
                            predictors = NULL,
                            f2 = NULL,
                            power = NULL,
                            alpha = 0.05,
                            drop = 0,
                            population = NULL,
                            plot = TRUE) {
  
  # Enhanced error checking
  if (is.null(predictors) || !is.numeric(predictors) || predictors < 1)
    stop("Number of predictors must be specified and at least 1")
  
  if (!is.null(population)) {
    if (!is.numeric(population))
      stop("Population size must be numeric")
    if (population <= 0)
      stop("Population size must be positive")
    if (population %% 1 != 0)
      stop("Population size must be a whole number")
    if (population < predictors + 2)
      stop(paste("Population size must be at least", predictors + 2, "for this analysis"))
  }
  
  if (!is.null(alpha) && (!is.numeric(alpha) || any(0 > alpha | alpha > 1)))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  
  if (!is.null(power) && (!is.numeric(power) || any(0 > power | power > 1)))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  
  if (!is.null(drop) && (!is.numeric(drop) || any(0 > drop | drop >= 1)))
    stop(sQuote("dropout rate"), " must be numeric in [0, 1)")
  
  # Initialize lists for results
  p_list <- list()
  n_list <- list()
  f2_list <- list()
  a_list <- list()
  
  # When n is NULL, solve for sample size
  if (is.null(n)) {
    for (i in seq_along(power)) {
      n_list2 <- list()
      for (j in seq_along(alpha)) {
        n_list3 <- list()
        for (k in seq_along(f2)) {
          n_list3[[k]] <- power_linreg(
            alpha = alpha[j],
            n = NULL,
            f2 = f2[k],
            power = power[i],
            predictors = predictors,
            population = population
          )
        }
        names(n_list3) <- f2
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
    
    # Calculate actual power
    res$Actual_Power <- mapply(function(n_adj, es, alpha) {
      as.numeric(power_linreg(
        n = n_adj,
        f2 = es,
        power = NULL,
        alpha = alpha,
        predictors = predictors,
        population = population
      ))
    }, n_adj = if(drop > 0) res$`n'` else res$n,
    es = res$Effect_Size,
    alpha = res$Alpha)
  }
  
  # When power is NULL, calculate power given n
  else if (is.null(power)) {
    for (i in seq_along(n)) {
      p_list2 <- list()
      for (j in seq_along(alpha)) {
        p_list3 <- list()
        for (k in seq_along(f2)) {
          p_list3[[k]] <- power_linreg(
            alpha = alpha[j],
            n = n[i],
            f2 = f2[k],
            power = NULL,
            predictors = predictors,
            population = population
          )
        }
        names(p_list3) <- f2
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
    }
  }
  
  # When f2 is NULL, calculate effect size
  else if (is.null(f2)) {
    for (i in seq_along(n)) {
      f2_list2 <- list()
      for (j in seq_along(power)) {
        f2_list3 <- list()
        for (k in seq_along(alpha)) {
          f2_list3[[k]] <- power_linreg(
            n = n[i],
            f2 = NULL,
            power = power[j],
            alpha = alpha[k],
            predictors = predictors,
            population = population
          )
        }
        names(f2_list3) <- alpha
        f2_list2[[j]] <- f2_list3
      }
      names(f2_list2) <- power
      f2_list[[i]] <- f2_list2
    }
    names(f2_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(f2_list)), each = length(f2_list[[1]][[1]])),
      "Target_Power" = rep(as.numeric(names(f2_list[[1]])), times = length(f2_list)),
      "Alpha" = as.numeric(names(f2_list[[1]][[1]])),
      "Effect_Size" = as.numeric(unlist(f2_list)),
      stringsAsFactors = FALSE
    )
    
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    }
    
    res$Actual_Power <- res$Target_Power
  }
  
  # When alpha is NULL, calculate alpha
  else if (is.null(alpha)) {
    for (i in seq_along(n)) {
      a_list2 <- list()
      for (j in seq_along(power)) {
        a_list3 <- list()
        for (k in seq_along(f2)) {
          a_list3[[k]] <- power_linreg(
            n = n[i],
            f2 = f2[k],
            power = power[j],
            alpha = NULL,
            predictors = predictors,
            population = population
          )
        }
        names(a_list3) <- f2
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
    
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    }
    
    res$Actual_Power <- res$Target_Power
  }
  
  # Plotting code
  if (plot) {
    if (is.null(power) && !is.null(n)) {
      title = "Power vs Effect Size by Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = Effect_Size, y = Actual_Power, col = as.factor(Alpha))) +
        geom_point(aes(x = Effect_Size, y = Actual_Power, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Effect Size (f²)",
             y = "Power",
             color = "Alpha") +
        theme_minimal()
    } else if (is.null(f2)) {
      title = "Effect Size vs n by Power and Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = n, y = Effect_Size, col = as.factor(Alpha))) +
        geom_point(aes(x = n, y = Effect_Size, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Sample Size (n)",
             y = "Effect Size (f²)",
             color = "Alpha") +
        facet_grid(~ Target_Power) +
        theme_minimal()
    } else {
      title = ifelse(length(unique(res$Effect_Size)) > 1,
                     "Power vs n by Alpha and Effect Size",
                     "Power vs n by Alpha")
      
      if (length(unique(res$Effect_Size)) > 1) {
        res$Effect_Size = round(as.numeric(as.character(res$Effect_Size)), 3)
        
        p <- ggplot(data = res) +
          geom_line(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
          geom_point(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
          labs(title = title,
               x = "Sample Size (n)",
               y = "Power",
               color = "Alpha") +
          facet_grid(~ Effect_Size) +
          theme_minimal()
      } else {
        p <- ggplot(data = res) +
          geom_line(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
          geom_point(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
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
  } else if (!is.null(alpha) & !is.null(power) & is.null(f2)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else if (is.null(alpha)) {
    cat("Solve For ...........................", bold("Alpha"), "\n")
  } else if (is.null(power) && !is.null(n)) {
    cat("Solve For ...........................", bold("Power"), "\n")
  }
  
  cat("Test ................................", 
      bold("Multiple Linear Regression"), "\n")
  
  cat("Number of Predictors ................", bold(predictors), "\n")
  
  if (drop > 0) {
    cat("Dropout Rate ........................", 
        bold(paste0(drop * 100, "%")), "\n")
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
  
  if (!is.null(f2)) {
    cat("Effect size (f²) ....................", bold(f2), "\n")
  }
  
  if (!is.null(population)) {
    cat("Population Size .....................", bold(population), "\n")
  }
  
  cat("Degrees of Freedom ..................", 
      bold(paste0(
        predictors,
        ", ",
        if(!is.null(n)) {
          n - predictors - 1
        } else {
          "calculated based on final n"
        }
      )), "\n")
  
  cat(".............................................................\n")
  cat("\n")
  
  # Final result messages
  if (!is.null(power) & is.null(n)) {
    if (nrow(res) == 1) {
      if (drop > 0) {
        cat(bold("The required sample size, n' =", res$`n'`), "\n")
        cat(bold("This includes an adjustment for a", paste0(drop * 100, "%"), "dropout rate."), "\n")
        cat(bold("The sample size before dropout adjustment, n =", res$n), "\n")
      } else {
        cat(bold("The required sample size, n =", res$n), "\n")
      }
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (!is.null(alpha) & !is.null(power) & is.null(f2)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, f² =", res$Effect_Size), "\n")
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
  } else if (is.null(power) && !is.null(n)) {
    if (nrow(res) == 1) {
      cat(bold("The actual power of the test, p =", res$Actual_Power), "\n")
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
}