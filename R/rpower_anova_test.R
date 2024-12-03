# Load necessary libraries
library(reshape2)
library(crayon)
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)

# Core power calculation function
power_anova_test <- function(n = NULL,
                             groups = NULL,
                             f = NULL,
                             power = NULL,
                             alpha = 0.05,
                             drop = 0,
                             type = c("between", "within", "mixed"),
                             repeated_measures = 1,
                             correlation = 0.5,
                             population = NULL) {
  
  type <- match.arg(type)
  
  # Validate inputs
  if (is.null(groups) || groups < 2) {
    stop("Number of groups must be specified and at least 2")
  }
  
  if (type != "between" && (repeated_measures < 2)) {
    stop("Number of repeated measures must be at least 2 for within or mixed designs")
  }
  
  # Calculate degrees of freedom based on type
  if (type == "between") {
    # Between-subjects ANOVA
    df1 <- groups - 1  # numerator df
    if (!is.null(n)) df2 <- groups * (n - 1)  # denominator df
  } else if (type == "within") {
    # Repeated measures ANOVA
    df1 <- repeated_measures - 1
    if (!is.null(n)) df2 <- (repeated_measures - 1) * (n - 1)
  } else {
    # Mixed type
    df1 <- (groups - 1) * (repeated_measures - 1)
    if (!is.null(n)) df2 <- groups * n * (repeated_measures - 1)
  }
  
  # Define power calculation based on type
  if (type == "between") {
    p = quote({
      df2 <- groups * (n - 1)
      lambda <- n * groups * f^2
      crit <- qf(1 - alpha, df1, df2)
      pf(crit, df1, df2, lambda, lower.tail = FALSE)
    })
  } else if (type == "within") {
    p = quote({
      df2 <- (repeated_measures - 1) * (n - 1)
      lambda <- n * repeated_measures * f^2 / (1 - correlation)
      crit <- qf(1 - alpha, df1, df2)
      pf(crit, df1, df2, lambda, lower.tail = FALSE)
    })
  } else {
    p = quote({
      df2 <- groups * n * (repeated_measures - 1)
      lambda <- n * groups * repeated_measures * f^2 / (1 - correlation)
      crit <- qf(1 - alpha, df1, df2)
      pf(crit, df1, df2, lambda, lower.tail = FALSE)
    })
  }
  
  # Solve based on which parameter is NULL
  if (is.null(power)) {
    res = eval(p)
  } else if (is.null(n)) {
    # Set bounds for sample size search
    lower_bound <- if(type == "between") 2 else 3
    upper_bound <- if(!is.null(population)) {
      floor(population/groups)
    } else {
      1e+06
    }
    
    tryCatch({
      res <- uniroot(
        function(n) {
          pow <- try(eval(p), silent = TRUE)
          if (inherits(pow, "try-error")) return(1)
          return(pow - power)
        },
        c(lower_bound, upper_bound),
        tol = 1e-6,
        extendInt = "yes"
      )$root
    }, error = function(e) {
      if (!is.null(population)) {
        stop("Could not achieve desired power with given population size constraints")
      } else {
        stop("Could not find valid sample size. Try adjusting effect size or power.")
      }
    })
  } else if (is.null(f)) {
    res <- uniroot(function(f) eval(p) - power, c(0.01, 2))$root
  } else if (is.null(alpha)) {
    res <- uniroot(function(alpha) eval(p) - power, c(1e-10, 1 - 1e-10))$root
  }
  
  return(res)
}

# Main function
rpower_anova_test <- function(n = NULL,
                              groups = NULL,
                              f = NULL,
                              power = NULL,
                              alpha = 0.05,
                              drop = 0,
                              population = NULL,
                              type = c("between", "within", "mixed"),
                              repeated_measures = 1,
                              correlation = 0.5,
                              plot = TRUE) {
  
  type <- match.arg(type)
  
  # Input validation
  if (is.null(groups) || !is.numeric(groups) || groups < 2)
    stop("Number of groups must be specified and at least 2")
  
  if (type != "between" && (repeated_measures < 2))
    stop("Number of repeated measures must be at least 2 for within or mixed designs")
  
  if (!is.null(population)) {
    if (!is.numeric(population))
      stop("Population size must be numeric")
    if (population <= 0)
      stop("Population size must be positive")
    if (population %% 1 != 0)
      stop("Population size must be a whole number")
    if (population < groups * 2)
      stop("Population size must be at least twice the number of groups")
  }
  
  if (!is.null(alpha) && (!is.numeric(alpha) || any(0 > alpha | alpha > 1)))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  
  if (!is.null(power) && (!is.numeric(power) || any(0 > power | power > 1)))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  
  if (!is.null(drop) && (!is.numeric(drop) || any(0 > drop | drop >= 1)))
    stop(sQuote("dropout rate"), " must be numeric in [0, 1)")
  
  if (type != "between") {
    if (!is.numeric(correlation) || correlation < -1 || correlation > 1)
      stop("Correlation must be between -1 and 1")
  }
  
  # Initialize lists for results
  p_list <- list()
  n_list <- list()
  f_list <- list()
  a_list <- list()
  
  # When n is NULL, solve for sample size
  if (is.null(n)) {
    for (i in seq_along(power)) {
      n_list2 <- list()
      for (j in seq_along(alpha)) {
        n_list3 <- list()
        for (k in seq_along(f)) {
          n_list3[[k]] <- power_anova_test(
            alpha = alpha[j],
            n = NULL,
            f = f[k],
            power = power[i],
            groups = groups,
            type = type,
            repeated_measures = repeated_measures,
            correlation = correlation,
            population = population
          )
        }
        names(n_list3) <- f
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
    
    # Calculate actual power with adjusted n
    res$Actual_Power <- mapply(function(n_adj, es, alpha) {
      as.numeric(power_anova_test(
        n = n_adj,
        f = es,
        power = NULL,
        alpha = alpha,
        groups = groups,
        type = type,
        repeated_measures = repeated_measures,
        correlation = correlation,
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
        for (k in seq_along(f)) {
          p_list3[[k]] <- power_anova_test(
            alpha = alpha[j],
            n = n[i],
            f = f[k],
            power = NULL,
            groups = groups,
            type = type,
            repeated_measures = repeated_measures,
            correlation = correlation,
            population = population
          )
        }
        names(p_list3) <- f
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
  
  # When f is NULL, calculate effect size
  else if (is.null(f)) {
    for (i in seq_along(n)) {
      f_list2 <- list()
      for (j in seq_along(power)) {
        f_list3 <- list()
        for (k in seq_along(alpha)) {
          f_list3[[k]] <- power_anova_test(
            n = n[i],
            f = NULL,
            power = power[j],
            alpha = alpha[k],
            groups = groups,
            type = type,
            repeated_measures = repeated_measures,
            correlation = correlation,
            population = population
          )
        }
        names(f_list3) <- alpha
        f_list2[[j]] <- f_list3
      }
      names(f_list2) <- power
      f_list[[i]] <- f_list2
    }
    names(f_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(f_list)), each = length(f_list[[1]][[1]])),
      "Target_Power" = rep(as.numeric(names(f_list[[1]])), times = length(f_list)),
      "Alpha" = as.numeric(names(f_list[[1]][[1]])),
      "Effect_Size" = as.numeric(unlist(f_list)),
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
        for (k in seq_along(f)) {
          a_list3[[k]] <- power_anova_test(
            n = n[i],
            f = f[k],
            power = power[j],
            alpha = NULL,
            groups = groups,
            type = type,
            repeated_measures = repeated_measures,
            correlation = correlation,
            population = population
          )
        }
        names(a_list3) <- f
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
  # Plotting code
  if (plot) {
    if (is.null(power) && !is.null(n)) {
      title = "Power vs Effect Size by Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = Effect_Size, y = Actual_Power, col = as.factor(Alpha))) +
        geom_point(aes(x = Effect_Size, y = Actual_Power, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Effect Size (f)",
             y = "Power",
             color = "Alpha") +
        theme_minimal()
    } else if (is.null(f)) {
      title = "Effect Size vs n by Power and Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = n, y = Effect_Size, col = as.factor(Alpha))) +
        geom_point(aes(x = n, y = Effect_Size, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Sample Size per Group (n)",
             y = "Effect Size (f)",
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
               x = "Sample Size per Group (n)",
               y = "Power",
               color = "Alpha") +
          facet_grid(~ Effect_Size) +
          theme_minimal()
      } else {
        p <- ggplot(data = res) +
          geom_line(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
          geom_point(aes(x = n, y = Actual_Power, col = as.factor(Alpha))) +
          labs(title = title,
               x = "Sample Size per Group (n)",
               y = "Power",
               color = "Alpha") +
          theme_minimal()
      }
    }
    
    # Add population size reference line if specified
    if (!is.null(population)) {
      p <- p + geom_vline(xintercept = population/groups,
                          linetype = "dashed",
                          color = "red",
                          alpha = 0.5) +
        annotate("text",
                 x = population/groups,
                 y = 0,
                 label = "Max Sample Size per Group",
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
    cat("Solve For ...........................", bold("Sample Size per Group"), "\n")
  } else if (!is.null(alpha) & !is.null(power) & is.null(f)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else if (is.null(alpha)) {
    cat("Solve For ...........................", bold("Alpha"), "\n")
  } else if (is.null(power) && !is.null(n)) {
    cat("Solve For ...........................", bold("Power"), "\n")
  }
  
  
  test_type <- switch(type,
                      "between" = "One-way Between-subjects ANOVA",
                      "within" = "One-way Repeated Measures ANOVA",
                      "mixed" = "Mixed Between-Within ANOVA")
  
  cat("Test ................................", bold(test_type), "\n")
  
  if (type != "between") {
    cat("Number of Measurements ..............", bold(repeated_measures), "\n")
    cat("Correlation .........................", bold(correlation), "\n")
  }
  
  cat("Number of Groups ....................", bold(groups), "\n")
  
  
  if (type == "within") {
    cat("Number of Measurements ..............", 
        bold(repeated_measures), "\n")
    cat("Correlation .........................", 
        bold(correlation), "\n")
  }
  
  cat("Number of Groups ....................", bold(groups), "\n")
  
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
    cat("n (Sample Size per Group) ...........", bold(n), "\n")
  }
  
  if (!is.null(f)) {
    cat("Effect size (f) .....................", bold(f), "\n")
  }
  
  if (!is.null(population)) {
    cat("Population Size .....................", bold(population), "\n")
  }
  
  cat("Degrees of Freedom ..................", 
      bold(paste0(
        if(type == "between") groups - 1 else repeated_measures - 1,
        ", ",
        if(type == "between" && !is.null(n)) {
          groups * (n - 1)
        } else if(type == "within" && !is.null(n)) {
          (repeated_measures - 1) * (n - 1)
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
        cat(bold("The required sample size per group, n' =", res$`n'`), "\n")
        cat(bold("This includes an adjustment for a", paste0(drop * 100, "%"), "dropout rate."), "\n")
        cat(bold("The sample size before dropout adjustment, n =", res$n), "\n")
        cat(bold("Total sample size =", res$`n'` * groups), "\n")
      } else {
        cat(bold("The required sample size per group, n =", res$n), "\n")
        cat(bold("Total sample size =", res$n * groups), "\n")
      }
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (!is.null(alpha) & !is.null(power) & is.null(f)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, f =", res$Effect_Size), "\n")
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