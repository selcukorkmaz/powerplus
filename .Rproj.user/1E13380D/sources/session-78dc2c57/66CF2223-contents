# Load necessary libraries
library(reshape2)
library(crayon)
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
library(dplyr)

# Core power calculation function
power_t_test <- function(n = NULL,
                         d = NULL,
                         power = NULL,
                         alpha = 0.05,
                         drop = 0,
                         type = c("two.sample", "one.sample", "paired"),
                         alternative = c("two.sided", "less", "greater"),
                         ratio = 1,
                         sd1 = 1,
                         sd2 = NULL,
                         population = NULL) {
  
  # Set sd2 to sd1 if not specified
  if (is.null(sd2)) sd2 <- sd1
  
  # Match arguments
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  
  # Check population size only when n is not NULL
  if (!is.null(population) && !is.null(n)) {
    if (n > population) {
      stop("Sample size cannot be larger than population size")
    }
  }
  
  hyptohesis <- switch(type,
                       "one.sample" = 1,
                       "two.sample" = 2,
                       "paired" = 1)
  
  if (alternative == "less") {
    p = quote({
      if (type == "two.sample") {
        n1 <- n / (1 + ratio)
        n2 <- n - n1
        v1 <- sd1^2 / n1
        v2 <- sd2^2 / n2
        df <- (v1 + v2)^2 / ((v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)))
        ncp <- d / sqrt(v1 + v2)
      } else {
        df <- (n - 1) * hyptohesis
        ncp <- sqrt(n / hyptohesis) * d
      }
      pt(
        qt(alpha, df, lower = TRUE),
        df,
        ncp = ncp,
        lower = TRUE
      )
    })
  }
  
  if (alternative == "two.sided") {
    p = quote({
      if (type == "two.sample") {
        n1 <- n / (1 + ratio)
        n2 <- n - n1
        v1 <- sd1^2 / n1
        v2 <- sd2^2 / n2
        df <- (v1 + v2)^2 / ((v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)))
        ncp <- d / sqrt(v1 + v2)
      } else {
        df <- (n - 1) * hyptohesis
        ncp <- sqrt(n / hyptohesis) * d
      }
      ta <- qt(alpha / 2, df, lower = FALSE)
      pt(ta, df, ncp = ncp, lower = FALSE) +
        pt(-ta, df, ncp = ncp, lower = TRUE)
    })
  }
  
  if (alternative == "greater") {
    p = quote({
      if (type == "two.sample") {
        n1 <- n / (1 + ratio)
        n2 <- n - n1
        v1 <- sd1^2 / n1
        v2 <- sd2^2 / n2
        df <- (v1 + v2)^2 / ((v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)))
        ncp <- d / sqrt(v1 + v2)
      } else {
        df <- (n - 1) * hyptohesis
        ncp <- sqrt(n / hyptohesis) * d
      }
      pt(
        qt(alpha, df, lower = FALSE),
        df,
        ncp = ncp,
        lower = FALSE
      )
    })
  }
  
  if (is.null(power)) {
    res = eval(p)
  } else if (is.null(n)) {
    # For n calculation, ensure we don't exceed population size if specified
    upper_bound <- if (!is.null(population)) population else 1e+09
    lower_bound <- if (type == "two.sample") 4 else 2
    
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
  } else if (is.null(d)) {
    if (alternative == "two.sided")
      res <- uniroot(function(d) eval(p) - power, c(1e-07, 10))$root
    if (alternative == "less")
      res <- uniroot(function(d) eval(p) - power, c(-10, 5))$root
    if (alternative == "greater")
      res <- uniroot(function(d) eval(p) - power, c(-5, 10))$root
  } else if (is.null(alpha)) {
    res <- uniroot(function(alpha) eval(p) - power, c(1e-10, 1 - 1e-10))$root
  }
  
  return(res)
}

# Main function
rpower_t_test <- function(n = NULL,
                          d = NULL,
                          power = NULL,
                          alpha = 0.05,
                          drop = 0,
                          population = NULL,
                          type = c("two.sample", "one.sample", "paired"),
                          alternative = c("two.sided", "less", "greater"),
                          ratio = 1,
                          plot = TRUE,
                          sd1 = 1,
                          sd2 = NULL) {
  # Enhanced error checking for population parameter
  if (!is.null(population)) {
    if (!is.numeric(population))
      stop("Population size must be numeric")
    if (population <= 0)
      stop("Population size must be positive")
    if (population %% 1 != 0)
      stop("Population size must be a whole number")
    if (type == "two.sample" && population < 4)
      stop("Population size must be at least 4 for two-sample tests")
    if (type != "two.sample" && population < 2)
      stop("Population size must be at least 2")
  }
  
  # Set sd2 to sd1 if not specified
  if (is.null(sd2)) sd2 <- sd1
  
  if (!is.null(alpha) && (!is.numeric(alpha) || any(0 > alpha | alpha > 1)))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  
  if (!is.null(power) && (!is.numeric(power) || any(0 > power | power > 1)))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  
  if (!is.null(drop) && (!is.numeric(drop) || any(0 > drop | drop >= 1)))
    stop(sQuote("dropout rate"), " must be numeric in [0, 1)")
  
  if (type == "two.sample") {
    if (!is.null(ratio) && (!is.numeric(ratio) || any(0.1 > ratio | ratio > 100)))
      stop(sQuote("ratio"), " must be numeric in [0.1, 100]")
  }
  
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  
  # Define lists for storing results
  p_list <- list()
  n_list <- list()
  d_list <- list()
  a_list <- list()
  
  # When n is NULL, solve for sample size
  if (is.null(n)) {
    # Existing code to solve for n and create res
    for (i in seq_along(power)) {
      n_list2 <- list()
      for (j in seq_along(alpha)) {
        n_list3 <- list()
        for (k in seq_along(d)) {
          n_list3[[k]] <- power_t_test(
            alpha = alpha[j],
            n = NULL,
            d = d[k],
            power = power[i],
            type = type,
            alternative = alternative,
            ratio = ratio,
            sd1 = sd1,
            sd2 = sd2,
            population = population
          )
        }
        names(n_list3) <- d
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
    } else {
      res$`n'` <- res$n
    }
    
    # Adjust n1 and n2 if it's a two-sample test
    if (type == "two.sample") {
      if (ratio == 1) {
        # Adjust n' to be even where necessary
        odd_indices <- which(res$`n'` %% 2 != 0)
        res$`n'`[odd_indices] <- res$`n'`[odd_indices] + 1
        
        res$n1 <- res$`n'` / 2
        res$n2 <- res$`n'` / 2
      } else if (ratio >= 1) {
        r2 <- ratio
        r1 <- 1
        res$n1 <- round(res$`n'` / (r1 + r2) * r1)
        res$n2 <- res$`n'` - res$n1
      } else {
        r2 <- 1
        r1 <- 1 / ratio
        res$n2 <- round(res$`n'` / (r1 + r2) * r2)
        res$n1 <- res$`n'` - res$n2
      }
    }
    
    # Recalculate actual power with adjusted n'
    res$Actual_Power <- mapply(function(n_adj, es, alpha) {
      as.numeric(power_t_test(
        n = n_adj,
        d = es,
        power = NULL,
        alpha = alpha,
        type = type,
        alternative = alternative,
        ratio = ratio,
        sd1 = sd1,
        sd2 = sd2,
        population = population
      ))
    }, n_adj = res$`n'`, es = res$Effect_Size, alpha = res$Alpha)
    
    # Add remaining columns
    res$D <- NA_real_
    
    # Reorder columns
    if (drop > 0) {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "n'", "Effect_Size", "Alpha", "n1", "n2", "Dropout_Rate", "D")
    } else {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "Effect_Size", "Alpha", "n1", "n2", "D")
    }
    
    for (col in expected_columns) {
      if (!(col %in% names(res))) {
        res[[col]] <- NA_real_
      }
    }
    res <- res[, expected_columns]
    
    # Remove columns with all NA values
    res <- res[, colSums(is.na(res)) < nrow(res)]
  }
  # When power is NULL, calculate power given n
  else if (is.null(power) && !is.null(n)) {
    # Code block to calculate power and create res
    p_list <- list()
    for (i in seq_along(n)) {
      p_list2 <- list()
      for (j in seq_along(alpha)) {
        p_list3 <- list()
        for (k in seq_along(d)) {
          p_list3[[k]] <- power_t_test(
            alpha = alpha[j],
            n = n[i],
            d = d[k],
            power = NULL,
            type = type,
            alternative = alternative,
            ratio = ratio,
            sd1 = sd1,
            sd2 = sd2,
            population = population
          )
        }
        names(p_list3) <- d
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
    
    # Adjust n for dropout rate if drop > 0
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    } else {
      res$`n'` <- res$n
    }
    
    # Adjust n1 and n2 if it's a two-sample test
    if (type == "two.sample") {
      if (ratio == 1) {
        # Adjust n' to be even where necessary
        odd_indices <- which(res$`n'` %% 2 != 0)
        res$`n'`[odd_indices] <- res$`n'`[odd_indices] + 1
        
        res$n1 <- res$`n'` / 2
        res$n2 <- res$`n'` / 2
      } else if (ratio >= 1) {
        r2 <- ratio
        r1 <- 1
        res$n1 <- round(res$`n'` / (r1 + r2) * r1)
        res$n2 <- res$`n'` - res$n1
      } else {
        r2 <- 1
        r1 <- 1 / ratio
        res$n2 <- round(res$`n'` / (r1 + r2) * r2)
        res$n1 <- res$`n'` - res$n2
      }
    }
    
    # Recalculate actual power with adjusted n'
    res$Actual_Power <- mapply(function(n_adj, es, alpha) {
      as.numeric(power_t_test(
        n = n_adj,
        d = es,
        power = NULL,
        alpha = alpha,
        type = type,
        alternative = alternative,
        ratio = ratio,
        sd1 = sd1,
        sd2 = sd2,
        population = population
      ))
    }, n_adj = res$`n'`, es = res$Effect_Size, alpha = res$Alpha)
    
    # Add remaining columns
    res$Target_Power <- NA_real_
    res$D <- NA_real_
    
    # Reorder columns
    if (drop > 0) {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "n'", "Effect_Size", "Alpha", "n1", "n2", "Dropout_Rate", "D")
    } else {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "Effect_Size", "Alpha", "n1", "n2", "D")
    }
    
    for (col in expected_columns) {
      if (!(col %in% names(res))) {
        res[[col]] <- NA_real_
      }
    }
    res <- res[, expected_columns]
    
    # Remove columns with all NA values
    res <- res[, colSums(is.na(res)) < nrow(res)]
  }
  # When d is NULL, calculate effect size given n and power
  else if (is.null(d)) {
    # Code to solve for d and create res
    d_list <- list()
    for (i in seq_along(n)) {
      d_list2 <- list()
      for (j in seq_along(power)) {
        d_list3 <- list()
        for (k in seq_along(alpha)) {
          d_list3[[k]] <- power_t_test(
            n = n[i],
            d = NULL,
            power = power[j],
            alpha = alpha[k],
            type = type,
            alternative = alternative,
            ratio = ratio,
            sd1 = sd1,
            sd2 = sd2,
            population = population
          )
        }
        names(d_list3) <- alpha
        d_list2[[j]] <- d_list3
      }
      names(d_list2) <- power
      d_list[[i]] <- d_list2
    }
    names(d_list) <- n
    
    # Create result data frame
    res <- data.frame(
      "n" = rep(as.numeric(names(d_list)), each = length(d_list[[1]][[1]])),
      "Target_Power" = rep(as.numeric(names(d_list[[1]])), times = length(d_list)),
      "Alpha" = as.numeric(names(d_list[[1]][[1]])),
      "Effect_Size" = as.numeric(unlist(d_list)),
      stringsAsFactors = FALSE
    )
    
    # Adjust n for dropout rate if drop > 0
    if (drop > 0) {
      res$Dropout_Rate <- drop
      res$`n'` <- ceiling(res$n / (1 - drop))
    } else {
      res$`n'` <- res$n
    }
    
    # Adjust n1 and n2 if it's a two-sample test
    if (type == "two.sample") {
      if (ratio == 1) {
        # Adjust n' to be even where necessary
        odd_indices <- which(res$`n'` %% 2 != 0)
        res$`n'`[odd_indices] <- res$`n'`[odd_indices] + 1
        
        res$n1 <- res$`n'` / 2
        res$n2 <- res$`n'` / 2
      } else if (ratio >= 1) {
        r2 <- ratio
        r1 <- 1
        res$n1 <- round(res$`n'` / (r1 + r2) * r1)
        res$n2 <- res$`n'` - res$n1
      } else {
        r2 <- 1
        r1 <- 1 / ratio
        res$n2 <- round(res$`n'` / (r1 + r2) * r2)
        res$n1 <- res$`n'` - res$n2
      }
    }
    
    # Recalculate actual power with adjusted n'
    res$Actual_Power <- res$Target_Power
    res$D <- NA_real_
    
    # Reorder columns
    if (drop > 0) {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "n'", "Effect_Size", "Alpha", "n1", "n2", "Dropout_Rate", "D")
    } else {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "Effect_Size", "Alpha", "n1", "n2", "D")
    }
    
    for (col in expected_columns) {
      if (!(col %in% names(res))) {
        res[[col]] <- NA_real_
      }
    }
    res <- res[, expected_columns]
    
    # Remove columns with all NA values
    res <- res[, colSums(is.na(res)) < nrow(res)]
  }
  # When alpha is NULL, calculate alpha given n, d, and power
  else if (is.null(alpha)) {
    # Code to solve for alpha and create res
    a_list <- list()
    for (i in seq_along(n)) {
      a_list2 <- list()
      for (j in seq_along(power)) {
        a_list3 <- list()
        for (k in seq_along(d)) {
          a_list3[[k]] <- power_t_test(
            n = n[i],
            d = d[k],
            power = power[j],
            alpha = NULL,
            type = type,
            alternative = alternative,
            ratio = ratio,
            sd1 = sd1,
            sd2 = sd2,
            population = population
          )
        }
        names(a_list3) <- d
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
    
    # Adjust n1 and n2 if it's a two-sample test
    if (type == "two.sample") {
      if (ratio == 1) {
        # Adjust n' to be even where necessary
        odd_indices <- which(res$`n'` %% 2 != 0)
        res$`n'`[odd_indices] <- res$`n'`[odd_indices] + 1
        
        res$n1 <- res$`n'` / 2
        res$n2 <- res$`n'` / 2
      } else if (ratio >= 1) {
        r2 <- ratio
        r1 <- 1
        res$n1 <- round(res$`n'` / (r1 + r2) * r1)
        res$n2 <- res$`n'` - res$n1
      } else {
        r2 <- 1
        r1 <- 1 / ratio
        res$n2 <- round(res$`n'` / (r1 + r2) * r2)
        res$n1 <- res$`n'` - res$n2
      }
    }
    
    # Recalculate actual power with adjusted n'
    res$Actual_Power <- res$Target_Power
    res$D <- NA_real_
    
    # Reorder columns
    if (drop > 0) {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "n'", "Effect_Size", "Alpha", "n1", "n2", "Dropout_Rate", "D")
    } else {
      expected_columns <- c("Target_Power", "Actual_Power", "n", "Effect_Size", "Alpha", "n1", "n2", "D")
    }
    
    for (col in expected_columns) {
      if (!(col %in% names(res))) {
        res[[col]] <- NA_real_
      }
    }
    res <- res[, expected_columns]
    
    # Remove columns with all NA values
    res <- res[, colSums(is.na(res)) < nrow(res)]
  }
  # Handle other scenarios as per existing code
  else {
    stop("Invalid input: please check the parameters provided.")
  }
  
  # Plotting code
  if (plot) {
    side = ifelse(alternative == "two.sided", "Two-sided", "One-sided")
    
    if (is.null(power) && !is.null(n)) {
      title = "Power vs Effect Size by Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = `Effect_Size`, y = `Actual_Power`, col = as.factor(Alpha))) +
        geom_point(aes(x = `Effect_Size`, y = `Actual_Power`, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Effect Size",
             y = "Power",
             color = "Alpha") +
        theme_minimal()
    } else if (is.null(d)) {
      title = "Effect Size vs n by Power and Alpha"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = n, y = `Effect_Size`, col = as.factor(Alpha))) +
        geom_point(aes(x = n, y = `Effect_Size`, col = as.factor(Alpha))) +
        labs(title = title,
             x = "Sample Size (n)",
             y = "Effect Size",
             color = "Alpha") +
        facet_grid(~ `Target_Power`) +
        theme_minimal()
    } else if (is.null(alpha)) {
      title = "Alpha vs n by Power and Effect Size"
      
      p <- ggplot(data = res) +
        geom_line(aes(x = n, y = Alpha, col = as.factor(`Effect_Size`))) +
        geom_point(aes(x = n, y = Alpha, col = as.factor(`Effect_Size`))) +
        labs(title = title,
             x = "Sample Size (n)",
             y = "Alpha",
             color = "Effect Size") +
        facet_grid(~ `Target_Power`) +
        theme_minimal()
    } else {
      title = ifelse(length(unique(res$`Effect_Size`)) > 1,
                     "Power vs n by Alpha and Effect Size",
                     "Power vs n by Alpha")
      
      if (length(unique(res$`Effect_Size`)) > 1) {
        res$`Effect_Size` = as.numeric(as.character(res$`Effect_Size`))
        res$`Effect_Size` = round(res$`Effect_Size`, 2)
        
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
  } else if (!is.null(alpha) & !is.null(power) & is.null(d)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else if (is.null(alpha)) {
    cat("Solve For ...........................", bold("Alpha"), "\n")
  } else if (is.null(power) && !is.null(n)) {
    cat("Solve For ...........................", bold("Power"), "\n")
  } else if (is.null(d)) {
    cat("Solve For ...........................", bold("Effect Size"), "\n")
  } else {
    cat("Solve For ...........................", bold("Power"), "\n")
  }
  
  if (type == "one.sample") {
    cat("Test ................................", bold("One-Sample t-Test"), "\n")
  }
  
  if (type == "two.sample") {
    cat("Test ................................", bold("Two-Samples t-Test"), "\n")
  }
  
  if (type == "paired") {
    cat("Test ................................", bold("Paired-Samples t-Test"), "\n")
  }
  
  if (alternative == "two.sided") {
    if (type == "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("Two-Sided (H1: μ ≠ μ₀)"), "\n")
    }
    if (type != "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("Two-Sided (H1: μ₁ ≠ μ₂)"), "\n")
    }
  }
  
  if (alternative == "less") {
    if (type == "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("One-Sided (H1: μ < μ₀)"), "\n")
    }
    if (type != "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("One-Sided (H1: μ₁ < μ₂)"), "\n")
    }
  }
  
  if (alternative == "greater") {
    if (type == "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("One-Sided (H1: μ > μ₀)"), "\n")
    }
    if (type != "one.sample") {
      cat("Alternative Hypothesis ..............",
          bold("One-Sided (H1: μ₁ > μ₂)"), "\n")
    }
  }
  
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
  
  if (!is.null(d)) {
    cat("Effect size .........................", bold(d), "\n")
  }
  
  if (!is.null(population)) {
    cat("Population Size .....................", bold(population), "\n")
  }
  
  if (type == "two.sample") {
    cat("SD group 1 ..........................", bold(sd1), "\n")
    if (!is.null(sd2)) {
      cat("SD group 2 ..........................", bold(sd2), "\n")
    }
  }
  
  if (type == "two.sample" & is.null(n)) {
    cat("Allocation ratio (n2/n1) ............", bold(ratio), "\n")
  }
  
  cat(".............................................................\n")
  cat("\n")
  
  # Final messages
  if (!is.null(power) & is.null(n)) {
    if (nrow(res) == 1) {
      if (drop > 0) {
        cat(bold("The required total sample size, n' =", res$`n'`), "\n")
        if (type == "two.sample") {
          cat(bold("This includes an adjustment for a", paste0(drop * 100, "%"), "dropout rate."), "\n")
          cat(bold("Group sizes: n1 =", res$n1, "and n2 =", res$n2), "\n")
        } else {
          cat(bold("This includes an adjustment for a", paste0(drop * 100, "%"), "dropout rate."), "\n")
        }
      } else {
        cat(bold("The required total sample size, n =", res$n), "\n")
        if (type == "two.sample") {
          cat(bold("Group sizes: n1 =", res$n1, "and n2 =", res$n2), "\n")
        }
      }
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else if (!is.null(alpha) & !is.null(power) & is.null(d)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, d =", res$`Effect_Size`), "\n")
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
  } else if (is.null(d)) {
    if (nrow(res) == 1) {
      cat(bold("The effect size, d =", res$`Effect_Size`), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  } else {
    if (nrow(res) == 1) {
      cat(bold("The actual power of the test, p =", res$`Actual_Power`), "\n")
      cat("\n")
    } else {
      cat(bold("Please see table for results."), "\n")
      cat("\n")
    }
  }
  
  # Remove columns with all NA values before returning
  res <- res[, colSums(is.na(res)) < nrow(res)]
  
  return(res)
}







