power_corr_test <- function(
                         n,
                         r,
                         power,
                         alpha,
                         drop,
                         alternative
                         ) {
  
  if (alternative == "two.sided" && !is.null(r)){ 
    r <- abs(r)
    
  }
  
  if (alternative == "less") {
    p <- quote({
      ttt <- qt(alpha, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3))
    })
  }
  
  if (alternative == "two.sided") {
    p = quote({
      r <- abs(-r)
      ttt <- qt(alpha, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3))
    })
  }
  
  if (alternative == "greater") {
    p = quote({
      ttt <- qt(alpha/2, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3)) + pnorm((-zr - zrc) * 
                                                sqrt(n - 3))
      
    })
  }
  
  if (is.null(power)) {
    res = eval(p)
  }
  
  else if (is.null(n)) {
    res <- uniroot(function(n) eval(p) - power, c(4 + 1e-10, 1e+09))$root
  }
  
  
  else if (is.null(r)) {
    if (alternative == "two.sided") 
      res <- uniroot(function(r) eval(p) - power, c(1e-10, 
                                                       1 - 1e-10))$root
    else res <- uniroot(function(r) eval(p) - power, c(-1 + 
                                                            1e-10, 1 - 1e-10))$root
  }
  
  else if (is.null(alpha)){ 
    res <- uniroot(function(alpha) eval(p) - 
                     power, c(1e-10, 1 - 1e-10))$root
  }
  return(res)
  
}

rpower_corr_test <-
  function(
    n = NULL,
    r = NULL,
    power = NULL,
    alpha = NULL,
    drop = 0,
    alternative = c("two.sided", "less", "greater"),
    plot = TRUE) {
    # if (sum(sapply(list(n, power, r), is.null)) != 1)
    #   stop("exactly one of n, power and r must be NULL")
    
    if (!is.null(alpha) &&
        !is.numeric(alpha) || any(0 > alpha | alpha > 1))
      stop(sQuote("alpha"), " must be numeric in [0, 1]")
    
    if (!is.null(power) &&
        !is.numeric(power) || any(0 > power | power > 1))
      stop(sQuote("power"), " must be numeric in [0, 1]")
    
    if (!is.null(drop) &&
        !is.numeric(drop) || any(0 > drop | drop >= 1))
      stop(sQuote("dropout rate"), " must be numeric in [0, 1)")
    

    alternative <- match.arg(alternative)
    
    p_list = list()
    n_list = list()
    r_list <- list()
    a_list <- list()
    
    if (is.null(n)) {
      
      for (i in 1:length(power)) {
        n_list2 <- list()
        
        for (j in 1:length(alpha)) {
          n_list3 <- list()
          
          for (k in 1:length(r)) {
            n_list3[[k]] <- power_corr_test(
              alpha = alpha[j],
              n = NULL,
              r = r[k],
              power = power[i],
              alternative = alternative
            )
          }
          names(n_list3) = r 
          
          
          n_list2[[j]] = n_list3
        }
        
        names(n_list2) = alpha 
        
        n_list[[i]] <- n_list2
        
        
        
      }
      names(n_list) = power
      
      for (i in 1:length(ceiling(unlist(n_list)))) {
        p_list2 <- list()
        
        for (j in 1:length(alpha)) {
          p_list3 <- list()
          
          for (k in 1:length(r)) {
            p_list3[[k]] <- power_corr_test(
              alpha = alpha[j],
              n = ceiling(unlist(n_list))[i],
              r = r[k],
              power = NULL,
              alternative = alternative
            )[[1]]
          }
          names(p_list3) = r 
          
          
          p_list2[[j]] = p_list3
        }
        
        names(p_list2) = alpha 
        
        p_list[[i]] <- p_list2
        
        
        
      }
      names(p_list) = ceiling(unlist(n_list))
      
      
      
    }
    
    if (is.null(power)){
      for (i in 1:length(n)) {
        p_list2 <- list()
        
        for (j in 1:length(alpha)) {
          p_list3 <- list()
          
          for (k in 1:length(r)) {
            p_list3[[k]] <- power_corr_test(
              alpha = alpha[j],
              n = n[i],
              r = r[k],
              power = NULL,
              alternative = alternative
            )
          }
          names(p_list3) = r 
          
          
          p_list2[[j]] = p_list3
        }
        
        names(p_list2) = alpha 
        
        p_list[[i]] <- p_list2
        
        
        
      }
      names(p_list) = n
    }
    
    if (is.null(r)){
      for (i in 1:length(n)) {
        d_list2 <- list()
        
        for (j in 1:length(power)) {
          d_list3 <- list()
          
          for (k in 1:length(alpha)) {
            d_list3[[k]] <- power_corr_test(
              alpha = alpha[k],
              n = n[i],
              r = NULL,
              power = power[j],
              alternative = alternative
            )
          }
          names(d_list3) = alpha 
          
          
          d_list2[[j]] = d_list3
        }
        
        names(d_list2) = power 
        
        r_list[[i]] <- d_list2
        
      }
      names(r_list) = n
    }
    
    if (is.null(alpha)){
      for (i in 1:length(n)) {
        a_list2 <- list()
        
        for (j in 1:length(power)) {
          a_list3 <- list()
          
          for (k in 1:length(r)) {
            a_list3[[k]] <- power_corr_test(
              r = r[k],
              n = n[i],
              alpha = NULL,
              power = power[j],
              alternative = alternative
            )
          }
          names(a_list3) = r 
          
          
          a_list2[[j]] = a_list3
        }
        
        names(a_list2) = power 
        
        a_list[[i]] <- a_list2
        
      }
      names(a_list) = n
    }
    
    if(!is.null(n_list) && is.null(n)){
      
      
      res <- data.frame(unlist(n_list))
      
      res$target_power <- str_sub( rownames(res), start = 1,
                                   end = rownames(res) %>%
                                     str_locate_all(pattern = '\\.') %>%
                                     map_int(extract(2))%>%
                                     subtract(1))
      
      
      res$es <- str_sub(rownames(res), start = rownames(res) %>%
                          str_locate_all(pattern = '\\.') %>%
                          map_int(extract(4)) %>%
                          add(1),
                        end = nchar(rownames(res)))
      
      
      res$alpha <- str_sub(rownames(res), start = rownames(res) %>%
                             str_locate_all(pattern = '\\.') %>%
                             map_int(extract(2)) %>%
                             add(1),
                           end = rownames(res) %>%
                             str_locate_all(pattern = '\\.') %>%
                             map_int(extract(4))%>%
                             subtract(1))
      
      
      
      names_p_list = gsub(" ", ".", str_c(ceiling(res$unlist.n_list.), res$alpha, res$es, sep = " "))
      
      res$actual_power <- unlist(p_list)[names(unlist(p_list)) %in% names_p_list]
      
      rownames(res) = NULL
      
      res <- res[,c(2,5,1,3,4)]
      res$drop = NA
      res$n_drop = NA
      res$D = NA
      
      colnames(res) = c(
        "Target Power",
        "Actual Power",
        "n",
        "Effect Size",
        "Alpha",
        "Dropout Rate",
        "n'",
        "D"
      )
    }
    
    if(!is.null(p_list) && is.null(power)){
      
      res <- data.frame(unlist(p_list))
      res$n <- sub("\\..*", "", rownames(res))
      
      
      res$es <- str_sub(rownames(res), start = rownames(res) %>%
                          str_locate_all(pattern = '\\.') %>%
                          map_int(extract(3)) %>%
                          add(1),
                        end = nchar(rownames(res)))
      
      res$alpha <- str_sub( rownames(res), start = rownames(res) %>%
                              str_locate_all(pattern = '\\.') %>%
                              map_int(extract(1)) %>%
                              add(1),
                            end = rownames(res) %>%
                              str_locate_all(pattern = '\\.') %>%
                              map_int(extract(3))%>%
                              subtract(1))
      
      
      rownames(res) = NULL
      
      colnames(res) = c("Actual Power", "n", "Effect Size", "Alpha")  
      
    }
    
    if(!is.null(r_list) && is.null(r)){
      
      res <- data.frame(unlist(r_list))
      res$n <- sub("\\..*", "", rownames(res))
      res$target_power <- str_sub( rownames(res), start = rownames(res) %>%
                                     str_locate_all(pattern = '\\.') %>%
                                     map_int(extract(1)) %>%
                                     add(1),
                                   end = rownames(res) %>%
                                     str_locate_all(pattern = '\\.') %>%
                                     map_int(extract(3))%>%
                                     subtract(1))
      
      res$alpha <- str_sub(rownames(res), start = rownames(res) %>%
                             str_locate_all(pattern = '\\.') %>%
                             map_int(extract(3)) %>%
                             add(1),
                           end = nchar(rownames(res)))
      
      
      rownames(res) = NULL
      
      res <- res[,c(3,2,1,4)]
      
      colnames(res) = c("Target Power", "n", "Effect Size", "Alpha")  
      
    }
    
    if(!is.null(a_list) && is.null(alpha)){
      
      res <- data.frame(unlist(a_list))
      res$n <- sub("\\..*", "", rownames(res))
      res$target_power <- str_sub( rownames(res), start = rownames(res) %>%
                                     str_locate_all(pattern = '\\.') %>%
                                     map_int(extract(1)) %>%
                                     add(1),
                                   end = rownames(res) %>%
                                     str_locate_all(pattern = '\\.') %>%
                                     map_int(extract(3))%>%
                                     subtract(1))
      
      res$es <- str_sub(rownames(res), start = rownames(res) %>%
                          str_locate_all(pattern = '\\.') %>%
                          map_int(extract(3)) %>%
                          add(1),
                        end = nchar(rownames(res)))
      
      
      rownames(res) = NULL
      
      res <- res[,c(3,2,4,1)]
      
      colnames(res) = c("Target Power", "n", "Effect Size", "Alpha")  
      
    }

    if (is.null(n)) {
      
      if (drop != 0){ 
        res[6] = paste0(drop * 100, "%")
        res[7] = ceiling(res$n / (1 - drop))
        res[8] = res[, 7] - res$n
      }
    }
    
    if (plot) {
      side = ifelse(alternative == "two.sided", "two-sided", "one-sided")
      title = ifelse(length(unique(res$`Effect Size`)) > 1,
                     "Power vs n by Alpha and Effect Size",
                     "Power vs n by Alpha")
      subt = ifelse(
        is.null(r),
        paste0(
          "\U03BC",
          "0=",
          m1,
          " ",
          "\U03BC",
          "1=",
          m2,
          " ",
          "\U03C3",
          "1=",
          sp,
          " ",
          side,
          " one-sample t-test"
        ),
        paste0("Based on effect size: " , side, " one-sample t-test")
      )
      
      if (length(unique(res$`Effect Size`)) > 1) {
        res$`Effect Size` = round(res$`Effect Size`, 2)
        
        p <-
          ggplot(data = res) + geom_line(aes(x = n, y = `Actual Power`, col = Alpha)) +
          geom_point(aes(x = n, y = `Actual Power`, col = Alpha)) +
          labs(title = title,
               subtitle = subt) +
          facet_grid( ~ `Effect Size`)
        
      } else{
        p <-
          ggplot(data = res) + geom_line(aes(x = n, y = `Actual Power`, col = Alpha)) +
          geom_point(aes(x = n, y = `Actual Power`, col = Alpha)) +
          labs(title = title,
               subtitle = subt)
        
      }
      print(p)
    }
    
    
    # cat(bold("Design\n"))
    # cat(".............................................................\n")
    if (is.null(n)) {
      cat("Solve For ...........................",
          bold("Sample Size"),
          "\n")
    } 
    
    else if(!is.null(alpha) &
            !is.null(power) & is.null(r)){
      cat("Solve For ...........................",
          bold("Effect Size"),
          "\n")
    }
    
    else if (is.null(alpha)) {
      cat("Solve For ...........................",
          bold("Alpha"),
          "\n")
    } else{
      cat("Solve For ...........................",
          bold("Power"),
          "\n")
    }
    
    cat("Test ................................",
        bold("Pearson's Correlation"),
        "\n")
    
    
    if (is.null(n)) {
      cat("Dropout Rate ........................",
          bold(paste0(drop * 100, "%")),
          "\n")
    }
    
    if (!is.null(power)) {
      cat("Power ...............................",
          bold(power),
          "\n")
    }
    
    if (!is.null(alpha)) {
      cat("Alpha ...............................", bold(alpha), "\n")
    }
    if (!is.null(n)) {
      cat("N (Sample Size) .....................", bold(n), "\n")
    }
    
    if (!is.null(r)) {
      cat("Effect size .........................", bold(r), "\n")
      
    }
    
    
    
    cat(".............................................................\n")
    cat("\n")
    
    
    
    if (!is.null(power) & is.null(n)){
      
      if(nrow(res)==1){  
        cat(bold("The required total sample size, n =", sample), "\n")
        cat("\n")
        
      }else{
        
        cat(bold("Please see table for results."), "\n")
        cat("\n")
        
      }
      
    }
    
    else if (!is.null(alpha) &
             !is.null(power) & is.null(r)){
      
      if(nrow(res) == 1){
        cat(bold("The effect size, r =", res$`Effect Size`), "\n")
        cat("\n")
      }else{
        
        cat(bold("Please see table for results."), "\n")
        cat("\n")
        
      }
      
    }
    
    else if (is.null(alpha)){
      
      if(nrow(res) == 1){
        cat(bold("The alpha level, alpha =", res$Alpha), "\n")
        cat("\n")
      }else{
        
        cat(bold("Please see table for results."), "\n")
        cat("\n")
        
      }
      
    }else{
      
      if(nrow(res) == 1){
        
        cat(bold("The actual power of the test, p =", res$`Actual Power`), "\n")
        cat("\n")
        
      }else{
        
        cat(bold("Please see table for results."), "\n")
        cat("\n")
        
      }
    }
    
    
    res <- res[, colSums(is.na(res)) < nrow(res)]
    return(res)
    
    
    
  }

rpower_corr_test(
  n = c(100,200),
  r = NULL,
  power = c(0.8),
  alpha = c(0.05),
  drop = 0.1,
  alternative = "two.sided",
  plot = F
)
