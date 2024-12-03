es_t <- function(n1 = NULL, n2 = NULL, m1 = NULL, m2 = NULL, s1 = NULL, s2 = NULL, r = NULL,
                 type = "two.samples", method = "Cohen"){
  
  if(type == "two.samples"){
    if(method == "Cohen"){
      md = abs((m1-m2)) 
      sd = sqrt((s1^2 + s2^2)/2)
      d = md/sd
    }
    
    if(method == "Glass"){
      
      md = abs((m1-m2)) 
      sd = s1
      d = md/sd
    }
    
    if(method == "Hedge"){
      sp =  sqrt(((s1^2)*(n1-1) + (s2^2)*(n2-1))/(n1+n2-2))
      md =  abs((m1-m2))
      sd = sp
      d = md/sd
    }
  }
  
  if(type == "one.sample"){
    md = abs((m1-m2)) 
    sd = s1
    d = md/sd
  }
  
  if(type == "paired"){
    if(is.null(m2) & is.null(s2)){
      md = abs((m1)) 
      sd = s1
      d = md/sd
    }else{
      md = abs((m1-m2))
      sd = sqrt((s1^2+s2^2-2*r))
      d = md/sd

    }
  }
  
  cat("Solve For ...........................",
      bold("Effect Size"),
      "\n")
  
  if (type == "one.sample") {
    cat("Test ................................",
        bold("One-Sample t-Test"),
        "\n")
  }
  
  if (type == "two.samples") {
    cat("Test ................................",
        bold("Two-Samples t-Test"),
        "\n")
  }
  
  if (type == "paired") {
    cat("Test ................................",
        bold("Paired-Samples t-Test"),
        "\n")
  }

  if (type == "two.samples") {
    if(method == "Cohen"){
    cat("Method ..............................",
        bold("Cohen's d"),
        "\n")
    }
    
    if(method == "Glass"){
      cat("Method ..............................",
          bold("Glass's g"),
          "\n")
    }
    
    if(method == "Hedge"){
      cat("Method ..............................",
          bold("Hedge's h"),
          "\n")
    }
  }
  
  if (type == "one.samples" | type == "paired") {
    if(method == "Cohen"){
      cat("Method ..............................",
          bold("Cohen's d"),
          "\n")
    }
  }
  
  
    cat("\n")
  

  
  cat(bold("The effect size, d =", round(d,4)), "\n")
  res <- cbind.data.frame(md, sd, d)
  colnames(res) <- c("Mean Diff.", "Standard Dev.", "Effect Size")
  cat("\n")
   return(res)
}

n1= c(70,90)
n2=50
m1=c(73.64)
m2=c(61.87)
s1=c(9.547)
s2=6.625
r=-0.207

es_t(n1, n2, m1, m2, s1, s2, r, type="two.samples")
