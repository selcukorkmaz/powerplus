
## t testi
rpower_t_test(
  d = 0.5,                   
  power = 0.8,              
  alpha = 0.05, 
  n = NULL,
  type = "two.sample",        
  alternative = "two.sided",
  drop = 0.10,
  plot = TRUE
)

## ki-kare testi
rpower_chisq_test(
  df = 4,
  n = NULL,
  w = 0.5,
  power = 0.8,
  alpha = 0.05,
  drop = 0,
  plot = F
)

## korelasyon testi
rpower_corr_test(
  n = 100,
  r = 0.5,
  power = 0.8,
  alpha = NULL,
  drop = 0,
  alternative = "two.sided",
  plot = F
)