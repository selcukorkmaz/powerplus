# Example 1: Basic two-sample t-test with single effect size
rpower_t_test(
  d = 0.5,                    # Medium effect size
  power = NULL,               # Standard power
  alpha = 0.05, 
  n = 200, # Standard alpha
  type = "two.sample",        
  alternative = "two.sided",
  drop = 0.10,
  plot = TRUE
)

# Example 2: Multiple effect sizes and power levels
rpower_t_test(
  d = c(0.2, 0.5, 0.8),      # Small, medium, large effect sizes
  power = c(0.80, 0.90),      # Two power levels
  alpha = 0.05,
  type = "two.sample",
  alternative = "two.sided",
  plot = TRUE
)

# Example 3: Different alpha levels
rpower_t_test(
  d = 0.5,                    
  power = 0.80,               
  alpha = c(0.01, 0.05, 0.10), # Three alpha levels
  type = "two.sample",
  alternative = "two.sided",
  plot = TRUE
)

# Example 4: One-sample t-test
rpower_t_test(
  d = 0.5,
  power = NULL,
  alpha = 0.05,
  n = 100,
  type = "one.sample",        # Changed to one-sample
  alternative = "two.sided",
  plot = TRUE
)

# Example 5: Different alternatives
example5 <- rpower_t_test(
  d = 0.5,
  power = 0.80,
  alpha = 0.05,
  type = "two.sample",
  alternative = "greater",    # One-sided test
  plot = TRUE
)

# Example 6: With population constraint
example6 <- rpower_t_test(
  d = 0.5,
  power = 0.80,
  alpha = 0.05,
  type = "two.sample",
  alternative = "two.sided",
  population = 1000,          # Fixed population size
  plot = TRUE
)

# Example 7: With dropout rate
example7 <- rpower_t_test(
  d = 0.5,
  power = 0.80,
  alpha = 0.05,
  drop = 0.2,                # 20% dropout rate
  type = "two.sample",
  alternative = "two.sided",
  plot = TRUE
)

# Example 8: Paired t-test
example8 <- rpower_t_test(
  d = 0.5,
  power = 0.80,
  alpha = 0.05,
  type = "paired",           # Paired test
  alternative = "two.sided",
  plot = TRUE
)

# Example 9: Unequal group sizes
example9 <- rpower_t_test(
  d = 0.5,
  power = 0.80,
  alpha = 0.05,
  type = "two.sample",
  alternative = "two.sided",
  ratio = 2,                 # Ratio of n2/n1 = 2
  plot = TRUE
)

# Example 10: Complex combination with multiple parameters
example10 <- rpower_t_test(
  d = c(0.3, 0.5, 0.7),
  power = c(0.80, 0.85, 0.90),
  alpha = c(0.01, 0.05),
  type = "two.sample",
  alternative = "two.sided",
  plot = TRUE
)