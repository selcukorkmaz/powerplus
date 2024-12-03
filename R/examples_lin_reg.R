# Example 1: Calculate required sample size
# 3 predictors, medium effect size (f² = 0.15), 80% power
rpower_linreg(
  n = NULL,
  predictors = 3,
  f2 = 0.15,
  power = 0.80,
  alpha = 0.05
)

# Example 2: Calculate power for given sample size
# 4 predictors, n = 100, medium effect size
rpower_linreg(
  n = 100,
  predictors = 4,
  f2 = 0.15,
  power = NULL,
  alpha = 0.05
)

# Example 3: Calculate minimum detectable effect size
# 2 predictors, n = 80, 80% power
rpower_linreg(
  n = 80,
  predictors = 2,
  f2 = NULL,
  power = 0.80,
  alpha = 0.05
)

# Example 4: With dropout rate
rpower_linreg(
  n = NULL,
  predictors = 3,
  f2 = 0.15,
  power = 0.80,
  alpha = 0.05,
  drop = 0.1
)

# Example 5: Multiple effect sizes and alphas
rpower_linreg(
  n = NULL,
  predictors = 3,
  f2 = c(0.02, 0.15, 0.35),  # Small, medium, large
  power = 0.80,
  alpha = c(0.01, 0.05)
)

# Example 6: Population constraint
rpower_linreg(
  n = NULL,
  predictors = 2,
  f2 = 0.15,
  power = 0.80,
  alpha = 0.05,
  population = 200
)

# Example 7: Calculate alpha level
rpower_linreg(
  n = 100,
  predictors = 3,
  f2 = 0.15,
  power = 0.80,
  alpha = NULL
)

# Example 8: Power analysis with many predictors
rpower_linreg(
  n = 200,
  predictors = 8,
  f2 = 0.15,
  power = NULL,
  alpha = 0.05
)

# Example 9: Sample size for high power
rpower_linreg(
  n = NULL,
  predictors = 4,
  f2 = 0.15,
  power = 0.95,
  alpha = 0.05
)

