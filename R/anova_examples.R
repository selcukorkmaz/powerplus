# One-way between-subjects ANOVA examples

# Example 1: Calculate required sample size per group
# 3 groups, medium effect size (f = 0.25), 80% power
rpower_anova_test(
  n = NULL,
  groups = 3,
  f = 0.25,
  power = 0.80,
  alpha = 0.05,
  type = "between"
)

# Example 2: Calculate power for given sample size
# 4 groups, n = 30 per group, large effect size (f = 0.40)
rpower_anova_test(
  n = 30,
  groups = 4,
  f = 0.40,
  power = NULL,
  alpha = 0.05,
  type = "between"
)

# Example 3: Calculate minimum detectable effect size
# 3 groups, n = 25 per group, 80% power
rpower_anova_test(
  n = 25,
  groups = 3,
  f = NULL,
  power = 0.80,
  alpha = 0.05,
  type = "between"
)

# Example 4: Between-subjects ANOVA with dropout rate
rpower_anova_test(
  n = NULL,
  groups = 3,
  f = 0.25,
  power = c(0.80,0.9),
  alpha = 0.05,
  drop = 0,
  type = "between"
)

# One-way repeated measures ANOVA examples

# Example 5: Calculate required sample size for repeated measures
# 4 measurements, medium effect size, moderate correlation
rpower_anova_test(
  n = NULL,
  groups = 2,
  f = 0.25,
  power = 0.80,
  alpha = 0.05,
  type = "within",
  repeated_measures = 4,
  correlation = 0.5
)

# Example 6: Power analysis for repeated measures with high correlation
rpower_anova_test(
  n = 20,
  groups = 2,
  f = 0.30,
  power = NULL,
  alpha = 0.05,
  type = "within",
  repeated_measures = 3,
  correlation = 0.8
)

# Example 7: Multiple effect sizes and alpha levels
rpower_anova_test(
  n = NULL,
  groups = 3,
  f = c(0.10, 0.25, 0.40),  # Small, medium, large effects
  power = 0.80,
  alpha = c(0.01, 0.05),
  type = "between"
)

# Example 8: Population size constraint
rpower_anova_test(
  n = NULL,
  groups = 2,
  f = 0.25,
  power = 0.80,
  alpha = 0.05,
  population = 200,
  type = "between"
)

# Example 9: Effect size calculation with repeated measures
rpower_anova_test(
  n = 30,
  groups = 3,
  f = NULL,
  power = 0.90,
  alpha = 0.05,
  type = "within",
  repeated_measures = 4,
  correlation = 0.6
)

# Example 10: Alpha level calculation
rpower_anova_test(
  n = 40,
  groups = 3,
  f = 0.25,
  power = 0.85,
  alpha = NULL,
  type = "between"
)
