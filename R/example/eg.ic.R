# Generate a trivial dataset, X has mean 0 and norm 1, y has mean 0
set.seed(11)
n = 20
p = 5
x = matrix(rnorm(n*p), nrow=n, ncol=p)
x = scale(x, center = colMeans(x))
x = scale(x, scale = sqrt(colSums(x^2)))
beta = c(1, 1, 0, 0, 0)
y = x%*%beta + scale(rnorm(20, sd=0.01), center = TRUE, scale = FALSE)

# run the model
boss_result = boss(x, y)
# what are the values of AICc for all candidates of BOSS?
print(boss_result$IC_boss$aicc)
# calculate them manually using the calc.ic function
print(calc.ic(x%*%boss_result$beta_boss, y, 'aicc', boss_result$hdf_boss))
# they shall match
