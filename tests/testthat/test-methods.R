test_that("Testing getCorM function", {
    library(Matrix)
    rand_vals <- sample(0:10,1e4,replace=T,p=c(0.99,rep(0.001,10)))
    x <- as.matrix(Matrix(rand_vals,ncol=5))
    w=replicate(nrow(x),1)

    expect_equal(getCorM("spearman", x, nCores=4), cor(x, method="spearman"))
    expect_equal(getCorM("pearson", x, nCores=4), cor(x, method="pearson"))
    expect_equal(getCorM("spearman", x, w, 4), cor(x, method="spearman"))
    expect_equal(getCorM("pearson", x, w, 4), cor(x, method="pearson"))
})

test_that("Testing colRanks_fast function", {
    library(Matrix)
    rand_vals <- sample(0:10,1e4,replace=T,p=c(0.99,rep(0.001,10)))
    x <- as.matrix(Matrix(rand_vals,ncol=5))
    expect_equal(colRanks_fast(x, 4), apply(x, 2, rank))
})

