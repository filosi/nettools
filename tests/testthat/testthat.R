context("Check nettools consistency.")


set.seed(0)
a <- matrix(rnorm(1000),ncol=100)
set.seed(1)
b <- matrix(rnorm(1000),ncol=100)

test_that("Test MINE:", {
    aadj <- mat2adj(a, infer.method="MINE", measure="MIC", n.cores=1)
    dima <- dim(aadj)
    expect_equal(dima[1], 100)
    expect_equal(dima[2], 100)    
})

test_that("Test WGCNA:",{
    aadj <- mat2adj(a, infer.method="WGCNA", P=6, n.cores=1, ciccio=3)
    dima <- dim(aadj)
    expect_equal(dima[1], 100)
    expect_equal(dima[2], 100)
})

test_that("Test bicorFDR:",{
    skip("Do not test multicore...")
    aadj <- mat2adj(a, infer.method="bicorFDR", P=6, FDR=1e-2, n.cores=1, ciccio=3)
    dima <- dim(aadj)
    expect_equal(dima[1], 100)
    expect_equal(dima[2], 100)
})

test_that("Test c3net:", {
    aadj <- mat2adj(a, infer.method="c3net", n.cores=1, ciccio=2)
    dima <- dim(aadj)
    expect_equal(dima[1], 100)
    expect_equal(dima[2], 100)    
})

test_that("Test empty-full distance:", {
    e <- matrix(0, ncol=100, nrow=100)
    f <- e + 1
    diag(f) <- 0
    dd <- netdist(e, f, d="HIM", components=FALSE, n.cores=1)
    expect_equal(as.numeric(dd), 1)
})


test_that("Test distance:", {
    aadj <- mat2adj(a, infer.method="WGCNA", n.cores=1)
    badj <- mat2adj(b, infer.method="WGCNA", n.cores=1)
    d1 <- netdist(aadj,badj, d="HIM", components=FALSE, n.cores=1)
    d2 <- netdist(list(aadj,badj), d="HIM", components=FALSE, n.cores=1)
})

test_that("Test SI:", {
    ssind1 <- netSI(a, adj.method="ARACNE", n.cores=1, use="pair", symm=TRUE)
    ssind2 <- netSI(a, n.cores=1, adj.method="ARACNE", use="pair", symm=TRUE)
    expect_equal(length(ssind1$Sw), length(ssind2$Sw))
})


