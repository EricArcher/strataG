context("alleleSplit validation")

test_that("missing loci preserve shape and labels", {
  x <- matrix(c(NA_character_,NA_character_),2,1,dimnames=list(c("a","b"),"L"))
  expect_identical(alleleSplit(x,sep="/"),matrix(NA_character_,2,2,
    dimnames=list(c("a","b"),c("L.1","L.2"))))
  expect_equal(dim(alleleSplit(x[FALSE,,drop=FALSE],sep="/")),c(0L,2L))
  expect_equal(dim(alleleSplit(x[,FALSE,drop=FALSE],sep="/")),c(2L,0L))
})

test_that("malformed separated calls cannot become homozygotes", {
  for(value in c("T","T/","/T","A/G/T","A//G","")) {
    expect_error(alleleSplit(matrix(c("A/G",value),ncol=1),sep="/"),
                 "two non-empty alleles",fixed=TRUE)
  }
})

test_that("separators are literal and zero alleles remain real", {
  for(sep in c("/","|",".","::")) {
    x <- matrix(c(paste0("A",sep,"G"),paste0("0",sep,"0")),ncol=1)
    expect_equal(unname(alleleSplit(x,sep=sep)),matrix(c("A","0","G","0"),2,2))
  }
})

test_that("legacy fixed width missing codes work", {
  x <- matrix(c("000000","000145","145095",NA,"NANA",""),ncol=1)
  expected <- cbind(c(NA,NA,"145",NA,NA,NA),c(NA,"145","095",NA,NA,NA))
  expect_equal(unname(alleleSplit(x)),expected)
  expect_error(alleleSplit(matrix("123",1,1)),"even width",fixed=TRUE)
  expect_equal(unname(alleleSplit(matrix(c("AG","GG"),ncol=1))),
               matrix(c("A","G","G","G"),2,2))
})
