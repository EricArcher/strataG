context("genind conversion safeguards")

test_that("diploid genotypes, names and populations survive conversion", {
  x <- adegenet::df2genind(data.frame(a=c("1/1","1/2",NA,"2/2"),
                                    b=c("2/2","1/2","1/1","1/2")),
                          sep="/",ind.names=c("c","a","b","d"),
                          pop=c("B","A","B","A"))
  adegenet::locNames(x) <- c("marker.1","marker_1")
  g <- genind2gtypes(x)
  expect_setequal(getLociNames(g), adegenet::locNames(x))
  z <- gtypes2genind(g)
  expect_setequal(adegenet::locNames(z), adegenet::locNames(x))
  expect_setequal(adegenet::indNames(z), adegenet::indNames(x))
  expect_equal(z@tab[rownames(x@tab),colnames(x@tab)],x@tab)
  expect_equal(as.character(adegenet::pop(z))[match(adegenet::indNames(x),adegenet::indNames(z))],
               as.character(adegenet::pop(x)))
  expect_true(methods::validObject(z))
})

test_that("haploid locus names retain no extra suffix", {
  x <- adegenet::df2genind(data.frame(L=c("1","2","1")),ploidy=1)
  g <- genind2gtypes(x)
  expect_identical(getLociNames(g),"L")
  expect_equal(getPloidy(g),1)
  expect_equal(gtypes2genind(g)@tab,x@tab)
})

test_that("mixed ploidy and presence absence are rejected", {
  x <- adegenet::df2genind(data.frame(L=c("1","1/2","2/2")),
                          sep="/",ploidy=c(1,2,2))
  expect_error(genind2gtypes(x),"mixed ploidy",fixed=TRUE)
  pa <- adegenet::df2genind(data.frame(a=c(0,1,0),b=c(1,0,1)),type="PA",ncode=1)
  expect_error(genind2gtypes(pa),"Presence/absence",fixed=TRUE)
  g <- df2gtypes(data.frame(id=c("a","b"),pop="P",L=c(0,1)),ploidy=1)
  expect_error(gtypes2genind(g,type="PA"),"Presence/absence",fixed=TRUE)
})
