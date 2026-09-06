context("loci conversion safeguards")

loci_fixture <- function() {
  df2gtypes(data.frame(id=c("a","b","c"),pop=c("A","B","A"),
    L.1=c("1","2",NA),L.2=c("1","2",NA),
    M.1=c("0","1","1"),M.2=c("0","2","2")),ploidy=2)
}

test_that("missing calls stay missing and real zero alleles survive", {
  g <- loci_fixture()
  l <- gtypes2loci(g)
  expect_true(is.na(l[["L"]][3]))
  expect_false("0/0" %in% levels(l[["L"]]))
  expect_true("0/0" %in% levels(l[["M"]]))
  z <- loci2gtypes(l)
  expect_true(all(is.na(z@data$allele[z@data$id=="c" & z@data$locus=="L"])))
  expect_identical(z@data$allele[z@data$id=="a" & z@data$locus=="M"],c("0","0"))
})

test_that("unknown populations do not discard genotyped samples", {
  g <- loci_fixture()
  g@data$stratum[g@data$id=="b"] <- NA_character_
  l <- gtypes2loci(g)
  expect_identical(rownames(l),c("a","b","c"))
  expect_true(is.na(l$population[2]))
  z <- loci2gtypes(l)
  expect_true(all(is.na(z@data$stratum[z@data$id=="b"])))
})

test_that("population column is located rather than assumed", {
  l <- pegas::as.loci(data.frame(L=factor(c("1/1","1/2")),
                       population=factor(c("A","B"))),col.pop=2)
  z <- loci2gtypes(l)
  expect_identical(unname(getStrata(z)),c("A","B"))
  l <- pegas::as.loci(data.frame(L=factor(c("1/1","1/2"))))
  expect_true(all(getStrata(loci2gtypes(l))=="Default"))
})

test_that("literal separators round trip without changing alleles", {
  g <- loci_fixture()
  for(sep in c(":","|",".")) {
    z <- loci2gtypes(gtypes2loci(g,sep=sep),sep=sep)
    expect_equal(z@data,g@data)
  }
})

test_that("unsupported ploidy is rejected explicitly", {
  g <- df2gtypes(data.frame(id=c("a","b"),pop="P",L=c("A","G")),ploidy=1)
  expect_error(gtypes2loci(g),"Only diploid",fixed=TRUE)
  l <- pegas::as.loci(data.frame(L=factor(c("A","G"))))
  expect_error(loci2gtypes(l),"Only diploid",fixed=TRUE)
})
