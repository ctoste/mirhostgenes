## testing the value matching...
test_doPairData <- function(){
    ## Simple test:
    A <- 1:10
    names(A) <- c("b", "a", "c", "b", "a", "a", "e", "f", "a", "c")

    B <- c(2, 5, 3, 2, 5, 8, 4, 20)
    names(B) <- c("a", "b", "a", "c", "c", "c", "n", "f")

    Res <- doPairData(A, B, chooseFunX=chooseAll, chooseFunY=chooseAll)
    expectx <- c(1, 4, 2, 2, 5, 5, 6, 6, 9, 9, 3, 3, 3, 10, 10, 10, 8, 7, NA)
    expecty <- c(5, 5, 2, 3, 2, 3, 2, 3, 2, 3, 2, 5, 8, 2, 5, 8, 20, NA, 4)
    checkEquals(Res$x, expectx)
    checkEquals(Res$y, expecty)

    ## now let's use chooseAll for x and chooseOrderedValue for y
    expectx <- c(1, 4, 2, 5, 6, 9, 3, 10, 8, 7, NA)
    expecty <- c(5, 5, 3, 3, 3, 3, 8, 8, 20, NA, 4)
    Res <- doPairData(A, B, chooseFunY=chooseOrderedValue)
    checkEquals(Res$x, expectx)
    checkEquals(Res$y, expecty)

    ## same but for chooseOrderedValue for x:
    Res <- doPairData(A, B, chooseFunX=chooseOrderedValue)
    expectx <- c(4, 9, 9, 10, 10, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 2, 3, 2, 5, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## now with chooseOrderedValue for both.
    Res <- doPairData(A, B, chooseFunX=chooseOrderedValue, chooseFunY=chooseOrderedValue)
    expectx <- c(4, 9, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 3, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## ## OK, now with closest value for y, all for x
    ## ## what does this: selects all values from A and selects the one value from B for
    ## ## each matching name that is closest to any of the values with the same name in A
    ## ## Not very meaningfull though...
    ## Res <- doPairData(A, B, chooseFunY=chooseClosestBetween)
    ## expectx <- c(1, 4, 2, 5, 6, 9, 3, 10, 8, 7, NA)
    ## checkEquals(Res$x, expectx)
    ## expecty <- c(5, 5, 2, 2, 2, 2, 2, 2, 20, NA, 4)
    ## checkEquals(Res$y, expecty)

    ## ## reverse thing:
    ## ## first select the one value from A closest to any of the values in B with the same
    ## ## name and then match that to all values in B; also not that meaningful
    ## Res <- doPairData(A, B, chooseFunX=chooseClosestBetween)
    ## expectx <- c(4, 2, 2, 3, 3, 3, 8, 7, NA)
    ## checkEquals(Res$x, expectx)
    ## expecty <- c(5, 2, 3, 2, 5, 8, 20, NA, 4)
    ## checkEquals(Res$y, expecty)

    ## OK, more meaningful:
    ## all in x and for each in x the closest from y with the same name.
    ## fails with second xselect
    Res <- doPairData(A, B, chooseFunY=chooseClosestForEach)
    expectx <- c(1, 4, 2, 5, 6, 9, 3, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 5, 2, 3, 3, 3, 2, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## reverse
    ## all in y and choose for each the closest in x
    Res <- doPairData(A, B, chooseFunX=chooseClosestForEach)
    expectx <- c(4, 2, 2, 3, 3, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 2, 3, 2, 5, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## ## now something meaningful:
    ## ## select first the largest value in x and then select the closest in y
    ## ## NOT meaningful, since it depends on the order the functions are applied!
    ## Res <- doPairData(A, B, chooseFunX=chooseOrderedValue, chooseFunY=chooseClosestBetween)
    ## expectx <- c(4, 9, 10, 8, 7, NA)
    ## checkEquals(Res$x, expectx)
    ## expecty <- c(5, 3, 8, 20, NA, 4)
    ## checkEquals(Res$y, expecty)

    ## ## not meaningful... other way round.
    ## ## select first the value best matching in y and then the largest in y
    ## ## NOT meaningful, since it depends on the order the functions are applied!
    ## Res <- doPairData(A, B, chooseFunX=chooseClosestBetween, chooseFunY=chooseOrderedValue)
    ## expectx <- c(4, 2, 3, 8, 7, NA)
    ## checkEquals(Res$x, expectx)
    ## expecty <- c(5, 3, 8, 20, NA, 4)
    ## checkEquals(Res$y, expecty)

    ## use the chooseClosestForEach insted... should give the same result than chooseClosestBetween
    ## i.e. first select the largest value in x and then select the closest from y matching
    ## THAT selected value in x.
    Res <- doPairData(A, B, chooseFunX=chooseOrderedValue, chooseFunY=chooseClosestForEach)
    expectx <- c(4, 9, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 3, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## reverse order; should also make sense!
    ## select for each in y the closest in x, then select the largest in y.
    ## Note: this is different than using chooseOrderedValue for x and chooseClosestForEach for y
    Res <- doPairData(A, B, chooseFunX=chooseClosestForEach, chooseFunY=chooseOrderedValue)
    expectx <- c(4, 2, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 3, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## now, chooseClosestForEach for both: selects the closest pairs of data.
    ## this selects first the closest in x to the values in y, then, for all in y the
    ## closest in x (already subsetted to those that are close to values in y).
    Res <- doPairData(A, B, chooseFunX=chooseClosestForEach, chooseFunY=chooseClosestForEach)
    ## b: A x(1, 4), B y(5):
    ##    chooseX: (4)
    ##    chooseY: (5)
    ## a: A x(2, 5, 6, 9), B y(2, 3)
    ##    chooseX: (2, 2)
    ##    chooseY: (2, 2)
    ## c: A x(3, 10), B y(2, 5, 8)
    ##    chooseX: (3, 3, 10)
    ##    chooseY: (2, 2, 8)
    ## f: A x(8), B y(20)
    ##    chooseX: (8)
    ##    chooseY: (20)
    ## e: A x(7), B y(NA)
    ## n: A x(NA), B y(4)
    ## At last unique is called on the resulting data.frame that reduced the values for
    ## a and c.
    expectx <- c(4, 2, 3, 10, 8, 7, NA)
    checkEquals(Res$x, expectx)
    expecty <- c(5, 2, 2, 8, 20, NA, 4)
    checkEquals(Res$y, expecty)

    ## second example:
    C <- c(3, 8, 7, 1)
    names(C) <- rep("a", length(C))
    D <- c(9, 2, 4, 20, 3)
    names(D) <- rep("a", length(D))
    Res <- doPairData(C, D, chooseFunX=chooseClosestForEach,
                      chooseFunY=chooseClosestForEach)

    ## same length for some
    A <- 1:4
    names(A) <- c("a", "a", "b", "b")
    Res <- doPairData(A, A)
    checkEquals(Res$x, rep(1:4, each=2))
    checkEquals(Res$y, c(1, 2, 1, 2, 3, 4, 3, 4))
}

notrun_test_doMatchMir <- function(){
    library(MirhostDb.Hsapiens.v81.v21)
    object <- MirhostDb.Hsapiens.v81.v21
    ## that's now for miRNA matching...
    ## Test data.
    ## have mat miRNAs from pre-miRNAs miR-15b, 16-1 and 223
    miRNA.exprs <- c(8, 8.3, 5.6, 9.5, 4.6, 13.1)
    names(miRNA.exprs) <- c("hsa-miR-15b-3p", "hsa-miR-15b-5p", "hsa-miR-16-5p", "hsa-miR-16-1-3p", "hsa-miR-223-3p", "hsa-miR-223-5p")
    ## expression data
    tx.exprs <- data.frame(tx_id=c("NR_002612.1", "ENST00000344722", "ENST00000344722",
                                   "ENST00000468653", "ENSESTT00000026526"),
                           probeset_id=c("1564443_at", "201663_s_at", "201664_at",
                                         "215623_x_at", "229934_at"),
                           pre_mirna=c("hsa-mir-16-1", "hsa-mir-15b/16-2", "hsa-mir-15b/16-2",
                                       "hsa-mir-15b/16-2", "hsa-mir-223"),
                           exprs=c(8.9, 5.4, 8.8, 12.1, 4.3))
    rownames(tx.exprs) <- tx.exprs$probeset_id

    xNamesAre <- "mat_mirna_name"
    yNamesAre <- "probeset_id"
    suppNames <- c("pre_mirna_id", "pre_mirna_name", "tx_id", "mat_mirna_id", "mat_mirna_name",
                   "gene_id", "gene_name", "probeset_id")
    ## define the table name for each...
    names(suppNames) <- c("pre_mirna", "pre_mirna", "host_tx", "mat_mirna", "mat_mirna",
                          "host_gene", "host_gene", "array_feature")
    xNamesAre <- match.arg(xNamesAre, suppNames)
    yNamesAre <- match.arg(yNamesAre, suppNames)
    ## just check if we do have probesets available... if needed...
    if(xNamesAre == "probeset_id" | yNamesAre == "probeset_id"){
        if(!object@have_array_features)
            stop("No probe set definitions available in the database! Thus 'probeset_id' can not be used for arguments 'xNamesAre' or 'yNamesAre'.")
    }
    ## use .getWhat (without filter), but need to define start.table...
    allData <- mirhostgenes:::.getWhat(object, columns=c(xNamesAre, yNamesAre),
                                       start.table=names(suppNames)[suppNames == xNamesAre])
    ## OK, that's now the mapping between names in x and names in y.
    ## check for how many in x and in y we don't have values...
    ## I run into a problem here... same-named values is nice, but what if we have a situation like here?
    head(allData)
    ## assume, have a value for let-7a-5p, and one each for the probe sets.
    ## 1) map first to pre-miRNAs? i.e. select for each pre-miRNA the one mat_mirna, or select
    ##    all... use doSelectData using the pre-miRNA id replacing the mat_mirna.
    ## mapping:
    ## mat_miRNA -> host_tx/probe set:
    ## for each pre-miRNA, select mat_miRNAs (all or just one.)
    ## common name is the pre-miRNA!

    ## Test it!
    ## first get from mat miRNA to pre-miRNA
    mat2pre <- mirhostgenes:::.getWhat(object, columns=c(xNamesAre, "pre_mirna_id"),
                                       start.table=names(suppNames)[suppNames == xNamesAre])
    mat2pre <- mat2pre[mat2pre[, xNamesAre] %in% names(miRNA.exprs), ]
    Test <- miRNA.exprs[mat2pre[, xNamesAre]]
    names(Test) <- mat2pre[, "pre_mirna_id"]
    pres <- doSelectData(Test, chooseOrderedValue)
    predata <- pres$x
    names(predata) <- pres$name
    ## that way I would have values for pre-miRNAs, that I could match then to host genes.
    ## mapping pre-miRNAs to probe sets is n:m
    pre2ps <- mirhostgenes:::.getWhat(object, columns=c("pre_mirna_id", "probeset_id"),
                                      start.table="pre_mirna")
    pre2ps <- pre2ps[pre2ps$pre_mirna_id %in% names(predata), ]
}


test_transferValues <- function(){
    ## load data
    library(MirhostDb.Hsapiens.v75.v20)
    mhg <- MirhostDb.Hsapiens.v75.v20
    ## define values
    miRNA.exprs <- c(8, 8.3, 5.6, 9.5, 4.6, 13.1)
    names(miRNA.exprs) <- c("hsa-miR-15b-3p", "hsa-miR-15b-5p", "hsa-miR-16-5p",
                            "hsa-miR-16-1-3p", "hsa-miR-223-3p", "hsa-miR-223-5p")
    ## map mat miRNAs to pre-miRNAs, select all.
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
                          toNames="pre_mirna_name")
    ResX <- Res$x
    names(ResX) <- Res$mat_mirna_name
    checkEquals(miRNA.exprs, ResX[names(miRNA.exprs)])
    ## what if we map to host genes?
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
                          toNames="gene_name")
    ResX <- Res$x
    names(ResX) <- Res$mat_mirna_name
    checkEquals(miRNA.exprs, ResX[names(miRNA.exprs)])
    ## use a solveFun
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
                          toNames="pre_mirna_name", solveFun=chooseOrderedValue)
    ResX <- Res$x
    names(ResX) <- Res$mat_mirna_name
    checkEquals(miRNA.exprs[names(ResX)], ResX)
    ## use a solveFun, host gene
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name",
                          toNames="gene_name", solveFun=chooseOrderedValue)
    ResX <- Res$x
    names(ResX) <- Res$mat_mirna_name
    checkEquals(miRNA.exprs[names(ResX)], ResX)

    ## now: use a seqstrand filter, that should only return values for pre-miRNA-16-2
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name", toNames="pre_mirna_name",
                          filter=SeqstrandFilter(1), na.rm=TRUE)
    checkEquals(any(Res$pre_mirna_name == "hsa-mir-16-1"), FALSE)

    ## check if the filter works for gene_biotype
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name", toNames="gene_id",
                          na.rm=TRUE, filter=DatabaseFilter("core"))
    Res <- hostgenes(mhg, filter=GeneidFilter(Res$gene_id))
    checkEquals(unique(Res$database), "core")
    ## biotype
    Res <- transferValues(miRNA.exprs, mhg, xNamesAre="mat_mirna_name", toNames="gene_id",
                          na.rm=TRUE, filter=GenebiotypeFilter("miRNA", condition="!="))
    Res <- hostgenes(mhg, filter=GeneidFilter(Res$gene_id))
    checkEquals(any(Res$gene_biotype == "miRNA"), FALSE)
}



