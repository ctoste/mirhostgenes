MA <- function(x, slide=1, main, control.spots=TRUE){
    if(control.spots){
        idx.spike <- grep(x$genes[, "Name"], pattern="^spike_control")
        idx.spike.ex <- grep(x$genes[, "Name"], pattern="extended_spike")
        idx.neg <- grep(x$genes[, "Name"], pattern="negative")
        bm.empty <- x$genes[, "Name"]=="Empty"
        bm.hy3 <- x$genes[, "Name"]=="Hy3"
        Col.con <- brewer.pal(9, "Set1")[c(3, 4, 1, 9, 8)]
        names(Col.con) <- c("Spike", "Spike.ex", "Neg", "Empty", "Hy3")
    }
    if(missing(main))
        main <- paste("MA plot slide", slide)
    if(class(x)=="RGList"){
        M <- log2(x$R[, slide]) - log2(x$G[, slide])
        A <- 0.5*(log2(x$R[, slide]) + log2(x$G[, slide]))
    }
    else{
        M <- x$M[, slide]
        A <- x$A[, slide]
    }
    plot(3, 3, xlab="A", ylab="M", main=main, xlim=c(min(A), max(A)), ylim=c(min(M), max(M)), pch=NA)
    Colors <- densCols(A, M, colramp=colorRampPalette(rev(brewer.pal(9, "Blues")[ -c(1,2) ])))
#  points(A[ bm.empty ], M[ bm.empty ], col=Col.con[ "Empty" ], cex=0.5, pch=16)
  #Colors[ bm.unknown ] <- Col.con[ "Unknown" ]
    points(A, M, col=Colors, pch=16, cex=0.5)
    lines(lowess(A, M), col="grey")

    if(control.spots){
        ## adding the control features...
        points(A[ bm.empty ], M[ bm.empty ], col=Col.con[ "Empty" ], cex=0.5, pch=16)
        points(A[ idx.spike ], M[ idx.spike ], col=Col.con[ "Spike" ], cex=0.5, pch=16)
        points(A[ idx.spike.ex ], M[ idx.spike.ex ], col=Col.con[ "Spike.ex" ], cex=0.5, pch=16)
        points(A[ idx.neg ], M[ idx.neg ], col=Col.con[ "Neg" ], cex=0.5, pch=16)
        points(A[ bm.hy3 ], M[ bm.hy3 ], col=Col.con[ "Hy3" ], cex=0.5, pch=16)
        legend("top", legend=names(Col.con), col=Col.con, pch=16, horiz=TRUE) ## or ncol=2
    }
}

## values for replicates are always averaged within
averageReps <- function(x, fast=TRUE, id.col="Gene ID"){
    Genes <- x$genes
    UIDs <- unique(Genes[, id.col])
    Genes.unique <- data.frame(matrix(ncol=2, nrow=length(UIDs)))
    colnames(Genes.unique) <- c(id.col, "Name")
    Genes.unique[ , id.col] <- UIDs
    if(class(x)=="RGList"){
        R <- x$R
        G <- x$G
    }
    if(class(x)=="MAList"){
        R <- ma2r(x$M, x$A)
        G <- ma2g(x$M, x$A)
    }
    if(!fast){
        R.unique <- matrix(ncol=ncol(R), nrow=length(UIDs), 0)
        colnames(R.unique) <- colnames(R)
        G.unique <- R.unique
        for(i in 1:length(UIDs)){
            bm <- Genes[ , id.col ]==UIDs[ i ]
            Genes.unique[ i, "Name" ] <- paste(unique(Genes[ bm, "Name" ]), collapse=";")
            if(sum(bm) > 1){
                R.unique[ i,  ] <- apply(R[ bm, ], MARGIN=2, mean)
                G.unique[ i,  ] <- apply(G[ bm, ], MARGIN=2, mean)
            }
            else{
                R.unique[ i,  ] <- R[ bm, ]
                G.unique[ i,  ] <- G[ bm, ]
            }
        }
    }else{
        ## using aggregate...
        tmp <- aggregate(R, by=list(Genes[, id.col]), FUN=mean)
        ## using aggregate on the Genes too...
        tmp.genes <- aggregate(Genes[, "Name"], by=list(Genes[, id.col]),
                               FUN=function(x){
                                   return(paste(unique(x), collapse=";"))
                               })
        colnames(tmp.genes) <- c(id.col, "Name")
        ##idx <- match(tmp[, 1], UIDs)
        idx <- match(UIDs, tmp[, 1])
        R.unique <- as.matrix(tmp[idx, -1])
        tmp <- aggregate(G, by=list(Genes[, id.col]), FUN=mean)
        G.unique <- as.matrix(tmp[idx, -1])
        Genes.unique <- tmp.genes[idx, ]
    }
    X <- x
    X$genes <- Genes.unique
    if(class(x)=="RGList"){
        X$R <- R.unique
        X$G <- G.unique
    }
    if(class(x)=="MAList"){
        X$M <- log2(R.unique)-log2(G.unique)
        X$A <- 0.5*(log2(R.unique)+log2(G.unique))
    }
    return(X)
}

## m=r-g, 2*a=r+g
## r=m+g, g=2*a-r -> r=m+2*a-r -> r=m/2+a, g=2*a-m/2-a=a-m/2
ma2r <- function(m, a){
    return(((m/2)+a))
}
ma2g <- function(m, a){
    return((a-(m/2)))
}

## Test objects.
## Data is from GSE25320.org
##idx <- which(miRNA.mapping$ps_count>0 & miRNA.mapping$pre_mirna_count>1)
##idx2 <- which(miRNA.mapping$ps_count>1)
##mirna.exp <- c( Slides.norm.sub$M[rownames(miRNA.mapping)[1:20],1], Slides.norm.sub$M[idx[1:5],1], Slides.norm.sub$M[idx2[20:25],1])
##names(mirna.exp) <- miRNA.mapping[names(mirna.exp), "Accession"]
#mirna.exp <- mirna.exp[-grep(names(mirna.exp))]


####
## the idea of that function:
## pair miRNA and transcript expression data:
## 1) for each pre-miRNA alignment, select the representative tx or probe set.
## 2) for each pre-miRNA alignment, select the 3p or 5p strand.
## 3) for each mature miRNA, select the pre-miRNA if it is encoded in more than one.
## What is still lacking is a way to select one pre- or mature miRNA for a miRNA cluster,
## or at least get that information.
## Details:
## we do have "special cases":
## a) more than one probe set or transcript per pre-miRNA alignment (allows to choose
##    the probe set id or tx using the chooseTxFun).
## b) two mature miRNAs per pre-miRNA (allows to choose the miRNA strand using the
##    chooseMatMirFun function).
## c) more than one pre-miRNA per mature miRNA (allows to choose the pre-miRNA-host tx
##    combination using the choosePreMirFun function).
    ## d) can not find mature miRNA (we just don't return anything).
## e) no host tx mapped to pre-miRNA (we just don't return anything).
## f) no probe set for host tx (we just don't return anything).
##
## object: a MirhostDb object
## mirna: named numeric vector with mature miRNA expression values
## tx: named numeric vector with host transcript expression values
## mirnaNamesAre: what are the names of mirna: mat_mirna_id, mat_mirna_name
## txNamesAre: what are the names of tx: probeset_id or tx_id
## array the microarray type
##
## chooseXXXFun: a function that chooses which XXX to be selected. See chooseAll, chooseOnValue
##               and chooseClosest
## returns a data.frame with:
pairData <- function(object, mirna, tx, mirnaNamesAre="mat_mirna_name", txNamesAre="probeset_id",
                     chooseTxFun=chooseAll, chooseMatMirFun=chooseAll, choosePreMirFun=chooseAll
                     ){
    ## TODO: allow tx and mirna to be also data.frame or matrix with rownames.
    if(missing(mirna) | missing(tx))
        stop("Both, mirna and tx are required!")
    if(!class(mirna) %in% c("numeric", "integer"))
        stop("mirna has to be a numeric vector!")
    if(!class(tx) %in% c("numeric", "integer"))
        stop("tx has to be a numeric vector!")
    if(is.null(names(mirna)))
        stop("mirna has to be a named vector with the names corresponding to mature miRNA identifiers!")
    if(is.null(names(tx)))
        stop("tx has to be a named vector with the names corresponding to either Ensembl transcript IDs or probe set IDs!")
    mirnaNamesAre <- match.arg(mirnaNamesAre, c("mat_mirna_name", "mat_mirna_id"))
    ## define the filter to be used for the miRNA...
    if(mirnaNamesAre == "mat_mirna_name"){
        MirFilter <- MatmirnaFilter
    }else{
        MirFilter <- MatmirnaidFilter
    }
    txNamesAre <- match.arg(txNamesAre, c("tx_id", "probeset_id"))
    ##returnMirna <- match.arg(c("mat_mirna", "pre_mirna"))
    ## check if array exists if we're going to query probe sets
    Filts <- NULL
    ## which column should be used for splitting of the pre-miRNA later on?
    premir <- "pre_mirna_algn_id"
    ## map the mirna names to pre_mirnas AND txNamesAre
    mirna2txFun <- function(x){
        res <- matmirnas(object, columns=unique(c(premir, "pre_mirna_name",
                                     "pre_mirna_id", txNamesAre)),
                         filter=list(MirFilter(x)), return.type="data.frame")
        return(res)
    }
    mir2tx <- mirna2txFun(x=names(mirna))
    if(nrow(mir2tx)==0)
        return(mir2tx)
    ## count for each mature miRNA the number of pre-miRNAs in which it is
    ## encoded.
    splittedMir2tx <- split(mir2tx, f=mir2tx[, mirnaNamesAre])
    mir2tx <- do.call(rbind,
                      lapply(splittedMir2tx, function(z){
                                 return(
                                     data.frame(z,
                                           pre_mir_count=rep(length(unique(z[, premir])),
                                               nrow(z)))
                                     )
                                 })
                      )
    ## rownames are confusing...
    rownames(mir2tx) <- NULL
    mir2tx <- mir2tx[which(mir2tx[, txNamesAre] %in% names(tx)),]
    if(nrow(mir2tx)==0)
        return(mir2tx)
    ## fill with data...
    mir2tx <- cbind(mir2tx, mir_exp=mirna[mir2tx[, mirnaNamesAre]],
                    tx_exp=tx[mir2tx[, txNamesAre]])
    ## 1) split by pre-miRNA and eventually reduce a) multiple probe sets/
    ##    pre-miRNAs to one per pre-miRNA and b) mature miRNAs, i.e.
    ##    select one of the 3p and 5p strand.
    splittedMir2tx <- split(mir2tx, f=mir2tx[, premir])
    mir2tx <- do.call(rbind,
                      lapply(splittedMir2tx, function(z){
                                 if(nrow(z)==1)
                                     return(z)
                                 ## choose the "best" probe set/tx; since we can have
                                 ## the 5' and the 3' miRNA here it is bettern to return
                                 ## all entries that are identified with the probe set or
                                 ## tx id
                                 id <- z[chooseTxFun(z[, c("tx_exp", "mir_exp")]), txNamesAre]
                                 z <- z[z[, txNamesAre] %in% id, ]
                                 ## next we can return reduce for mature miRNAs
                                 id <- z[chooseMatMirFun(z[, c("mir_exp", "tx_exp")]), mirnaNamesAre]
                                 return(z[z[, mirnaNamesAre] %in% id, ])
                             })
                      )
    if(nrow(mir2tx)==0)
        return(mir2tx)
    ## 2) what if a mature miRNA is encoded by several pre-miRNAs? Split
    ##    by mirnaNamesAre and call the choosePreMirFun on that.
    splittedMir2tx <- split(mir2tx, f=mir2tx[, mirnaNamesAre])
    mir2tx <- do.call(rbind,
                      lapply(splittedMir2tx, function(z){
                                 if(nrow(z)==1)
                                     return(z)
                                 return(z[choosePreMirFun(z[, c("tx_exp", "mir_exp")]), ])
                             })
                      )
    rownames(mir2tx) <- NULL
    return(mir2tx)
}



#### choose functions:
## these functions take a matrix or data.frame as input and are supposed
## to return numeric values representing the index of the values in the
## submitted x that were selected.


## Takes a matrix or data.frame and returns a sequence from 1:nrow(x), i.e.
## selects all.
chooseAllOld <- function(x){
    ##return(rownames(x))
    return(seq(1, nrow(x)))
}
## Takes a matrix or data.frame and returns the first element of the
## result from the orderFun. The orderFun will be applied to the FIRST! column
## of the submitted matrix/data.frame.
chooseOnValueOld <- function(x, orderFun=function(z){order(z, decreasing=TRUE)}){
    if(nrow(x)==1){
        return(1)
    }
    return(orderFun(x[, 1])[1])
}

## Takes a matrix or data.frame calculates the difference between the elements
## in columns 1 and 2 and returns the index of the smallest difference.
chooseClosestOld <- function(x){
    Diffs <- x[,1] - x[,2]
    return(order(abs(Diffs))[1])
}



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

mirna <- miRNA.exprs
tx <- tx.exprs[, "exprs", drop=TRUE]
names(tx) <- rownames(tx.exprs)
mirnaNamesAre <- "mat_mirna_name"
txNamesAre <- "probeset_id"
chooseTxFun <- chooseAll
chooseMatMirFun <- chooseAll
choosePreMirFun <- chooseAll

## testing the function below...
## library(MirhostDb.Hsapiens.v75.v20)
##MhgDb <- MirhostDb.Hsapiens.v75.v20
library(MirhostDb.Hsapiens.v81.v21)
MhgDb <- MirhostDb.Hsapiens.v81.v21
## get all
All <- pairData(MhgDb, mirna=mirna, tx=tx, chooseTxFun=chooseAll, chooseMatMirFun=chooseAll,
                choosePreMirFun=chooseAll)
All
## OK
## get all mat miRNAs, all pre-miRNAs but only one tx per pre-miRNA
OneTx <- pairData(MhgDb, mirna=mirna, tx=tx, chooseTxFun=chooseOnValue,
                  chooseMatMirFun=chooseAll, choosePreMirFun=chooseAll)
## OK
## get one mat miRNA, all pre-miRNAs and one tx per per-miRNA
OneMatOneTx <- pairData(MhgDb, mirna=mirna, tx=tx, chooseTxFun=chooseOnValue,
                        chooseMatMirFun=chooseOnValue, choosePreMirFun=chooseAll)
## OK
## get the mat miRNA just from one pre-miRNA
OnePre <- pairData(MhgDb, mirna=mirna, tx=tx, chooseTxFun=chooseAll,
                   chooseMatMirFun=chooseAll, choosePreMirFun=chooseOnValue)
## NOPE!!! does reduce also per tx!!!
