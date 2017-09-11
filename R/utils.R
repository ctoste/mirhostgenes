## Some small utility functions.

## checks whether x is a class of type BasicFilter or a list
## of BasicFilters.
## returns a list of BasicFilter
checkIsBasicFilter <- function(x){
    if(class(x) != "list")
        x <- list(x)
    notBF <- unlist(lapply(x, function(z){
        !inherits(z, what="AnnotationFilter")
    }))
    if(any(notBF))
        stop("filter should be a list of filter objects or a single filter object (i.e. an object extending AnnotationFilter)!")
    return(x)
}





## generic version...
## input: x, y are two named numeric vectors with names at least partially matching.
## The aim of the function is to match entries in input vector x with values in input
## vector y based on their names, i.e. values with the same name are matched with each
## other and returned in the same row in the result table. If x and/or y each contains
## entries with the same name, values are repeated such that each value in x is matched
## to each value in y with the same name (i.e. if x contains 2 values with the name a,
## and y 3 with that name, the function matches each of the 2 values in x with each of the
## 3 values in y resulting in 6 rows in the result table for name a). This default
## behaviour can be changed by specifying a function different than chooseAll for argument
## chooseFunX or chooseFunY.
## The function first applies the chooseFunX and chooseFunY on vectors x and y to
## eventually select one value for each elements with the same name (by default, i.e.
## function chooseAll all are selected). Next the values are matched between vectors,
## i.e. each value in x and y are matched with each other by eventually repeating values
## IMPORTANT NOTE: chooseFunX is applied before chooseFunY!
doPairData <- function(x, y, chooseFunX=chooseAll, chooseFunY=chooseAll){
    if(!is.numeric(x) | !is.numeric(y))
        stop("x and y have to be numeric!")
    if(is.null(names(x)) | is.null(names(y)))
        stop("x and y have to be named numeric vectors!")
    ## names have to match...
    comNames <- intersect(names(x), names(y))
    ## looks tedious and it is. all just because we want to return
    ## also the index of the value from each vector.
    xSplit <- split(data.frame(x, idx=1:length(x), name=names(x),
                               stringsAsFactors=FALSE), names(x))
    ySplit <- split(data.frame(y, idx=1:length(y), name=names(y),
                               stringsAsFactors=FALSE), names(y))
    xSplitProc <- lapply(xSplit[comNames], function(zx){
        zy <- ySplit[[unique(zx$name)]]
        ## repeat entries so that lengths match between zx and zy.
        ## that way we also ensure to map each with each for the
        ## chooseAll case
        lenzy <- nrow(zy)
        lenzx <- nrow(zx)
        zx <- zx[rep(1:nrow(zx), each=lenzy), ]
        zy <- zy[rep(1:nrow(zy), times=lenzx), ]
        ## OK, let's select values based in chooseFun:
        idxzx <- chooseFunX(zx$x, zy$y)
        ## Next the same for zy, given we did something to zx
        idxzy <- chooseFunY(zy$y, zx$x[idxzx])
        ## select x idx a second time, given idxzx and idxzy.
        idxzx <- idxzx[chooseFunX(zx$x[idxzx], zy$y[idxzy])]
        if(length(idxzx) != length(idxzy)){
        ## ALWAYS do that; it can be that both length are the same by chance (same input value)
            lenzy <- length(idxzy)
            lenzx <- length(idxzx)
            idxzx <- idxzx[rep(1:length(idxzx), each=lenzy)]
            idxzy <- idxzy[rep(1:length(idxzy), times=lenzx)]
        }
        ## that way we can merge...
        tmpx <- zx[idxzx, ]
        tmpy <- zy[idxzy, ]
        rownames(tmpx) <- NULL
        rownames(tmpy) <- NULL
        return(cbind(x=tmpx, y=tmpy,
                     stringsAsFactors=FALSE))
    })
    ## add the remaining ones.
    xSplit <- xSplit[!(names(xSplit) %in% comNames)]
    if(length(xSplit) > 0){
        xSplit <- lapply(xSplit, function(z){
            z <- cbind(z, y.y=NA, y.idx=NA, y.name=NA)
            colnames(z)[1:3] <- c("x.x", "x.idx", "x.name")
            return(z)
        })
        xSplitProc <- c(xSplitProc, xSplit)
    }
    ySplit <- ySplit[!(names(ySplit) %in% comNames)]
    if(length(ySplit) > 0){
        ySplit <- lapply(ySplit, function(z){
            z <- cbind(x.x=NA, x.idx=NA, x.name=NA, z)
            colnames(z)[4:6] <- c("y.y", "y.idx", "y.name")
            return(z)
        })
        xSplitProc <- c(xSplitProc, ySplit)
    }
    res <- do.call(rbind, xSplitProc)
    Names <- apply(res[, c("x.name", "y.name")], MARGIN=1, function(nn){
        return(unique(unlist(nn[!is.na(nn)])))
    })
    res <- cbind(name=Names, res[, -grep(colnames(res), pattern="name")])
    colnames(res) <- c("name", "x", "x.idx", "y", "y.idx")
    res <- unique(res)
    rownames(res) <- NULL
    return(res)
}

doPairDataOld <- function(x, y, chooseFunX=chooseAll, chooseFunY=chooseAll){
    if(!is.numeric(x) | !is.numeric(y))
        stop("x and y have to be numeric!")
    if(is.null(names(x)) | is.null(names(y)))
        stop("x and y have to be named numeric vectors!")
    ## names have to match...
    comNames <- intersect(names(x), names(y))
    ## looks tedious and it is. all just because we want to return
    ## also the index of the value from each vector.
    xSplit <- split(data.frame(x, idx=1:length(x), name=names(x),
                               stringsAsFactors=FALSE), names(x))
    ySplit <- split(data.frame(y, idx=1:length(y), name=names(y),
                               stringsAsFactors=FALSE), names(y))
    xSplitProc <- lapply(xSplit[comNames], function(zx){
        zy <- ySplit[[unique(zx$name)]]
        ## got now values z and z2...
        ## OK, now we're matching stuff...
        ## repeat values to match numbers between values
        ## each value in x (i.e. zx) should be matched to
        ## each value in y (i.e. zy). create indices that allow
        ## to extract values for such matched cases.
        idxzx <- rep(1:nrow(zx), times=nrow(zy))  ## use that to extract values from zx
        idxzy <- rep(1:nrow(zy), each=nrow(zx))  ## use that to extract values from zy
        ## first: select value in zx, i.e. select value for x
        idxzx <- idxzx[chooseFunX(zx[idxzx, "x"], zy[idxzy, "y"])]
        if(length(idxzx) != length(idxzy)){
            ## re-create the idxzy; that keeps it updated with idxzx
            idxzy <- rep(1:nrow(zy), each=length(idxzx))
        }
        ## second: select value in z2
        idxzy <- idxzy[chooseFunY(zy[idxzy, "y"], zx[idxzx, "x"])]
        ## OK, so, now we have the indices for the guys we want to match. The chooseFun
        ## is expected to select either all or only one value, thus only if different
        ## chooseFuns were selected for x and y we can have indexes of different length
        ## (supposedly either the same length or 1), and only in that case we have to
        ## repeat one of them to match values again.
        if(length(idxzx) != length(idxzy)){
            ## nope, rep on the index, not on the data.frame!!!
            idxzx <- idxzx[rep(1:length(idxzx), times=length(idxzy))]
            idxzy <- idxzy[rep(1:length(idxzy), each=length(idxzx))]
            ## idxzx <- idxzx[rep(1:length(idxzx), times=length(idxzy))]
            ## idxzy <- idxzy[rep(1:length(idxzy), times=length(idxzx))]
        }
        ## that way we can merge...
        tmpx <- zx[idxzx, ]
        tmpy <- zy[idxzy, ]
        rownames(tmpx) <- NULL
        rownames(tmpy) <- NULL
        return(cbind(x=tmpx, y=tmpy,
                     stringsAsFactors=FALSE))
    })
    ## add the remaining ones.
    xSplit <- xSplit[!(names(xSplit) %in% comNames)]
    if(length(xSplit) > 0){
        xSplit <- lapply(xSplit, function(z){
            z <- cbind(z, y.y=NA, y.idx=NA, y.name=NA)
            colnames(z)[1:3] <- c("x.x", "x.idx", "x.name")
            return(z)
        })
        xSplitProc <- c(xSplitProc, xSplit)
    }
    ySplit <- ySplit[!(names(ySplit) %in% comNames)]
    if(length(ySplit) > 0){
        ySplit <- lapply(ySplit, function(z){
            z <- cbind(x.x=NA, x.idx=NA, x.name=NA, z)
            colnames(z)[4:6] <- c("y.y", "y.idx", "y.name")
            return(z)
        })
        xSplitProc <- c(xSplitProc, ySplit)
    }
    res <- do.call(rbind, xSplitProc)
    Names <- apply(res[, c("x.name", "y.name")], MARGIN=1, function(nn){
        return(unique(unlist(nn[!is.na(nn)])))
    })
    res <- cbind(name=Names, res[, -grep(colnames(res), pattern="name")])
    colnames(res) <- c("name", "x", "x.idx", "y", "y.idx")
    res <- unique(res)
    rownames(res) <- NULL
    return(res)
}

## choose all values
## x: numeric vector of length xlen
## y: numeric vector of length ylen
## xlen and ylen do not have to match. the
## returns an index
chooseAll <- function(x, y){
    return(seq(1, length(x)))
}

## orders the values in x and returns the index of the first.
chooseOrderedValue <- function(x, y, orderFun=function(z){order(z, decreasing=TRUE)}){
    if(length(x)==1){
        return(1)
    }
    return(orderFun(x)[1])
}

## selects the index of the value in x with the smallest difference to any of the
## values in y.
chooseClosestBetween <- function(x, y, abs=TRUE){
    Diffs <- x-y
    if(abs){
        Diffs <- abs(Diffs)
    }
    return(order(Diffs)[1])
}

## returns the index in x which is closest to each value in y.
## resulting vector is thus equal to the length of y!
## selects for each entry in y the closest value in x
chooseClosestForEach <- function(x, y, abs=TRUE){
    xidx <- sapply(y, function(z){
        Diffs <- x-z
        if(abs){
            Diffs <- abs(Diffs)
        }
        return(order(Diffs)[1])
    })
    return(xidx)
}

## basically, that takes a large input and selects
## x: named vector.
doSelectData <- function(x, chooseFunX=chooseAll){
    if(is.null(names(x)))
        stop("'x' has to be a named numeric vector!")
    xSplit <- split(data.frame(x=x, x.idx=1:length(x), name=names(x),
                               stringsAsFactors=FALSE), names(x))
    xSplit <- lapply(xSplit, function(y){
        return(y[chooseFunX(y$x), ])
    })
    res <- do.call(rbind, xSplit)
    rownames(res) <- NULL
    res <- res[order(res$x.idx),]
    return(res)
}

