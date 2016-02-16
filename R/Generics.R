##***********************************************************************
##
##     Generic methods
##
##***********************************************************************
if(!isGeneric("condition"))
    setGeneric("condition", function(x, ...)
        standardGeneric("condition"))
if(!isGeneric("hosttx")){
    setGeneric("hosttx", function(x, ...)
        standardGeneric("hosttx"))
}
if(!isGeneric("hosttxBy")){
    setGeneric("hosttxBy", function(x, ...)
        standardGeneric("hosttxBy"))
}
if(!isGeneric("hostgenes")){
    setGeneric("hostgenes", function(x, ...)
        standardGeneric("hostgenes"))
}
if(!isGeneric("hostgenesBy")){
    setGeneric("hostgenesBy", function(x, ...)
        standardGeneric("hostgenesBy"))
}
if(!isGeneric("listArrays")){
    setGeneric("listArrays", function(x, ...)
        standardGeneric("listArrays"))
}
if(!isGeneric("listColumns")){
    setGeneric("listColumns", function(x, ...)
        standardGeneric("listColumns"))
}
if(!isGeneric("listDatabases")){
    setGeneric("listDatabases", function(x, ...)
        standardGeneric("listDatabases"))
}
if(!isGeneric("listTables")){
setGeneric("listTables", function(x, ...)
    standardGeneric("listTables"))
}
if(!isGeneric("matmirnas")){
    setGeneric("matmirnas", function(x, ...)
        standardGeneric("matmirnas"))
}
if(!isGeneric("matmirnasInMultiplePremirnas")){
    setGeneric("matmirnasInMultiplePremirnas", function(x, ...)
        standardGeneric("matmirnasInMultiplePremirnas"))
}
if(!isGeneric("matmirnasBy")){
    setGeneric("matmirnasBy", function(x, ...)
        standardGeneric("matmirnasBy"))
}
if(!isGeneric("mirbaseVersion"))
    setGeneric("mirbaseVersion", function(x)
        standardGeneric("mirbaseVersion"))
if(!isGeneric("pairData"))
    setGeneric("pairData", function(x, y, object, chooseFunX, chooseFunY,
                                    xNamesAre, yNamesAre, ...)
        standardGeneric("pairData"))
if(!isGeneric("premirnas")){
    setGeneric("premirnas", function(x, ...)
        standardGeneric("premirnas"))
}
if(!isGeneric("premirnasWithMultipleAlignments")){
    setGeneric("premirnasWithMultipleAlignments", function(x, ...)
        standardGeneric("premirnasWithMultipleAlignments"))
}
if(!isGeneric("premirnasBy")){
    setGeneric("premirnasBy", function(x, ...)
        standardGeneric("premirnasBy"))
}
if(!isGeneric("probesets")){
    setGeneric("probesets", function(x, ...)
        standardGeneric("probesets"))
}
if(!isGeneric("probesetsBy")){
    setGeneric("probesetsBy", function(x, ...)
        standardGeneric("probesetsBy"))
}
if(!isGeneric("seqinfo"))
    setGeneric("seqinfo", function(x)
        standardGeneric("seqinfo"))
## if(!isGeneric("transferValues"))
##     setGeneric("transferValues", function(x, object, xNamesAre, toNames, solveFun,
##                                           filter, na.rm, ...)
##         standardGeneric("transferValues"))
if(!isGeneric("transferValues"))
    setGeneric("transferValues", function(x, object, ...)
        standardGeneric("transferValues"))
## if(!isGeneric("organism")){
##     setGeneric("organism", function(x)
##         standardGeneric("organism"))
## }
if(!isGeneric("version")){
    setGeneric("version", function(object, ...)
               standardGeneric("version"))
}

####
## private methods
setGeneric("cleanColumns", function(x, columns, ...)
    starndardGeneric("cleanColumns"))
setGeneric("tablesForColumns", function(x, columns, ...)
    standardGeneric("tablesForColumns"))
##if(!isGeneric("tablesByDegree")){
setGeneric("tablesByDegree", function(x, ...)
    standardGeneric("tablesByDegree"))
##}
