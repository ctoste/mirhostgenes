## The purpose of this file is to evaluate whether the probe sets identified by the perl function
## do to some extend match the probe sets defined for the genes by the Bioconductor packages.
## read the "new" probe set file.
## read one of the old probe set files.
## test for each miRNA... is there a overlap of probesets.

testProbesets <- function( old, new ){
    Old <- read.table( old, sep="\t", as.is=TRUE, header=TRUE )
    ## Old has a nasty formatting: will just reduce to columns transcript_id and probesets
    New <- read.table( new, sep="\t", as.is=TRUE, header=TRUE )
    ## New has tx_id and feature_id
    ## what i want to have in the end is a list with the names corresponding to the
    ## transcript IDs, the elements being the probe set IDs:
    ## let's start with Old:
    cat("re-formatting old file...")
    Old <- unique( Old[ , c( "transcript_id", "probesets" ) ] )
    Old <- Old[ Old$probesets!="", ]
    Old <- split( Old[ , "probesets" ], f=Old$transcript_id )
    Old <- lapply( Old, function( z ){
        return( unique( unlist( strsplit( z, split=";" ) ) ) )
    } )
    cat( "done.\n" )
    cat("re-formatting new file...")
    New <- New[ New$array_id == "HG-U133_Plus_2", c( "tx_id", "feature_id" )]
    New <- split( New[ , "feature_id" ], f=New$tx_id )
    cat("done.\n")
    ## now we're going to compare...
    only.in.old <- names(Old)[ !(names(Old) %in% names(New)) ]
    only.in.new <- names(New)[ !names(New) %in% names(Old) ]
    in.both <- names(New)[names(New) %in% names(Old)]
    cat("\n SUMMARIES:\n ----------\n")
    cat(" | Number of host transcripts with probesets in both files: ... ", length(in.both), ".\n")
    cat(" | Number of host transcripts with probesets only in old (that's bad): ... ", length(only.in.old), ".\n")
    cat(" | Number of host transcripts with probesets only in new (that's good): ... ", length(only.in.new), ".\n")
    cat(" | Total number of unique probe sets in the old file: ... ",
        length(unique(unlist(Old, use.names=FALSE))), ".\n")
    cat(" | Total number of unique probe sets in the new file: ... ",
        length(unique(unlist(New, use.names=FALSE))), ".\n")
    ## now processing those for which we do have data in both and comparing the probe sets...
    New <- New[in.both]
    Old <- Old[in.both]
    Results <- data.frame(matrix(nrow=length(New), ncol=5), stringsAsFactors=FALSE)
    colnames(Results) <- c("tx_id", "no_new", "no_old", "no_new_in_old", "no_old_in_new")
    Results[ , "tx_id" ] <- names(New)
    for(i in 1:nrow(Results)){
        Results[i, "no_new"] <- length(New[[i]])
        Results[i, "no_old"] <- length(Old[[i]])
        Results[i, "no_new_in_old"] <- sum(New[[i]] %in% Old[[i]])
        Results[i, "no_old_in_new"] <- sum(Old[[i]] %in% New[[i]])
    }
    cat(" | Number of host transcripts for which more probesets were defined in the old file (that's bad): ... ", sum(Results$no_old > Results$no_new), ".\n")
    cat(" | Number of host transcripts for which more probesets are defined in the new file (that's good): ... ", sum(Results$no_new > Results$no_old), ".\n")
    cat(" | Number of host transcripts for which some of the probesets in the new file are NOT in the old file: ... ", sum((Results$no_new - Results$no_new_in_old) > 0 ), ".\n")
    cat(" | Number of host transcripts for which some of the probesets in the old file are NOT in the new file: ... ", sum((Results$no_old - Results$no_old_in_new) > 0 ), ".\n")
    cat(" ----------\n\n")
    list( old=Old, new=New )
}



















