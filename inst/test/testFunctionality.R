## This file lists a large number of possible method calls

detachem <- function(x){
    NS <- loadedNamespaces()
    if(any(NS==x)){
        pkgn <- paste0("package:", x)
        detach(pkgn, unload=TRUE, character.only=TRUE)
    }
}
Pkgs <- c("MirhostDb.Hsapiens.v75.v20", "mirhostgenes")
tmp <- sapply(Pkgs, detachem)
tmp <- sapply(Pkgs, library, character.only=TRUE)

DB <- MirhostDb.Hsapiens.v75.v20
DB
##***************************************
##
##  mature miRNAs
##
##***************************************
cat("Checking mature miRNA methods...")

## simply get all mature miRNAs; the result is however not a unique list of miRNAs,
## since miRNAs from pre-miRNAs with multiple genomic alignments are listed in
## mulitple rows.
MatMir <- matmirnas(DB)
MatMir
length(unique(MatMir$mat_mirna_id))
length(unique(MatMir$pre_mirna_algn_id))

## get all mature miRNAs as GRanges
MatMir <- matmirnas(DB, return.type="GRanges")
MatMir
seqinfo(MatMir)

## get mature miRNAs and add also pre-miRNA names to the table
MatMir <- matmirnas(DB, columns=c(listColumns(DB, "mat_mirna"),
                                  "pre_mirna_name"))
MatMir

## get mat_mirna and pre_mirna entries for mature miRNA MIMAT0000062
MatMir <- matmirnas(DB, columns=unique(c(listColumns(DB, "mat_mirna"),
                             listColumns(DB, "pre_mirna"))),
                    filter=list(MatmirnaidFilter("MIMAT0000062")))
MatMir
## the same mature miRNA is encoded in 3 different pre-miRNAs.

## get all mature miRNAs along with their pre-miRNAs in which they are encoded
MatMir <- matmirnas(DB, columns=c("mat_mirna_id", "mat_mirna_name",
                             "pre_mirna_name", "seq_name"))
MatMir
length(unique(MatMir$mat_mirna_id))
length(unique(MatMir$pre_mirna_name))

## gete all mature miRNAs along with the potential host gene in which they are encoded.
MatMir <- matmirnas(DB, columns=c("mat_mirna_id", "mat_mirna_name",
                             "seq_name", "gene_id", "gene_name", "gene_biotype"))
MatMir
## the mature miRNAs present in host genes.
MatMir.inhg <- MatMir[ !is.na(MatMir$gene_id), ]
MatMir.nohg <- MatMir[ is.na(MatMir$gene_id), ]

MatMir.inhg
## however, a considerable number of "host genes" are actually the pre-miRNAs, which some of them
## are stored in the Ensembl database as "gene" with the biotype "miRNA"
sort( table(MatMir.inhg$gene_biotype), decreasing=TRUE )

## now, get all mature miRNAs for which the gene_biotype!=miRNA
MatMir <- matmirnas(DB, columns=c("mat_mirna_id", "mat_mirna_name",
                             "seq_name", "gene_id", "gene_name", "gene_biotype"),
                    filter=list(GenebiotypeFilter("miRNA", condition="!=")))
MatMir
sum(is.na(MatMir$gene_biotype))
sort( table(MatMir$gene_biotype), decreasing=TRUE )

## get only unique mature miRNAs:
MatMir <- matmirnas(DB, columns=c("mat_mirna_id", "mat_mirna_name"))


## get Mature miRNAs along with probe sets ids.
MatMir <- matmirnas(DB, columns=c("mat_mirna_id", "mat_mirna_name", "probeset_id", "array_id"))
MatMir

##***************************
## matmirnasBy
## get all mature miRNAs grouped by pre-miRNA alignment
matmirnasBy(DB)

## get all mature miRNAs grouped by pre-miRNA
matmirnasBy(DB, by="pre_mirna")

matmirnasBy(DB, by="pre_mirna", use.names=TRUE)

matmirnasBy(DB, by="host_gene")

## get a warnings because some gene names are NA
matmirnasBy(DB, by="host_gene", use.names=TRUE)

MMbyHT <- matmirnasBy(DB, by="host_tx")

MMbyHT.all <- matmirnasBy(DB, by="host_tx", ifEmptyBy="empty")

matmirnasBy(DB, by="host_tx", use.names=TRUE)

Res <- matmirnasBy(DB, by="probeset")

## returning also mat mirnas for which we don't have a probe set
Res <- matmirnasBy(DB, by="probeset", ifEmptyBy="no_probeset")

Res$no_probeset

## get all mature miRNAs as GRanges...
Res <- matmirnasBy(DB, by="pre_mirna", return.type="GRanges")

head(Res)

## get mature miRNAs for pre-miRNA miR-16-1 and miR-16-2
matmirnasBy(DB, filter=list(PremirnaFilter(c("hsa-mir-16-2", "hsa-mir-16-1"))))

matmirnasBy(DB, by="database")
cat("done\n\n")




##***************************************
##
##  pre-miRNAs
##
##***************************************
cat("Checking pre-miRNA methods...")


## get all pre-miRNAs
PreMir <- premirnas(DB)
PreMir
length(unique(PreMir$pre_mirna_name))

## get all pre-miRNAs as GRanges
PreMir <- premirnas(DB, return.type="GRanges")
PreMir

## get all pre-miRNAs along with their miRNA family and their sequence. Since we don't ask for the
## pre_mirna_seq_start and end we get a unique table of pre-miRNAs.
PreMir <- premirnas(DB, columns=c("pre_mirna_name", "mirfam_name", "sequence"))
PreMir

## we have some pre-miRNAs without family
sum(is.na(PreMir$mirfam_name))
## but none without sequence.
sum(is.na(PreMir$sequence))

## get all exonic pre-miRNAs (that are NOT in a host gene of miRNA biotype):
premirnas(DB, columns=c("pre_mirna_id", "pre_mirna_name"),
          filter=list(PositionFilter("exonic"),
              GenebiotypeFilter("miRNA", condition="!=")))

## get all intronic pre-miRNAs (that are NOT in a host gene of miRNA biotype):
premirnas(DB, columns=c("pre_mirna_id", "pre_mirna_name"),
          filter=list(PositionFilter("intronic"),
              GenebiotypeFilter("miRNA", condition="!=")))


##***************************
## premirnasBy
## get the pre-miRNAs by mature mirna
PB <- premirnasBy(DB)
PB

premirnasBy(DB, by="mirfam")

premirnasBy(DB, by="mirfam", use.names=TRUE)

premirnasBy(DB, by="mat_mirna")

premirnasBy(DB, by="mat_mirna", use.names=TRUE)

premirnasBy(DB, by="host_gene")

premirnasBy(DB, by="host_gene", use.names=TRUE)

premirnasBy(DB, by="host_tx")

premirnasBy(DB, by="host_tx", use.names=TRUE)

premirnasBy(DB, by="database")


## add also additional stuff and fetch all pre-miRNAs for host gene SMC4:
premirnasBy(DB, columns=c("pre_mirna_name", "sequence",
                     "mirfam_name", "mat_mirna_name"),
            filter=list(GenenameFilter("SMC4")))

## get all pre-miRNAs by host_gene SMC4
premirnasBy(DB, by="host_gene", filter=list(GenenameFilter("SMC4")))

## get all pre-miRNAs by host_tx of SMC4 and list in which exon/intron the miRNAs are encoded.
premirnasBy(DB, by="host_tx", columns=c("pre_mirna_name",
                                   "in_intron", "in_exon"),
            filter=list(GenenameFilter("SMC4")))

## get all pre-miRNAs by mirfam
premirnasBy(DB, by="mirfam")

## get all pre-miRNAs by mirfam as GRanges
premirnasBy(DB, by="mirfam", return.type="GRanges")

## by database
premirnasBy(DB, by="database")

cat("done\n")


##***************************************
##
##  host transcripts
##
##***************************************
cat("Checking host transcript methods...")

## get all host transcripts from the database.
HT <- hosttx(DB)
HT
nrow(HT)
## the same host_tx might be the host for multiple miRNAs, thus we do have non-unique tx_ids.
length(unique(HT$tx_id))

hosttx(DB, filter=list(TxbiotypeFilter("miRNA")))

## get a unique table of host transcripts
HT <- hosttx(DB, columns=c("tx_id", "tx_biotype", "gene_id"))
HT
nrow(HT)
length(unique(HT$tx_id))

## get the host transcripts along with the corresponding gene.
HT <- hosttx(DB, columns=c("tx_id", "in_intron", "in_exon", "gene_id",
                      "gene_name", "entrezid", "database"))
HT
## in what databases are these transcripts defined?
table(HT$database)
nrow(HT)

## include now also the pre_mirna ids.
HT <- hosttx(DB, columns=c("tx_id", "in_intron", "in_exon", "gene_id",
                      "gene_name", "entrezid", "database", "pre_mirna_id", "pre_mirna_name"))
HT
nrow(HT)
## we have now more rows, since different pre-miRNAs might be associated with the same host_tx
length(unique(HT$tx_id))



##***************************
## hosttxBy
## get the host transcripts by the pre-miRNA
## this will drop automatically empty entries, i.e. pre-miRNAs for which no host transcript
## was defined.
HT <- hosttxBy(DB, by="pre_mirna", columns=c("tx_id", "tx_biotype",
                                        "in_intron", "in_exon", "pre_mirna_name"))
HT

HT <- hosttxBy( DB, by="host_gene", use.names=TRUE )
HT

HT <- hosttxBy( DB, by="host_gene", use.names=TRUE,
               filter=list(GeneidFilter("ENSG%", "like")) )
HT


## to get all of them we scan set drop.empty=FALSE
HT <- hosttxBy(DB, by="pre_mirna", columns=c("tx_id", "tx_biotype", "in_intron",
                                        "in_exon", "pre_mirna_name"), drop.empty=FALSE)
HT

HT <- hosttxBy(DB, by="pre_mirna", columns=c("tx_id", "tx_biotype", "in_intron",
                                       "in_exon", "pre_mirna_name"),
               drop.empty=FALSE,
               use.names=TRUE)
HT


## there are however also some without any entries:
empties <- unlist(lapply(HT, function(z){ return(all(is.na(z$tx_id))) }))
sum(empties)
HT[ empties ]

## host transcripts by gene
HT <- hosttxBy(DB, by="host_gene")
HT

hosttxBy(DB, by="database")

cat("done.\n")


##***************************************
##
##  host genes
##
##***************************************
cat("Checking host gene methods...")

## with the host genes it is just the same as above.
HG <- hostgenes(DB)
HG
length(unique(HG$gene_id))
nrow(HG)
## that's the only unique table...
hostgenes(DB)

hostgenes(DB,columns=c("gene_id", "tx_id"))

hostgenes(DB,columns=c("gene_id", "tx_id", "gene_name"))

hostgenes(DB,columns=c("gene_id", "tx_id"),
          filter=list(DatabaseFilter("core")))



##***************************
## hostgenesBy
## get the host genes by the pre-miRNA
HG <- hostgenesBy(DB, by="pre_mirna")
HG

HG <- hostgenesBy(DB, by="pre_mirna",
                  filter=list(GenebiotypeFilter("miRNA")),
                  use.names=TRUE)
HG

HG <- hostgenesBy(DB, by="pre_mirna", filter=list(GeneidFilter("OTTHUMG00000150446")))
HG

HG <- hostgenesBy(DB, by="pre_mirna", filter=list(DatabaseFilter("core")))
HG

HG <- hostgenesBy(DB, by="pre_mirna", filter=list(DatabaseFilter("otherfeatures")))
HG

HG <- hostgenesBy(DB, by="pre_mirna", filter=list(DatabaseFilter("vega")))
HG

HG <- hostgenesBy(DB, by="pre_mirna", filter=list(DatabaseFilter(c("vega", "otherfeatures"))))
HG

HG <- hostgenesBy(DB, by="pre_mirna",
                  filter=list(DatabaseFilter(c("core", "otherfeatures"))),
                  use.names=TRUE )
HG


HG <- hostgenesBy(DB, by="pre_mirna", use.names=FALSE)
HG

HG <- hostgenesBy(DB, by="pre_mirna", use.names=TRUE)
HG


## get host genes by mirfam
HG <- hostgenesBy(DB, by="mirfam", columns=c("gene_id", "gene_name", "mirfam_name"),
                 use.names=FALSE)
HG

HG <- hostgenesBy(DB, by="mirfam", columns=c("gene_id", "gene_name", "mirfam_name"),
                  use.names=TRUE)
HG


## get host genes by mirfam but only for mirfam mir-15
### CHECK THIS!!!
HG <- hostgenesBy(DB, by="mirfam", columns=c("gene_id", "gene_name", "mirfam_name"),
                  filter=list(MirfamFilter("mir-15")))

mirhostgenes:::.buildQuery( DB, columns=c( "gene_id", "gene_name", "mirfam_name" ) )
## OK
mirhostgenes:::.buildQuery( DB, columns=c( "gene_id", "gene_name", "mirfam_name" ),
                           filter=list(MirfamFilter("mir-15")) )

HG

hostgenesBy(DB, by="database")
cat("done.\n")

##*************************************************
##
## Probe sets stuff...
##
##*************************************************
cat("checking probe set related stuff...")

listArrays(DB)

HT <- hosttx(DB, columns=c(listColumns(DB, "host_tx"), "array_id", "probeset_id"))
nrow(HT)
## have more rows than tx_ids...
length(unique(HT$tx_id))
HT
## we also have some host tx without probe sets
sum(is.na(HT$feature_id))

## that's much slower...
HT <- hosttx(DB, columns=c(listColumns(DB, "host_tx"), "array_id", "probeset_id"),
             filter=list(DatabaseFilter("core")))
nrow(HT)
## have more rows than tx_ids...
length(unique(HT$tx_id))
HT
## we also have some host tx without probe sets
sum(is.na(HT$probeset_id))

HT[HT$pre_mirna_algn_id==1234, ]

## what with pre-miRNAs?
PM <- premirnas(DB, columns=c(listColumns(DB, "pre_mirna"), "array_id", "probeset_id"))
PM

PM.na <- PM[is.na(PM$probeset_id), ]
PM.notna <- PM[!is.na(PM$probeset_id), ]
PM.na
PM.notna

## by="probeset"
listArrays(DB)
HT <- hosttxBy(DB, by="probeset", filter=list(ArrayFilter("HG-U133_Plus_2")))
HT

HG <- hostgenesBy(DB, by="probeset", filter=list(ArrayFilter("HG-U133_Plus_2")))
HG

PM <- premirnasBy(DB, by="probeset", filter=list(ArrayFilter("HG-U133_Plus_2")))
PM

MM <- matmirnasBy(DB, by="probeset", filter=list(ArrayFilter("HG-U133_Plus_2")))
MM

## TODO: add a probesets and probesetsBy.
PS <- probesets(DB)
PS

PS <- probesets(DB, columns=c(listColumns(DB, "array_feature"), "gene_id",
                        "gene_name", "pre_mirna_name"))
PS

## probesetsBy:
PS <- probesetsBy(DB)
PS

PS <- probesetsBy(DB, by="pre_mirna", use.names=TRUE)
PS

PS <- probesetsBy(DB, by="pre_mirna", use.names=TRUE, filter=list(ArrayFilter("PrimeView")))
PS

PS <- probesetsBy(DB, by="mirfam", use.names=TRUE)
PS

PS <- probesetsBy(DB, by="mirfam", use.names=TRUE, drop.empty=FALSE)
PS

cat("done\n")


## testing the nocase stuff...
MF <- MatmirnaFilter("hsa-mir-15b-5p", match.case=TRUE)
MF2 <- MatmirnaFilter("hsa-mir-15b-5p", match.case=FALSE)

matmirnas(DB, filter=list(MF))
matmirnas(DB, filter=list(MF2))

premirnas(DB, filter=list(MF))
premirnas(DB, filter=list(MF2))
