## this file contains all required functions (most being wrappers for external perl scripts)
## to create a hostgenedb for a species.


## this function uses the internal get-mirbase.pl function.
downloadMirbase <- function(version, path=".", force.download=FALSE){
    if(missing(version)){
        version = ""
    }
    if(missing(path)){
        path = "."
    }
    func <- system.file("perl", "get-mirbase.pl", package="mirhostgenes")
    cmd <- paste0("perl ", func, " -p ", path)
    if(!missing(version))
        cmd <- paste0(cmd, " -v ", version)
    if(force.download)
        cmd <- paste0(cmd, " -f")
    ## call the stuff...
    ##Sys.setenv(ENS=ensemblapi)
    system(cmd)
    ##Sys.unsetenv("ENS")
}

.ENSEMBL_DATABASES <- c("core", "cdna", "otherfeatures", "vega")

## this function calls the define_mirna_host_genes.pl script.
## gff: the GFF file (from miRBase) containing the genomic alignments
## of the pre-miRNAs (and mature miRNAs).
## database: character vector with the database name(s) that should be queried.
## host: the host name of the Ensembl database to query. Defaults to the
## main site (ensembldb.ensembl.org).
## user: the user name for the Ensembl database host (default anonymous).
## password: the password for the Ensembl database host.
## ensemblapi: the path to the Ensembl Perl API corresponding to the Ensembl version to query.
## v: be verbose.
defineMirhostgenes <- function(gff,
                               database=c("core"),
                               host="ensembldb.ensembl.org",
                               user="anonymous",
                               pass,
                               ensemblapi,
                               verbose=FALSE){
    if(missing(gff)){
        stop("The GFF file name has to be submitted!")
    }
    ## check databases:
    bm <- database %in% .ENSEMBL_DATABASES
    if(sum(bm)==0)
        stop(paste0("The specified database is not valid! valid entries are ",
                    paste(.ENSEMBL_DATABASES, collapse=", ")))
    database <- database[ bm ]

    func <- system.file("perl", "define_mirna_host_genes.pl", package="mirhostgenes")
    cmd <- paste0("perl ", func,
                  " -g ", gff,
                  " -H ", host,
                  " -U ", user)
    if(!missing(pass)){
        cmd <- paste0(cmd, " -P ", pass)
    }
    cmd <- paste0(cmd, " -d ", paste(database, collapse=",", sep=""))
    if(verbose){
        cmd <- paste0(cmd, " -v")
    }

    if(!missing(ensemblapi)){
        Sys.setenv(ENS=ensemblapi)
    }
    system(cmd)
    if(!missing(ensemblapi)){
        Sys.unsetenv("ENS")
    }
}

## prop.probes: proportion of probes per probe set that have to match the tx.
## max.mm: maximal number of mismatches of the alignment.
## min.probe.algn: minimal required length of the probe alignment; a value
## of 24 means that all nucleotides of a 25nt long probe have to align.
getArrayFeaturesForTx <- function(species,
                                  arrays=c("HG-U133_Plus_2", "PrimeView"),
                                  prop.probes=0.8,
                                  max.mm=0,
                                  min.probe.algn=24,
                                  host="ensembldb.ensembl.org",
                                  user="anonymous",
                                  pass,
                                  ensemblapi,
                                  verbose=FALSE){
    if(missing(species)){
        stop("The species is required!")
    }
    if(prop.probes <= 0 | prop.probes > 1){
        stop("Parameter prop.probes has to be between 0 and 1!\n")
    }
    ## otpions:
    ## -a , separated character string of array names.
    ## -e , ensembl version.
    ## -p, numeric, proportion of probes of a probe set that have to map to the transcript.
    ## -s, species

    func <- system.file("perl", "get_array_features_for_host_tx.pl", package="mirhostgenes")
    cmd <- paste0("perl ", func,
                  " -s ", species,
                  " -p ", prop.probes,
                  " -m ", max.mm,
                  " -M ", min.probe.algn,
                  " -a ", paste0(arrays, collapse=","),
                  " -H ", host,
                  " -U ", user)
    if(!missing(pass)){
        cmd <- paste0(cmd, " -P ", pass)
    }
    if(verbose){
        cmd <- paste0(cmd, " -v")
    }

    if(!missing(ensemblapi)){
        Sys.setenv(ENS=ensemblapi)
    }
    system(cmd)
    if(!missing(ensemblapi)){
        Sys.unsetenv("ENS")
    }
}


## create a SQLite database containing the information defined in the txt files.
makeHostgeneSQLiteFromTables <- function(path="."){
    ## check if we have all files...
    in_files <- c("pre_mirna.txt", "mat_mirna.txt",
                  "host_tx.txt", "host_gene.txt", "metadata.txt",
                  "chromosome.txt")
    add_files <- c("pre_mirna_sequence.txt", "mirfam.txt", "array_features.txt")
    ## check if we have all files...
    all_files <- dir(path, pattern="txt")
    if(sum(in_files %in% all_files)!=length(in_files))
        stop("Something went wrong! I'm missing some of the txt files the perl script should have generated.")
    ## read information
    info <- read.table(paste0(path, "/metadata.txt"), sep="\t", as.is=TRUE, header=TRUE)
    ## define the database name:
    ## it will be e.g. Hsapiens.Ensembl75 for Homo sapiens and Ensembl version 75
    species <- .organismName(info[ info$name=="Organism", "value" ])
    ##substring(species, 1, 1) <- toupper(substring(species, 1, 1))
    ## dbname <- paste0(substring(species, 1, 1),
    ##                  unlist(strsplit(species, split=" ", fixed=TRUE))[ 2 ],
    ##                  ".Ensembl", info[ info$key=="ensembl_version", "value" ],
    ##                  ".mirhostgenes.sqlite")
    dbname <- paste0(.makePackageName(species=info[ info$name=="Organism", "value" ],
                                      ensembl_version=info[ info$name=="ensembl_version", "value" ],
                                      mirbase_version=info[ info$name=="mirbase_version", "value" ]
                                      ), ".sqlite")
    con <- dbConnect(dbDriver("SQLite"), dbname=paste0(path, "/", dbname))

    ## process the mature miRNAs table
    tmp <- read.table(paste0(path, "/mat_mirna.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="mat_mirna", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index mm_pre_mirna_algn_id_idx on mat_mirna (pre_mirna_algn_id);")
    dbGetQuery(con, "create index mm_mat_mirna_name_idx on mat_mirna (mat_mirna_name);")

    ## process pre-miRNA:
    tmp <- read.table(paste0(path, "/pre_mirna.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="pre_mirna", tmp, row.names=FALSE)
    premirnas <- tmp
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index pm_pre_mirna_id_idx on pre_mirna (pre_mirna_id);")
    dbGetQuery(con, "create index pm_pre_mirna_algn_id_idx on pre_mirna (pre_mirna_algn_id);")
    dbGetQuery(con, "create index pm_pre_mirna_name_idx on pre_mirna (pre_mirna_name);")

    ## process host_tx:
    tmp <- read.table(paste0(path, "/host_tx.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="host_tx", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index ht_pre_mirna_algn_id_idx on host_tx (pre_mirna_algn_id);")
    dbGetQuery(con, "create index ht_tx_id_idx on host_tx (tx_id);")
    dbGetQuery(con, "create index ht_gene_id_idx on host_tx (gene_id);")

    ## process host_gene:
    tmp <- read.table(paste0(path, "/host_gene.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="host_gene", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index hg_gene_id_idx on host_gene (gene_id);")
    dbGetQuery(con, "create index hg_database_idx on host_gene (database);")

    ## process chromosome:
    tmp <- read.table(paste0(path, "/chromosome.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="chromosome", tmp, row.names=FALSE)
    rm(tmp)

    ## check additional files.
    add_files <- add_files[ add_files %in% all_files ]
    if(length(add_files) > 0){
        ## do I have the sequence file?
        if(any(add_files=="pre_mirna_sequence.txt")){
            message("Found pre_mirna_sequence.txt, will add support for pre-miRNA sequences.")
            Seq <- read.table(paste0(path, "/pre_mirna_sequence.txt"),
                              sep="\t", as.is=TRUE, header=TRUE)
            ## maybe we should also check that the miRNA ids match the one in pre_mirna.txt...
            if(any(!(Seq$pre_mirna_id %in% premirnas$pre_mirna_id))){
                warning("Ah, the pre_mirna_id from the pre_mirna_sequence.txt and pre_mirna.txt do not match! Drop support for pre-miRNA sequences!")
            }else{
                dbWriteTable(con, name="pre_mirna_sequence", Seq, row.names=FALSE)
                dbGetQuery(con,
                           "create index seq_pre_mirna_id on pre_mirna_sequence (pre_mirna_id);")
            }
        }
        ## do I have a miRNA family definition file?
        if(any(add_files=="mirfam.txt")){
            message("Found mirfam.txt, will add support for pre-miRNA families.")
            MF <- read.table(paste0(path, "/mirfam.txt"), sep="\t", as.is=TRUE, header=TRUE)
            ## maybe we should also check that the miRNA ids match the one in pre_mirna.txt...
            if(any(!(MF$pre_mirna_id %in% premirnas$pre_mirna_id))){
                warning("Ah, the pre_mirna_id from the mirfam.txt and pre_mirna.txt do not match! Drop support for miRNA families!")
            }else{
                dbWriteTable(con, name="mirfam", MF, row.names=FALSE)
                dbGetQuery(con, "create index mf_pre_mirna_id on mirfam (pre_mirna_id);")
            }
        }
        ## do I have a microarray probe set table?
        if(any(add_files=="array_features.txt")){
            message("found array_features.txt, will add support for microarray probe sets.")
            ## first check if ensembl version and species fits...
            firstlines <- readLines("array_features.txt", n=10)
            af_ensembl_version <- unlist(strsplit(
                firstlines[ grep(firstlines, pattern="# ensembl_version", fixed=TRUE) ], split="="))[ 2 ]
            af_species <- unlist(strsplit(
                firstlines[ grep(firstlines, pattern="# species", fixed=TRUE) ], split="="))[ 2 ]
            af_prop_probes <- unlist(strsplit(
                firstlines[ grep(firstlines, pattern="# prop_probes", fixed=TRUE) ], split="="))[ 2 ]
            af_max_mm <- unlist(strsplit(
                firstlines[ grep(firstlines, pattern="# max_mm", fixed=TRUE) ], split="="))[ 2 ]
            af_min_probe_algn <- unlist(
                strsplit(firstlines[ grep(firstlines, pattern="# min_probe_algn", fixed=TRUE) ], split="="))[ 2 ]
            if(af_species!=info[ info$name=="Organism", "value" ] | af_ensembl_version!=info[ info$name=="ensembl_version", "value" ]){
                stop("Found a microarray feature table but the Ensembl version or species does not match! Please either delete file array_features.txt and try again.")
            }
            info <- rbind(info, data.frame(name=c("prop_probes", "max_mm", "min_probe_algn"), value=c(af_prop_probes, af_max_mm, af_min_probe_algn), stringsAsFactors=FALSE))
            ## processing the table
            AF <- read.table(paste0(path, "/array_features.txt"), sep="\t", as.is=TRUE, header=TRUE)
            dbWriteTable(con, name="array_feature", AF, row.names=FALSE)
            rm(AF)
            ## make index
            dbGetQuery(con, "create index af_tx_id_idx on array_feature (tx_id);")
            dbGetQuery(con, "create index af_probeset_id_idx on array_feature (probeset_id);")
        }
    }
    ## write information table
    dbWriteTable(con, name="metadata", info, row.names=FALSE)

    dbDisconnect(con)
    ## done.
    return(dbname)
}


.cleanOrganismName <- function(x){
    return(gsub(.organismName(x), pattern="_", replacement=" ", fixed=TRUE))
}

.organismName <- function(x){
    substring(x, 1, 1) <- toupper(substring(x, 1, 1))
    return(x)
}

.abbrevOrganismName <- function(organism){
    spc <- unlist(strsplit(organism, "_", fixed=TRUE))
    ## this assumes a binomial nomenclature has been maintained.
    return(paste0(substr(spc[[1]], 1, 1), spc[[2]]))
}

## x has to be the connection to the database.
.makePackageName <- function(species, ensembl_version, mirbase_version){
    pkgName <- paste0("MirhostDb.", .abbrevOrganismName(.organismName(species)),
                      ".v", ensembl_version, ".", mirbase_version)
    ##    pkgName <- paste0(.abbrevOrganismName(.organismName(species)),
    ##                      ".Ensembl", ensembl_version, ".mirhostgenes")
    return(pkgName)
}

.makePackageNameFromDb <- function(x){
    species <- .getMetaDataValue(x, "Organism")
    ensembl_version <- .getMetaDataValue(x, "ensembl_version")
    mirbase_version <- .getMetaDataValue(x, "mirbase_version")
    .makePackageName(species=species, ensembl_version=ensembl_version, mirbase_version)
}

## ensdb should be a connection to an SQLite database, or a character string...
makeMirhostgenesPackage <- function(db,
                                    version,
                                    maintainer,
                                    author,
                                    destDir=".",
                                    license="Artistic-2.0"){
    if(class(db)!="character")
        stop("db has to be the name of the SQLite database!")
    dbfile <- db
    db <- MirhostDb(x=dbfile)
    con <- db@con
    pkgName <- .makePackageNameFromDb(con)
    species <- .getMetaDataValue(con, "Organism")
    ensembl_version <- .getMetaDataValue(con, "ensembl_version")
    ## there should only be one template
    template_path <- system.file("pkg-template",package="mirhostgenes")
    ## We need to define some symbols in order to have the
    ## template filled out correctly.
    symvals <- list(
        PKGTITLE=paste0(species, " miRNA host gene definitions"),
        PKGDESCRIPTION=paste0("Contains miRNA host gene definitions for all miRNAs in miRBase ",
            .getMetaDataValue(con, "mirbase_version")),
        PKGVERSION=version,
        AUTHOR=author,
        MAINTAINER=maintainer,
        LIC=license,
        ORGANISM=.organismName(species),
        SPECIES=.organismName(species),
        PROVIDER="miRBase, Ensembl",
        PROVIDERVERSION=paste0("miRBase ",
            .getMetaDataValue(con, "mirbase_version"),
            ", Ensembl ",
            ensembl_version),
        RELEASEDATE= paste0("miRBase release date: ", .getMetaDataValue(con ,'mirbase_date')),
        ORGANISMBIOCVIEW=gsub(" ","_",.organismName(.getMetaDataValue(con ,'Organism'))),
        SOURCEURL="http://www.mirbase.org, http://www.ensembl.org"
        )
    ## Should never happen
    if (any(duplicated(names(symvals)))) {
        str(symvals)
        stop("'symvals' contains duplicated symbols")
    }
    createPackage(pkgname=pkgName,
                  destinationDir=destDir,
                  originDir=template_path,
                  symbolValues=symvals)
    ## then copy the contents of the database into the extdata dir
    sqlfilename <- unlist(strsplit(dbfile, split=.Platform$file.sep))
    sqlfilename <- sqlfilename[ length(sqlfilename) ]
    dir.create(paste(c(destDir, pkgName, "inst", "extdata"),
                     collapse=.Platform$file.sep), showWarnings=FALSE, recursive=TRUE)
    db_path <- file.path(destDir, pkgName, "inst", "extdata",
                         paste(pkgName,"sqlite",sep="."))
    file.copy(dbfile, to=db_path)
}

## x... the hairpin.fa (or hairpin.fa.zip) file
## premirna.file... the pre_mirna.txt file generated by the functions above.
createPremirnaSequenceTable <- function(x, path="."){
    ## process the hairpin.fa file to extract the sequences...
    ## read the
    premirna.file <- paste0(path, .Platform$file.sep, "pre_mirna.txt")
    if(!file.exists(premirna.file))
        stop(premirna.file, " not found!")
    Premirnas <- read.table(premirna.file, sep="\t", as.is=TRUE, header=TRUE)
    if(!any(colnames(Premirnas)=="pre_mirna_id"))
        stop("Required column pre_mirna_id not found in file", premirna.file, "!")
    Premirnas <- unique(Premirnas$pre_mirna_id)
    Seq.table <- matrix(ncol=2, nrow=length(Premirnas))
    colnames(Seq.table) <- c("pre_mirna_id", "sequence")
    Seq.table[ , 1 ] <- Premirnas
    rownames(Seq.table) <- Premirnas
    ## so that's the pre_mirnas for which we want to get the sequence...
    if(length(grep(x, pattern="zip$")) > 0){
        FN <- unlist(strsplit(x, .Platform$file.sep, fixed=TRUE))
        FN <- FN[ length(FN) ]
        FN <- sub(FN, pattern=".zip", replacement="", fixed=TRUE)
        File <- unz(x, filename=FN, open="r")
    }else{
        File <- file(x, open="r")
    }
    on.exit(close(File))
    current_sequence <- ""
    current_mirna <- ""
    addme <- FALSE
    repeat{
        Line <- scan(File, what="character", nlines=1, quiet=TRUE)
        if(length(Line)==0){
            break
        }
        ## that's a header
        if(length(grep(Line, pattern="^>")) > 0){
            if(addme){
                ## add the previous entry...
                Seq.table[ current_mirna, "sequence" ] <- current_sequence
                current_sequence <- ""
                addme <- FALSE
            }
            ## check if we would add this one.
            if(any(Premirnas==Line[ 2 ])){
                addme <- TRUE
                current_mirna <- Line[ 2 ]
                current_sequence <- ""
            }else{
                addme <- FALSE
            }
        }else{
            ## that's the sequence part
            if(addme){
                current_sequence <- paste0(current_sequence, Line[ 1 ])
            }
        }
    }
    ## to be on the save side... add the last one (again)
    Seq.table[ current_mirna, "sequence" ] <- current_sequence
    ##close(File)
    write.table(Seq.table, file=paste0(path, .Platform$file.sep, "pre_mirna_sequence.txt")
              , sep="\t", row.names=FALSE)
}

## x... the miFam.dat file from miRBase
createMirfamTable <- function(x, path="."){
    ## scan the file line by line
    ## AC field: accession/mirfam_id
    ## ID field: mirfam_name
    ## MI entries: pre_mirna_id pre_mirna_name
    ## // end of entry; also very last line in the file.
    ## means: we record the AC and ID, each time we find a MI that is also in the
    ## pre_mirna.txt we write the line to the output file.
    ## check input files:
    premirna.file <- paste0(path, .Platform$file.sep, "pre_mirna.txt")
    if(!file.exists(premirna.file))
        stop(premirna.file, " not found!")
    Premirnas <- read.table(premirna.file, sep="\t", as.is=TRUE, header=TRUE)
    if(!any(colnames(Premirnas)=="pre_mirna_id"))
        stop("Required column pre_mirna_id not found in file", premirna.file, "!")
    Premirnas <- unique(Premirnas$pre_mirna_id)
    ## check the miFam file: is it a zip file?
    if(length(grep(x, pattern="zip$")) > 0){
        FN <- unlist(strsplit(x, .Platform$file.sep, fixed=TRUE))
        FN <- FN[ length(FN) ]
        FN <- sub(FN, pattern=".zip", replacement="", fixed=TRUE)
        miFam <- unz(x, filename=FN, open="r")
    }else{
        miFam <- file(x, open="r")
    }
    on.exit(close(miFam), add=TRUE)
    ## open the output file:
    outFile <- file(paste0(path, .Platform$file.sep, "mirfam.txt"), open="w")
    cat("mirfam_id\tmirfam_name\tpre_mirna_id\n", file=outFile, append=FALSE)
    on.exit(if(isOpen(outFile)){ close(outFile) }, add=TRUE)
    current_mirfam_id <- ""
    current_mirfam_name <- ""
    repeat{
        Line <- scan(miFam, what="character", nlines=1, quiet=TRUE)
        if(length(Line)==0){
            break
        }
        if(Line[ 1 ]=="AC")
            current_mirfam_id <- Line[ 2 ]
        if(Line[ 1 ]=="ID")
            current_mirfam_name <- Line[ 2 ]
        if(Line[ 1 ]=="MI"){
            ## check if we have this miRNA:
            if(any(Premirnas==Line[ 2 ])){
                cat(paste0(current_mirfam_id, "\t",
                           current_mirfam_name, "\t", Line[ 2 ], "\n"),
                    file=outFile, append=TRUE)
            }
        }
    }
    ## if(isOpen(miFam))
    ##     close(miFam)
    ## if(isOpen(outFile))
    ##     close(outFile)
}

## basically, we're reading additional stuff from the miRBase and storing them into
## the (existing) database tables.
## read miRNA_high_conf.dat.gz
## read high_conf_hairpin.fa.gz to get high confidence pre-miRNAs.
## read high_conf_mature.fa.gz to get high confidence mature miRNAs.
## read database_files/mature_read_count.txt.gz: auto_mature, mature_acc, read_count, experiment_count
## read database_files/mirna_read_count.txt.gz: auto_mirna, mirna_acc, read_count, experiment_count
## mirbase.path: path to the local miRBase
## path: path to the directory where the mirhostgenes tables can be found.
fetchAdditionalInformation <- function(mirbase.path=".", path=".", verbose=FALSE){
    ## read the pre-miRNA table.
    did.pre.conf <- FALSE
    did.pre.read <- FALSE
    did.pre.exp <- FALSE
    did.mat.conf <- FALSE
    did.mat.read <- FALSE
    did.mat.exp <- FALSE
    premirfile <- paste0(path, .Platform$file.sep, "pre_mirna.txt")
    if(!file.exists(premirfile))
        stop(paste0("Can not find 'pre_mirna.txt' in ", path,
                    "! Please set 'path' to point to the directory where the files generated by defineMirhostgenes are located!"))
    ## message(paste0("Reading pre-miRNA definitions from: ", premirfile, "..."), appendLF=FALSE)
    premirs <- read.table(premirfile, sep="\t", as.is=TRUE, header=TRUE)
    ## message("OK")
    ## first check if we can get the high_conf fa files:
    inFile <- paste0(mirbase.path, .Platform$file.sep, "high_conf_hairpin.fa.gz")
    if(file.exists(inFile)){
        highConfAcc <- .parseHighConfFA(inFile)
    }else{
        inFile <- paste0(mirbase.path, .Platform$file.sep, "high_conf_hairpin.fa")
        if(file.exists(inFile)){
            highConfAcc <- .parseHighConfFA(inFile)
        }else{
            highConfAcc <- NULL
        }
    }
    if(is.null(highConfAcc)){
        ## check if we have the high_conf.dat file.
        inFile <- paste0(mirbase.path, .Platform$file.sep, "miRNA_high_conf.dat.gz")
        if(file.exists(inFile)){
            highConfAcc <- .parseHighConfDat(inFile)
        }else{
            ## was it eventually unzipped?
            inFile <- paste0(mirbase.path, .Platform$file.sep, "miRNA_high_conf.dat")
            if(file.exists(inFile)){
                highConfAcc <- .parseHighConfDat(inFile)
            }else{
                highConfAcc <- NULL
            }
        }
    }
    if(!is.null(highConfAcc)){
        if(verbose)
            message("Adding pre-miRNA confidence data...", appendLF=FALSE)
        ## adding confidence data to pre-miRNAs.
        did.pre.conf <- TRUE
        if(!any(colnames(premirs) == "pre_mirna_confidence"))
            premirs <- cbind(premirs, pre_mirna_confidence=rep(0, nrow(premirs)),
                             stringsAsFactors=FALSE)
        premirs[premirs$pre_mirna_id %in% highConfAcc, "pre_mirna_confidence"] <- 1
        if(verbose)
            message("OK")
    }else{
        warning("Did not find required files to read confidence data. High confidence information will thus not be available!")
    }
    ## processing pre-miRNA read count:
    inFile <- paste0(mirbase.path, .Platform$file.sep, "database_files",
                     .Platform$file.sep, "mirna_read_count.txt.gz")
    if(file.exists(inFile)){
        readData <- .parseReadCountData(inFile)
    }else{
        inFile <- paste0(mirbase.path, .Platform$file.sep, "database_files",
                         .Platform$file.sep, "mirna_read_count.txt")
        if(file.exists(inFile)){
            readData <- .parseReadCountData(inFile)
        }else{
            readData <- NULL
        }
    }
    ## adding read count data.
    if(!is.null(readData)){
        if(verbose)
            message("Adding read count data for pre-miRNAs...", appendLF=FALSE)
        did.pre.read <- TRUE
        did.pre.exp <- TRUE
        if(!any(colnames(premirs) == "pre_mirna_read_count"))
            premirs <- cbind(premirs, pre_mirna_read_count=rep(0, nrow(premirs)),
                             stringsAsFactors=FALSE)
        if(!any(colnames(premirs) == "pre_mirna_experiment_count"))
            premirs <- cbind(premirs, pre_mirna_experiment_count=rep(0, nrow(premirs)),
                             stringsAsFactors=FALSE)
        readData <- readData[readData$acc %in% premirs$pre_mirna_id, ]
        rownames(readData) <- readData$acc
        ## now add the data:
        premirsHaveData <- premirs$pre_mirna_id %in% readData$acc
        premirs[premirsHaveData, "pre_mirna_read_count"] <-
            readData[premirs$pre_mirna_id[premirsHaveData], "read_count"]
        premirs[premirsHaveData, "pre_mirna_experiment_count"] <-
            readData[premirs$pre_mirna_id[premirsHaveData], "exp_count"]
        rm(readData)
        if(verbose)
            message("OK")
    }
    ## save the premirnas again
    write.table(premirs, file=premirfile, row.names=FALSE, quote=FALSE, sep="\t")
    ##
    ## now processing mature miRNAs.
    matmirfile <- paste0(path, .Platform$file.sep, "mat_mirna.txt")
    if(!file.exists(matmirfile))
        stop(paste0("Can not find 'mat_mirna.txt' in ", path,
                    "! Please set 'path' to point to the directory where the files generated by defineMirhostgenes are located!"))
    ## message(paste0("Reading pre-miRNA definitions from: ", premirfile, "..."), appendLF=FALSE)
    matmirs <- read.table(matmirfile, sep="\t", as.is=TRUE, header=TRUE)
    ## message("OK")
    ## first check if we can get the high_conf fa files:
    inFile <- paste0(mirbase.path, .Platform$file.sep, "high_conf_mature.fa.gz")
    if(file.exists(inFile)){
        highConfAcc <- .parseHighConfFA(inFile)
    }else{
        inFile <- paste0(mirbase.path, .Platform$file.sep, "high_conf_mature.fa")
        if(file.exists(inFile)){
            highConfAcc <- .parseHighConfFA(inFile)
        }else{
            highConfAcc <- NULL
        }
    }
    if(is.null(highConfAcc)){
        ## check if we have the high_conf.dat file.
        inFile <- paste0(mirbase.path, .Platform$file.sep, "miRNA_high_conf.dat.gz")
        if(file.exists(inFile)){
            highConfAcc <- .parseHighConfDat(inFile)
        }else{
            ## was it eventually unzipped?
            inFile <- paste0(mirbase.path, .Platform$file.sep, "miRNA_high_conf.dat")
            if(file.exists(inFile)){
                highConfAcc <- .parseHighConfDat(inFile)
            }else{
                highConfAcc <- NULL
            }
        }
    }
    if(!is.null(highConfAcc)){
        if(verbose)
            message("Adding mature miRNA confidence data...", appendLF=FALSE)
        did.mat.conf <- TRUE
        ## adding confidence data to pre-miRNAs.
        if(!any(colnames(matmirs) == "mat_mirna_confidence"))
            matmirs <- cbind(matmirs, mat_mirna_confidence=rep(0, nrow(matmirs)),
                             stringsAsFactors=FALSE)
        matmirs[matmirs$mat_mirna_id %in% highConfAcc, "mat_mirna_confidence"] <- 1
        if(verbose)
            message("OK")
    }
    readData <- NULL
    ## processing mat-miRNA read count:
    inFile <- paste0(mirbase.path, .Platform$file.sep, "database_files",
                     .Platform$file.sep, "mature_read_count.txt.gz")
    if(file.exists(inFile)){
        readData <- .parseReadCountData(inFile)
    }else{
        inFile <- paste0(mirbase.path, .Platform$file.sep, "database_files",
                         .Platform$file.sep, "mature_read_count.txt")
        if(file.exists(inFile)){
            readData <- .parseReadCountData(inFile)
        }else{
            readData <- NULL
        }
    }
    ## adding read count data.
    if(!is.null(readData)){
        if(verbose)
            message("Adding read count data for mat-miRNAs...", appendLF=FALSE)
        did.mat.read <- TRUE
        did.mat.exp <- TRUE
        if(!any(colnames(matmirs) == "mat_mirna_read_count"))
            matmirs <- cbind(matmirs, mat_mirna_read_count=rep(0, nrow(matmirs)),
                             stringsAsFactors=FALSE)
        if(!any(colnames(matmirs) == "mat_mirna_experiment_count"))
            matmirs <- cbind(matmirs, mat_mirna_experiment_count=rep(0, nrow(matmirs)),
                             stringsAsFactors=FALSE)
        readData <- readData[readData$acc %in% matmirs$mat_mirna_id, ]
        rownames(readData) <- readData$acc
        ## now add the data:
        matmirsHaveData <- matmirs$mat_mirna_id %in% readData$acc
        matmirs[matmirsHaveData, "mat_mirna_read_count"] <-
            readData[matmirs$mat_mirna_id[matmirsHaveData], "read_count"]
        matmirs[matmirsHaveData, "mat_mirna_experiment_count"] <-
            readData[matmirs$mat_mirna_id[matmirsHaveData], "exp_count"]
        rm(readData)
        if(verbose)
            message("OK")
    }
    ## save the matmirnas again
    write.table(matmirs, file=matmirfile, row.names=FALSE, quote=FALSE, sep="\t")
    ## do the rest...
    miFam <- paste0(mirbase.path, .Platform$file.sep, "miFam.dat.gz")
    if(file.exists(miFam)){
        if(verbose)
            message("Getting miRNA family definitions...", appendLF=FALSE)
        createMirfamTable(miFam, path=path)
        if(verbose)
            message("OK")
    }else{
        miFam <- paste0(mirbase.path, .Platform$file.sep, "miFam.dat")
        if(file.exists(miFam)){
        if(verbose)
            message("Getting miRNA family definitions...", appendLF=FALSE)
            createMirfamTable(miFam, path=path)
            if(verbose)
                message("OK")
        }else{
            warning(paste0("Did not find file 'miFam.dat.gz' in folder ",
                           mirbase.path,
                           ". Skipped generating the miRFam database table."))
        }
    }
    seqFile <- paste0(mirbase.path, .Platform$file.sep, "hairpin.fa.gz")
    if(file.exists(seqFile)){
        if(verbose)
            message("Getting pre-miRNA sequence data...", appendLF=FALSE)
        createPremirnaSequenceTable(seqFile, path=path)
        if(verbose)
            message("OK")
    }else{
        seqFile <- paste0(mirbase.path, .Platform$file.sep, "hairpin.fa")
        if(file.exists(seqFile)){
            if(verbose)
                message("Getting pre-miRNA sequence data...", appendLF=FALSE)
            createPremirnaSequenceTable(seqFile, path=path)
            if(verbose)
                message("OK")
        }else{
            warning(paste0("Did not find file 'hairpin.fa.gz' in folder ",
                           mirbase.path,
                           ". Skipped generating the miRNA sequence table."))
        }
    }
    ## adding info to the metadata.
    metaFile <- paste0(path, .Platform$file.sep, "metadata.txt")
    if(file.exists(metaFile)){
        meta <- read.table(metaFile, sep="\t", as.is=TRUE, header=TRUE)
        if(did.pre.conf)
            meta <- rbind(meta, c("have_premirna_confidence", "TRUE"))
        if(did.pre.read)
            meta <- rbind(meta, c("have_premirna_readcount", "TRUE"))
        if(did.mat.conf)
            meta <- rbind(meta, c("have_matmirna_confidence", "TRUE"))
        if(did.mat.read)
            meta <- rbind(meta, c("have_matmirna_readcount", "TRUE"))
        write.table(meta, file=metaFile, sep="\t", row.names=FALSE, quote=FALSE)
    }
}

## returns a character vector of accession numbers
.parseHighConfDat <- function(x){
    ##message(paste0("Reading miRNA confidence data from ", x, "..."), appendLF=FALSE)
    confData <- readLines(x)
    Accs <- confData[grep(confData, pattern="^AC")]
    Accs <- gsub(Accs, pattern="AC\\s+", replacement="")
    Accs <- gsub(Accs, pattern=";", replacement="")
    ##message("OK\n")
    return(Accs)
}

.parseReadCountData <- function(x){
    Data <- read.table(x, sep="\t", as.is=TRUE)
    colnames(Data)[1:4] <- c("id", "acc", "read_count", "exp_count")
    return(Data)
}

## read and parse accession numbers from FA file
.parseHighConfFA <- function(x){
    confData <- readLines(x)
    confData <- confData[grep(confData, pattern="^>")]
    tmp <- sapply(confData, function(z){
        return(unlist(strsplit(z, split=" ", fixed=TRUE), use.names=FALSE)[2])
    }, USE.NAMES=FALSE)
    return(tmp)
}


