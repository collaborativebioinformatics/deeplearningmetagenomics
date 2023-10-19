###### -- convert a metagenome into a table of domains ------------------------

suppressMessages(library(SynExtend))

# requirements: SynExtend R package -- from bioconductor
# hmmer -- from Sean Eddy's lab -- command line executable
# prodigal -- from dan hyatt, on github somewhere -- command line executable

###### -- arguments -----------------------------------------------------------
# take in fasta file
# a directory with HMMs in a somewhat specific format
# return two RData files that are named after the input file, into the same directory
# that the input file came from
# compressed!

ARGS <- commandArgs(trailingOnly = TRUE)
# example input:
# ARGS <- c("~/Downloads/nmdc/nmdc_mga0wn63_contigs.fna",
#           "~/tmp_models")
INPUT <- ARGS[1L]
MODEL_DIR <- ARGS[2L]
if (length(ARGS) < 3L) {
  LIM <- 10000
}
# will return the example output:
# ~/Downloads/nmdc/nmdc_mga0wn63_contigs.RData

# it is likely that we need to drop contigs smaller than X size
dna <- readDNAStringSet(INPUT)
dna <- unique(dna)
dna <- dna[width(dna) >= LIM]

tmp3 <- tempfile()
writeXStringSet(x = dna,
                filepath = tmp3)

# call prodigal genes
tmp1 <- tempfile()
tmp2 <- paste0(tmp1,
               ".gff")
# we can eventually turn different nobs here if we want
CALLGENES <- paste0("prodigal -p meta -f gff -o ",
                    tmp2,
                    " -i ",
                    tmp3)
tstart <- Sys.time()
system(command = CALLGENES)
tend <- Sys.time()
print(paste("Genes called in",
            tend - tstart,
            "for",
            length(dna),
            "contigs and",
            sum(width(dna)),
            "total nucleotides"))

prodigal_genes <- rtracklayer::import(tmp2)
# munge this into a DFrame that i have a convenient function for handling

z <- prodigal_genes@ranges
z1 <- vector(mode = "list",
             length = length(z))
for (m1 in seq_along(z)) {
  z1[[m1]] <- z[m1]
}
z2 <- IRangesList(z1)
z <- DataFrame("Contig" = as.character(prodigal_genes@seqnames),
               "Range" = z2,
               "ID" = prodigal_genes$ID,
               "Strand" = ifelse(test = as.character(prodigal_genes@strand) == "+",
                                 yes = 0,
                                 no = 1))
prots <- ExtractBy(x = z,
                   y = dna,
                   Verbose = TRUE)
prots <- suppressWarnings(translate(x = prots,
                                    if.fuzzy.codon = "solve"))


target_models <- list.files(path = MODEL_DIR,
                            full.names = TRUE)

modeldir <- paste0(paste0(tempdir(), "/models"))
hmmfile <- paste0(modeldir, "/models.hmm")
dir.create(path = modeldir)

ModFiles <- target_models
PBAR <- length(ModFiles)
pBar <- txtProgressBar(style = 1)
ModAcc <- vector(mode = "list",
                 length = length(ModFiles))
ModTypes <- vector(mode = "character",
                   length = length(ModFiles))
TSTART <- Sys.time()
for (m1 in seq_along(ModFiles)) {
  temp_env <- new.env()
  y <- load(file = ModFiles[m1],
            envir = temp_env,
            verbose = FALSE)
  x <- get(x = y,
           envir = temp_env)
  rm(temp_env)
  
  acc <- x[grepl(pattern = "ACC  ",
                 x = x)]
  acc <- substr(x = acc,
                start = 7,
                stop = nchar(acc))
  ModTypes[m1] <- y
  ModAcc[[m1]] <- acc
  
  cat(paste0(x, "\n"),
      file = hmmfile,
      append = TRUE)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TEND <- Sys.time()
print(TEND - TSTART)

system(command = paste("hmmpress",
                       hmmfile))

temp01 <- tempfile()
temp02 <- tempfile()
writeXStringSet(x = prots,
                filepath = temp01)
# these limits are specified on the table,
# but DomEL can be passed to the initial hmmer command if we so choose
DomLL <- 0.75
DomUL <- 1.25
DomEL <- 1e-3

HMMCOMMAND <- paste("hmmscan",
                    "--cpu 6",
                    "--noali",
                    "--domtblout",
                    temp02,
                    "--domE 1e-3",
                    hmmfile,
                    temp01)

T1 <- Sys.time()
system(command = HMMCOMMAND,
       ignore.stdout = TRUE,
       ignore.stderr = FALSE)
T2 <- Sys.time()
print(T2 - T1)

k2 <- readLines(temp02)
k3 <- k2[-grep(pattern = "^#",
               x = k2)]

k3 <- strsplit(x = k3,
               split = "[ ]+")
k4 <- sapply(k3,
             function(x) {
               c(x[1:22], paste0(x[23:length(x)],
                                 collapse = " "))
             },
             simplify = FALSE,
             USE.NAMES = FALSE)

k4 <- do.call(rbind,
              k4)
DomainTable <- data.frame("target_name" = k4[, 1L],
                          "target_accession" = k4[, 2L],
                          "tlen" = as.integer(k4[, 3L]),
                          "query_name" = k4[, 4L],
                          "query_accession" = k4[, 5L],
                          "qlen" = as.integer(k4[, 6L]),
                          "FULL_E_value" = as.numeric(k4[, 7L]),
                          "FULL_score" = as.numeric(k4[, 8L]),
                          "FULL_bias" = as.numeric(k4[, 9L]),
                          "Percent" = as.integer(k4[, 10L]),
                          "of" = as.integer(k4[, 11L]),
                          "THISDOMAIN_c_Evalue" = as.numeric(k4[, 12]),
                          "THISDOMAIN_i_Evalue" = as.numeric(k4[, 13]),
                          "THISDOMAIN_score" = as.numeric(k4[, 14]),
                          "THISDOMAIN_bias" = as.numeric(k4[, 15]),
                          "HMM_from" = as.integer(k4[, 16]),
                          "HMM_to" = as.integer(k4[, 17]),
                          "ALI_from" = as.integer(k4[, 18]),
                          "ALI_to" = as.integer(k4[, 19]),
                          "ENV_from" = as.integer(k4[, 20]),
                          "ENV_to" = as.integer(k4[, 21]),
                          "acc" = k4[, 22],
                          "desc" = k4[, 23],
                          stringsAsFactors = FALSE)

HMM_Footprint <- DomainTable[, "HMM_to"] - DomainTable[, "HMM_from"] + 1L
ENV_Footprint <- DomainTable[, "ENV_to"] - DomainTable[, "ENV_from"] + 1L

# RetainedDomains <- DomainTable[HMM_Footprint >= DomainTable$tlen * DomLL &
#                                  ENV_Footprint <= (HMM_Footprint * DomUL) &
#                                  ENV_Footprint >= (HMM_Footprint * DomLL) &
#                                  DomainTable$FULL_E_value < DomEL, ]
RetainedDomains <- DomainTable[HMM_Footprint >= DomainTable$tlen * DomLL &
                                 ENV_Footprint <= (HMM_Footprint * DomUL) &
                                 ENV_Footprint >= (HMM_Footprint * DomLL), ]

# save all domains
save(DomtainTable,
     file = gsub(x = INPUT,
                 pattern = "\\.fna$|\\.fasta$|\\.fa$",
                 replacement = "_alldomains.RData"),
     compress = "xz")

save(RetainedDomains,
     file = gsub(x = INPUT,
                 pattern = "\\.fna$|\\.fasta$|\\.fa$",
                 replacement = "_retaineddomains.RData"),
     compress = "xz")