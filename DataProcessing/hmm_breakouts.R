###### -- parse out hmms to specific type groups ------------------------------

suppressMessages(library(SynExtend))

# pfam files are available from ebi
# https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/

# if run from the terminal, expects three trailing arguments
# if run interactively, just supply ARGS in a manner similar to the commented 
# example below
# 1 == the pfam-a.hmm.dat file
# 2 == the pfam-a.hmm file
# 3 == a prefix for the compressed rdata files
ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) != 3L) {
  stop("if run from the terminal this script expects three trailing arguments")
}
# Example:
# ARGS <- c("~/Downloads/Pfam-A.hmm.dat.gz",
#           "~/Downloads/Pfam-A.hmm",
#           "~/Data/20230831_ModelsByType_")

v1 <- readLines(ARGS[1L])
v2 <- readLines(ARGS[2L])

STARTS <- grep(pattern = "# ",
               x = v1)
STARTS2 <- STARTS + 1L
ENDS <- grep(pattern = "//",
             x = v1)
ENDS2 <- ENDS - 1L
TYPE <- grep(pattern = "GF TP",
             x = v1)
ACC <- grep(pattern = "GF AC",
            x = v1)
LINECHAR <- nchar(v1)

pBar <- txtProgressBar(style = 1L)
PBAR <- length(STARTS)
MODELS <- vector(mode = "list",
                 length = PBAR)
M_NAMES <- M_TYPE <- M_ACC <- vector(mode = "character",
                                     length = PBAR)

TSTART <- Sys.time()
for (m1 in seq_along(STARTS2)) {
  MODELS[[m1]] <- c(substr(x = v1[STARTS2[m1]:ENDS2[m1]],
                           start = 11,
                           stop = LINECHAR[STARTS2[m1]:ENDS2[m1]]))
  M_NAMES[m1] <- MODELS[[m1]][1]
  M_TYPE[m1] <- substr(x = v1[TYPE[m1]],
                       start = 11,
                       stop = LINECHAR[TYPE[m1]])
  M_ACC[m1] <- substr(x = v1[ACC[m1]],
                      start = 11,
                      stop = LINECHAR[ACC[m1]])
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TEND <- Sys.time()
print(TEND - TSTART)

# MODELS <- do.call(rbind,
#                   MODELS)

ModelsByType <- tapply(X = M_NAMES,
                       INDEX = M_TYPE,
                       FUN = c)
AccessionByType <- tapply(X = M_ACC,
                          INDEX = M_TYPE,
                          FUN = c)
U_Types <- names(table(M_TYPE))

# this line will need to be edited if HMMER version changes!
HMMSTARTS <- grep(pattern = "HMMER3/f [3.1b2 | February 2015]",
                  x = v2,
                  fixed = TRUE)
HMMENDS <- grep(pattern = "//",
                x = v2)
HMMACC <- grep(pattern = "^ACC",
               x = v2)
HMMLINECHARS <- nchar(v2)
HMMAccessions <- substr(x = v2[HMMACC],
                        start = 7,
                        stop = HMMLINECHARS[HMMACC])

PBAR <- length(AccessionByType)
pBar <- txtProgressBar(style = 1)
TSTART <- Sys.time()
for (m1 in seq_along(AccessionByType)) {
  ph1 <- unlist(mapply(FUN = function(x, y) {
    v2[x:y]
  },
  x = HMMSTARTS[HMMAccessions %in% AccessionByType[[m1]]],
  y = HMMENDS[HMMAccessions %in% AccessionByType[[m1]]]))
  
  assign(x = U_Types[m1],
         value = ph1)
  
  save(list = U_Types[m1],
       file = paste0(ARGS[3L],
                     U_Types[m1],
                     ".RData"),
       compress = "xz")
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TEND <- Sys.time()
print(TEND - TSTART)

# save(list = U_Types,
#      file = ARGS[3L],
#      compress = "xz")


