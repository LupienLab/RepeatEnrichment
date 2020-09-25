### Objective : Generate binary matrix

## Implemented frm https://github.com/ManchesterBioinference/Scasat/blob/master/Deconvoluting_cell-types_Scasat_functions.ipynb


args <- commandArgs(trailingOnly = TRUE)

if (length(args)==5) {

  mergedPeakFile <- as.character(args[1])
  peakFolder <- as.character(args[2])
  peakFilePattern <- as.character(args[3])
  outputFolder <- as.character(args[4])
  outputPeakFileName <- as.character(args[5])

} else if (length(args)==0) {
  stop("Please provide the following arguments - Merged peak file, Directory containing all peak files, pattern to identify peak files,
       Output dir and Binary matrix file name", call.=FALSE)
}


######################################
### load dependencies
######################################
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))
suppressMessages(library("genomation"))

#devtools::install_github("walaj/roverlaps")
suppressMessages(library(data.table))
#suppressMessages(library(roverlaps))
#3 Switching to roverlaps for memory efficiency

######################################
### Create function
######################################

files = dir(path=peakFolder,  pattern = peakFilePattern, full.names=TRUE)
query <- suppressMessages(readBed(mergedPeakFile))
queryDF <- data.frame(query)
totalOverlap <- data.frame(seqnames = queryDF$seqnames, start = queryDF$start, end = queryDF$end)
cellName = basename(files)

for (i in 1:length(files)){
  #print(i)
  subject = suppressMessages(readBed(files[i]))
  hits = findOverlaps(query, subject)
  #print("overlapfinding done")
  hitsDF <- data.frame(hits)
  cellName[i] <- gsub(peakFilePattern, '', cellName[i])
  totalOverlap[hitsDF$queryHits, cellName[i]] <- 1
  totalOverlap[-hitsDF$queryHits, cellName[i]] <- 0
  #print("annotation for this repeat done")
}

print("overlap finding completely done, save file")
outputFile = paste0(outputFolder,"/",outputPeakFileName)
saveRDS(totalOverlap, paste0(outputFile, ".rds"))
print("this is it - done")
