# Load required libraries
library(dplyr)
library(stringr)

#read the SRS.gtf filtered from "remove_low_quality_transcripts"
SRS.gtf <- read.delim("/../assembly/SRS.gtf", header=FALSE, quote="'")

#change the gene and transcript prefix from "MSTRG." to "PB."
bff2nocap <- SRS.gtf %>%
  mutate(V2 = str_replace(V2, "PBRI", "PacBio")) %>%
  mutate(V9 = gsub('"MSTRG.', '"PB.', V9))

# write the modified data to a new file
write.table(bff2nocap, "/../NOCAP_MERGE/SRS.modified.gtf", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
