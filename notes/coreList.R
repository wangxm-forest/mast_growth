############Started by Mao##############
############April 22, 2026##############
############Extract file names##########

setwd("C:/PhD/scanning/MORA_2024/TO04")
file <- list.files(pattern=NULL, all.files=FALSE,
                   full.names=FALSE)
file <- do.call(rbind, strsplit(file, "_"))

colnames(file) <- c("Species", "Stand", "Tag")

df <- as.data.frame(file)
df
write.csv(df, "C:/PhD/Project/PhD_thesis/mast_growth/notes/TO04.csv")
