library(DESeq)
library(ggplot2)
library(geneplotter)
library(EDASeq)
library(RColorBrewer)
library(gplots)

# To use this script the working directory must be set in the same place as the counts file (and use the correct count file name). The path to the
# real cds must be specified too.

Counts <- read.delim("merged_counts_R2T.txt",
                     header = T, row.names = 1)
Counts <- round(Counts)

system("mkdir Results")
system("mkdir Filtered_results")

path.filtered <- "Filtered_results"
path.raw <- "Results"

real.cds <- as.vector(t(read.delim("/media/propietario/ddbf522f-f164-451b-ad73-2967952d697b/SequencesPGC/SizeFactors/R2_SizeFactors2R.txt", header=F)))

PGC11 <- c("PGC11F", "PGC11M")
PGC12 <- c("PGC12F", "PGC12M")
PGC13 <- c("PGC13F", "PGC13M")
SC11 <- c("SC11F", "SC11M")
SC12 <- c("SC12F", "SC12M")
SC13 <- c("SC13F", "SC13M")

conds <- factor(c(PGC11, PGC12, PGC13, SC11, SC12, SC13))
comparisons <- (combn(conds, 2))

filt_5counts <- apply(Counts, 1, function(x) { any(x > 5)})
Counts.ok <- Counts[filt_5counts , ]


cds <- newCountDataSet(Counts.ok, conds)

cds$sizeFactor <- real.cds

cds <- estimateDispersions( cds, method="blind", fitType="local", sharingMode="fit-only" )

write.table(counts(cds,normalized=TRUE), file = paste(path.raw, "/", "NormCounts.txt", sep=""), sep="\t", quote = FALSE)
write.table(summary(counts(cds, normalized = TRUE)), file = paste(path.raw, "/", "NormCountsSummary.txt", sep=""), sep="\t", quote = FALSE)

###################
# Sex comparisons #
###################

# PGC
res <- nbinomTest(cds, "PGC11F", "PGC11M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11FvsPGC11M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11FvsPGC11M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC12F", "PGC12M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC12FvsPGC12M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC12FvsPGC12M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC13F", "PGC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC13FvsPGC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC13FvsPGC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)


# Somatic
res <- nbinomTest(cds, "SC11F", "SC11M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC11FvsSC11M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC11FvsSC11M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC12F", "SC12M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC12FvsSC12M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC12FvsSC12M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC13F", "SC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC13FvsSC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC13FvsSC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)


####################
# Time comparisons #
####################

# Female PGCs
res <- nbinomTest(cds, "PGC11F", "PGC12F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11FvsPGC12F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11FvsPGC12F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC11F", "PGC13F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11FvsPGC13F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11FvsPGC13F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC12F", "PGC13F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC12FvsPGC13F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC12FvsPGC13F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

# Male PGCs
res <- nbinomTest(cds, "PGC11M", "PGC12M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11MvsPGC12M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11MvsPGC12M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC11M", "PGC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11MvsPGC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11MvsPGC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC12M", "PGC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC12MvsPGC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC12MvsPGC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)


# Female Somatics
res <- nbinomTest(cds, "SC11F", "SC12F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC11FvsSC12F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC11FvsSC12F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC11F", "SC13F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC11FvsSC13F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC11FvsSC13F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC12F", "SC13F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC12FvsSC13F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC12FvsSC13F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

# Male Somatics
res <- nbinomTest(cds, "SC11M", "SC12M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC11MvsSC12M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC11MvsSC12M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC11M", "SC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC11MvsSC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC11MvsSC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "SC12M", "SC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "SC12MvsSC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "SC12MvsSC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)


##################
# PGC vs Somatic #
##################
res <- nbinomTest(cds, "PGC11F", "SC11F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11FvsSC11F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11FvsSC11F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC12F", "SC12F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC12FvsSC12F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC12FvsSC12F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC13F", "SC13F")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC13FvsSC13F.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC13FvsSC13F_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC11M", "SC11M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC11MvsSC11M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC11MvsSC11M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC12M", "SC12M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC12MvsSC12M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC12MvsSC12M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

res <- nbinomTest(cds, "PGC13M", "SC13M")
res.05 <- res[(abs(res$log2FoldChange) >= 1) & (res$baseMeanA >=100 | res$baseMeanB >=100), ]
write.table(na.omit(res), file = paste(path.raw, "/", "PGC13MvsSC13M.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
write.table(na.omit(res.05), file = paste(path.filtered, "/", "PGC13MvsSC13M_filt.txt", sep=""), row.names = F, sep="\t", quote = FALSE)
