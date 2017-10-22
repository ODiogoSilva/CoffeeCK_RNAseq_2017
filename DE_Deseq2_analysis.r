
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("reshape")
library("plotly")
library("grid")
library("gridExtra")
library("plyr")
library("VennDiagram")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

c <- read.csv("gene_count_matrix_sorted.csv", row.names=1)
gcts <- as.matrix(c)
head(gcts)

c <- read.csv("transcript_count_matrix.csv", row.names=1)
tcts <- as.matrix(c)
head(tcts)

colcols <- c("exp", "cond", "hpi", "rep_bio", "rep_tec")

coldata <- read.csv("pheno_data.csv", row.names=1, strip.white=T)
coldata <- coldata[, colcols]

# Correct the technical replicate column for later collapsing
coldata$rep_tec <- c(rep(1:4, each=4), c(5,5,5), rep(6:24, each=4))
head(coldata)

# Check with the gene counts
all(rownames(coldata) %in% colnames(gcts))
# Check with the transcript counts
all(rownames(coldata) %in% colnames(tcts))

gdds <- DESeqDataSetFromMatrix(countData = gcts,
                               colData = coldata,
                               design = ~ cond)
tdds <- DESeqDataSetFromMatrix(countData = tcts,
                               colData = coldata,
                               design = ~ cond)
summary(gdds)
summary(tdds)

annot <- read.csv("annotations.txt", header=F)
colnames(annot) <- c("gene_id", "chr", "g_star", "g_end")
# Re-order annotation gene_ids according to the order in the data set
annot <- annot[match(rownames(gdds), annot$gene_id),]
head(annot)

gdds <- collapseReplicates(gdds, gdds$rep_tec)

mcols(gdds) <- cbind(mcols(gdds), annot)

ntd <- normTransform(gdds)
plotPCA(ntd, intgroup=c("exp", "cond"))

### subset conditions ###
pair_cond = list(
    # Resistant at 24HPI
    "rh24c" = colData(gdds)$hpi=="24h" & colData(gdds)$exp=="R",
    # Resistant at 48HPI
    "rh48c" = colData(gdds)$hpi=="48h" & colData(gdds)$exp=="R",
    # Resistant at 72HPI
    "rh72c" = colData(gdds)$hpi=="72h" & colData(gdds)$exp=="R",
    # Susceptible at 24HPI
    "sh24c" = colData(gdds)$hpi=="24h" & colData(gdds)$exp=="S",
    # Susceptible at 48HPI
    "sh48c" = colData(gdds)$hpi=="48h" & colData(gdds)$exp=="S",
    # Susceptible at 72HPI
    "sh72c" = colData(gdds)$hpi=="72h" & colData(gdds)$exp=="S",
    ######
    # Control at 24HPI
    "c24" = colData(gdds)$cond=="C" & colData(gdds)$hpi=="24h",
    # Control at 48HPI
    "c48" = colData(gdds)$cond=="C" & colData(gdds)$hpi=="48h",
    # Control at 72HPI
    "c72" = colData(gdds)$cond=="C" & colData(gdds)$hpi=="72h"
)

getpair_stats <- function(ds, cond) {
    pc <- unlist(pair_cond[cond])
    sdd <- subset(ds, select=pc)
    
    if (cond %in% c("c24", "c48", "c72")) {
        design(sdd) <- formula(~ exp)
    }
    
    sddd <- DESeq(sdd);
    sddres <- results(sddd, alpha=0.05, lfcThreshold=1, altHypothesis="greaterAbs")
    
    #print(sdd@colData)
    print(paste("Data set with ", nrow(sdd@colData), " samples and ", length(sdd@rowRanges), " genes."))
    #print(sddres)
    print(paste("Total of significant DE :", sum(sddres$padj <= 0.05, na.rm=T)))
    
    return(list("res"=sddres, "ds"=sdd))
}

rh24 <- getpair_stats(gdds, "rh24c")

rh48 <- getpair_stats(gdds, "rh48c")

rh72 <- getpair_stats(gdds, "rh72c")

sh24 <- getpair_stats(gdds, "sh24c")

sh48 <- getpair_stats(gdds, "sh48c")

sh72 <- getpair_stats(gdds, "sh72c")

pairsr <- list(
    "rh24" = rh24,
    "rh48" = rh48,
    "rh72" = rh72
)

pairss <- list(
    "sh24" = sh24,
    "sh48" = sh48,
    "sh72" = sh72
)

x <- data.frame(
    "x"=c("h24", "h48", "h72"),
    "r"=unlist(lapply(pairsr, function(x) sum(x$res$padj <= 0.05, na.rm=T))),
    "s"=unlist(lapply(pairss, function(x) sum(x$res$padj <= 0.05, na.rm=T)))
)
x.m <- melt(x)

p <- ggplot(x.m, aes(x, value)) +
    geom_bar(aes(fill=variable), stat = "identity", position="dodge") +
    xlab("HPI")+ ylab("Frequency") + ggtitle ("Number of DE genes for each phenotype across HPI") +
    scale_fill_discrete(name = "Phenotype")

ggplotly(p)
x

x <- data.frame(
    "x"=c("h24", "h48", "h72"),
    "rUP"=unlist(lapply(pairsr, function(x) sum(subset(x$res, x$res$padj <= 0.05)$log2FoldChange >= 1, na.rm=T))),
    "rDOWN"=unlist(lapply(pairsr, function(x) sum(subset(x$res, x$res$padj <= 0.05)$log2FoldChange <= -1, na.rm=T))),
    "sUP"=unlist(lapply(pairss, function(x) sum(subset(x$res, x$res$padj <= 0.05)$log2FoldChange >= 1, na.rm=T))),
    "sDOWN"=unlist(lapply(pairss, function(x) sum(subset(x$res, x$res$padj <= 0.05)$log2FoldChange <= -1, na.rm=T)))
)
x.m <- melt(x)

p <- ggplot(x.m, aes(variable, value,fill=variable)) +
    geom_bar(aes(), stat = "identity", position="dodge") +
    facet_wrap(~x, ncol=3) +
    xlab("")+ ylab("Frequency") + ggtitle("Number of UP and DOWN regulated DE genes for each phenotype across HPI") +
    scale_fill_discrete(name = "Phenotype")

ggplotly(p)
#x.m

# Subset the results of each time point and get the rownames, which contain the number of the gene.
rh24_gi <- rownames(subset(rh24$res, rh24$res$padj <= 0.05))
rh48_gi <- rownames(subset(rh48$res, rh48$res$padj <= 0.05))
rh72_gi <- rownames(subset(rh72$res, rh72$res$padj <= 0.05))

ven_rh <- list(
    rh24=rh24_gi,
    rh48=rh48_gi,
    rh72=rh72_gi
)

rh_ven <- venn.diagram(ven_rh, main="Resistant DE over time", filename=NULL,print.mode=c("raw", "percent"))

# Subset the results of each time point and get the rownames, which contain the number of the gene.
sh24_gi <- rownames(subset(sh24$res, sh24$res$padj <= 0.05))
sh48_gi <- rownames(subset(sh48$res, sh48$res$padj <= 0.05))
sh72_gi <- rownames(subset(sh72$res, sh72$res$padj <= 0.05))

ven_sh <- list(
    sh24=sh24_gi,
    sh48=sh48_gi,
    sh72=sh72_gi
)

sh_ven <- venn.diagram(ven_sh, main="Susceptible DE over time", filename=NULL, print.mode=c("raw", "percent"))

library(repr)
options(repr.plot.height=4)
grid.arrange(gTree(children=rh_ven), gTree(children=sh_ven), ncol=2, nrow=1, heights=c(1), widths=c(1,1))

# Subset the results of each time point and get the rownames, which contain the number of the gene.
rh24_gi <- rownames(subset(rh24$res, rh24$res$padj <= 0.05))
sh24_gi <- rownames(subset(sh24$res, sh24$res$padj <= 0.05))

ven_24h <- list(
    rh24=rh24_gi,
    sh24=sh24_gi
)

h24_ven <- venn.diagram(ven_24h, main="24HPI overlap between phenotypes", filename=NULL, print.mode=c("raw", "percent"))

# Subset the results of each time point and get the rownames, which contain the number of the gene.
rh48_gi <- rownames(subset(rh48$res, rh48$res$padj <= 0.05))
sh48_gi <- rownames(subset(sh48$res, sh48$res$padj <= 0.05))

ven_48h <- list(
    rh48=rh48_gi,
    sh48=sh48_gi
)

h48_ven <- venn.diagram(ven_48h, main="24HPI overlap between phenotypes", filename=NULL, print.mode=c("raw", "percent"))

# Subset the results of each time point and get the rownames, which contain the number of the gene.
rh72_gi <- rownames(subset(rh72$res, rh72$res$padj <= 0.05))
sh72_gi <- rownames(subset(sh72$res, sh72$res$padj <= 0.05))

ven_72h <- list(
    rh72=rh72_gi,
    sh72=sh72_gi
)

h72_ven <- venn.diagram(ven_72h, main="72HPI overlap between phenotypes", filename=NULL, imagetype="png",print.mode=c("raw", "percent"))

library(repr)
options(repr.plot.height=6)
blank<-rectGrob(gp=gpar(col="white"))
grid.arrange(gTree(children=h24_ven), blank, gTree(children=h48_ven), gTree(children=h72_ven), ncol=3, nrow=2, widths=c(0.8, 0.1, 0.8))

multiMA <- function(x) {
    plotMA(x$res)
}

par(mfrow = c(2,3), mar=c(7, 4, 3, 0))


y <- lapply(pairsr, multiMA)
mtext("Resistant", side = 3, line = -2, outer = TRUE)
mtext("Susceptible", side = 3, line = -28, outer = TRUE)
mtext("24HPI", at = 0.18, side = 3, line = -52, outer = TRUE)
mtext("48HPI", at = 0.53, side = 3, line = -52, outer = TRUE)
mtext("72HPI", at = 0.85, side = 3, line = -52, outer = TRUE)
y <- lapply(pairss, multiMA)

plotFC <- function(ds, title) {
    ds$chr <- mcols(gdds)$chr;
    mdf <- data.frame("gene" = rownames(ds),
                      "chr" = mcols(gdds)$chr ,
                      "fc" = ds$log2FoldChange,
                      "qval" = ds$padj)
    df <- subset(mdf, mdf$qval <= 0.05)
    dfc <- count(df, "chr")
    p <- ggplot(df, aes(as.numeric(rownames(df)), fc)) +
        geom_point(aes(color=fc)) +
        geom_hline(data = dfc,
                   aes(yintercept=freq * (max(df$fc) / max(dfc$freq)), alpha=.7, text=paste("Frequency: ", freq)),
                   colour="red") + 
        facet_grid(~chr, scales="free_x") +
        scale_colour_gradientn(colours=rainbow(7)) +
        xlab("Genome position") +
        ggtitle(paste(title, "(", sum(dfc$freq) ,")")) + 
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.spacing=unit(0.15, "lines"))
        #scale_y_continuous(name = "log2(FC)",
        #                   sec.axis = sec_axis(~ . * (max(df$fc) / max(dfc$freq)) , name = "Frequency"),
        #                   limits = c(min(df$fc), max(df$fc)))
    
    ggplotly(p)
    
}

plotFC(rh24$res, "Resistant 24HPI")
plotFC(rh48$res, "Resistant 48HPI")
plotFC(rh72$res, "Resistant 72HPI")

plotFC(sh24$res, "Susceptible 24HPI")
plotFC(sh48$res, "Susceptible 48HPI")
plotFC(sh72$res, "Susceptible 72HPI")

get_tcstats <- function(ds) {
    design(ds) <- formula(~ cond + hpi + cond:hpi);
    dsd <- DESeq(ds, test="LRT", reduced=~ cond + hpi);
    res <- results(dsd);
    return(list("res"=res, "ds"=dsd))
}

get_tcstats_cond <- function(ds) {
    design(ds) <- formula(~ exp + hpi + exp:hpi);
    dsd <- DESeq(ds, test="LRT", reduced=~ exp + hpi);
    res <- results(dsd, alpha=0.05);
    return(list("res"=res, "ds"=dsd))
}

# Subset the data to include only Q samples
rqt = colData(gdds)$exp=="R"
rds <- subset(gdds, select=rqt)

rtc <- get_tcstats(rds)

length(which(rtc$res$padj <= 0.05))

betas <- coef(rtc$ds)
topGenes <- head(order(rtc$res$padj), 83)

mat <- betas[topGenes, -c(1,2)]

thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=F)

options(warn=-1)
sprof <- function(gnID) {
    fiss <- plotCounts(rtc$ds, gnID, 
                       intgroup = c("hpi","cond"), returnData = TRUE)
    p <- ggplot(fiss,
      aes(x = as.numeric(hpi), y = count, color = cond, group = cond)) +
      geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()
    return(p)
}

gene27 <- head(order(rtc$res$padj), 27)

z <- lapply(gene27, sprof)

options(repr.plot.height=16)
multiplot(plotlist=z, cols=3)

# Subset the data to include only Q samples
sqt = colData(gdds)$exp=="S"
sds <- subset(gdds, select=sqt)

stc <- get_tcstats(sds)

length(which(stc$res$padj <= 0.05))

sbetas <- coef(stc$ds)
stopGenes <- head(order(stc$res$padj), 25)

smat <- betas[stopGenes, -c(1,2)]

thr <- 3
smat[smat < -thr] <- -thr
smat[smat > thr] <- thr

options(repr.plot.height=6)
pheatmap(smat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=F)

rde <- which(rtc$res$padj <= 0.05)
sde <- which(stc$res$padj <= 0.05)

tc_ds <- list(
    "rv" = rde,
    "sv" = sde
)
options(repr.plot.height=5)
tc_ven <- venn.diagram(tc_ds, filename=NULL,  imagetype="png",print.mode=c("raw", "percent"))
grid.draw(tc_ven)

plotFC(rtc$res, "DEG over time for Resistant")
plotFC(stc$res, "DEG over time for Susceptible")

hpi24_c <- colData(gdds)$hpi=="24h"
hpi48_c <- colData(gdds)$hpi=="48h"
hpi72_c <- colData(gdds)$hpi=="72h"

hpi_intest <- function(sel) {
    ss <- subset(gdds, select=sel)
    design(ss) <- formula(~ cond + exp + cond:exp);
    dsd <- DESeq(ss, test="LRT", reduced=~ cond + exp);
    res <- results(dsd, alpha=.05);
    print(length(which(res$padj <= 0.05)))
    return(list("res"=res, "ds"=dsd))
}

hpi24i <- hpi_intest(hpi24_c)

hpi48i <- hpi_intest(hpi48_c)

hpi72i <- hpi_intest(hpi72_c)

na_conv <- function(df) {
    ndf <- df;
    ndf[ndf > 0.05] <- NA;
    return(ndf)
}

ds_v <- list(rh24$res$padj,rh48$res$padj, rh72$res$padj,
             sh24$res$padj, sh48$res$padj, sh72$res$padj,
             rtc$res$padj, stc$res$padj,
             hpi24i$res$padj, hpi48i$res$padj, hpi72i$res$padj)

ds_p <- rapply(ds_v, na_conv, how="list")

main_df <- data.frame(
    "rh24"=ds_p[[1]],
    "rh48"=ds_p[[2]],
    "rh72"=ds_p[[3]],
    "sh24"=ds_p[[4]],
    "sh48"=ds_p[[5]],
    "sh72"=ds_p[[6]],
    "rtc"=ds_p[[7]],
    "stc"=ds_p[[8]],
    "hpi24i"=ds_p[[9]],
    "hpi48i"=ds_p[[10]],
    "hpi72i"=ds_p[[11]]
)

#Remove rows containing only NA
final_df <- main_df[rowSums(is.na(main_df)) != 11,]

nrow(final_df)

write.csv(final_df, file="DEG.csv")



c24 <- getpair_stats(gdds, "c24")

c48 <- getpair_stats(gdds, "c48")

c72 <- getpair_stats(gdds, "c72")

plotFC(c24$res, "este")
plotFC(c48$res, "este")
plotFC(c72$res, "este")






