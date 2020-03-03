suppressPackageStartupMessages({
    library(getopt)
    library(DESeq2)
    library(rtracklayer)
    library(GenomicFeatures)
    library(WriteXLS)
    library(dplyr)
})

basic_functions_file = "basic_functions.R"
source(basic_functions_file)


args_return=Sys.getenv("args_return")

helpVals = matrix(c(
	"t", "title", 1, "Title",
	"g","gtf", 1, "gtf file using which the count table was generated",
	"f","inputFile",1, "count table",
	"0","cmd0",1, "file containing comparisons",
	"s","sampleFile", 1, "file containing <name> <condn> <library type>",
	"O","outpt", 1, "output directory"
	),
	ncol = 4, byrow = T
)

quitMessage = function(){
 message("

SampleFile:  <name> <condn> <library type>
BCOR_shRNA1.star.ReadsPerGene.out.tab	BCOR_shRNA1	BCOR_shRNA	paired-end
BCOR_shRNA2.star.ReadsPerGene.out.tab	BCOR_shRNA2	BCOR_shRNA	paired-end
Control_shRNA1.star.ReadsPerGene.out.tab	Control_shRNA1	Control	paired-end
Control_shRNA2.star.ReadsPerGene.out.tab	Control_shRNA2	Control	paired-end
Note: <name> will not be used in this case.

cmd0: comparisons.txt
Control BCOR_shRNA

")
}
if(args_return > 0){
    Help(helpVals)
    quitMessage()
    quit(status = 2)
}
Help(helpVals)

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Required Args:

title: 
gtfFile: gtf file using which the count table was generated
inputFile: count table
sampleFile: file containing <name> <condn> <library type>
cmd0: file containing comparisons

SampleFile:  <name> <condn> <library type>
BCOR_shRNA1.star.ReadsPerGene.out.tab	BCOR_shRNA1	BCOR_shRNA	paired-end
BCOR_shRNA2.star.ReadsPerGene.out.tab	BCOR_shRNA2	BCOR_shRNA	paired-end
Control_shRNA1.star.ReadsPerGene.out.tab	Control_shRNA1	Control	paired-end
Control_shRNA2.star.ReadsPerGene.out.tab	Control_shRNA2	Control	paired-end
Note: <name> will not be used in this case.

cmd0: comparisons.txt
Control BCOR_shRNA

		")
    quit(status = 2)
}

title = Sys.getenv("title")
gtfFile=Sys.getenv("gtf")
countTable=Sys.getenv("inputFile")
sampleMatrix = Sys.getenv("sampleFile")
comparisonsFile = Sys.getenv("cmd0")
outpt = Sys.getenv("outpt")

# gtfFile = "~/Volumes/share/reference/gtfs/grcm38_mm10.gencode.vM16.annotation.gtf"
# countTable = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/counts.rev.frFS.Illum.tsv"
# sampleMatrix = "~/Volumes/scratch119/GONIA/2018.08.Setd1bMouse_Gonia_RNASeq/STAR/samples.counttable.txt"
# c1="Empty_vector_day5"
# c2="Ehbp1l1_gRNA_day5"
#
" Read contrast matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
sample.matrix = read.table(sampleMatrix)
colnames(sample.matrix) = c("bam","names", "condition", "libtype")
rownames(sample.matrix) = sample.matrix$names
sample.matrix = sample.matrix [, c("condition","libtype") ]
#head(sample.matrix)
" conditions "
conditions = as.character(unique((sample.matrix$condition)))
message("[ Info ] Conditions =  ", paste0(conditions, sep=" "))
sample.matrix

" Read in GTF file "
txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
gtfData <- import.gff2(gtfFile)
annot_data <- as.data.frame(gtfData)
annot_data <- annot_data[ annot_data$type == "gene", ]
annot_data <- data.frame("gene_name" = annot_data$gene_name,
                         "gene_id" = annot_data$gene_id,
                         "coords" = paste0(annot_data$seqnames, ":", annot_data$start, "-", annot_data$end), 
                         "strand" = annot_data$strand,
                         "gene_type" = annot_data$gene_type
)
annot_data = unique(annot_data)
genesToConsider = as.character(annot_data$gene_id)

" Read count table "
x = read.table(countTable, head=T, row.names=1)
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
"|--> Genes before filtering"
nrow(counttable)
counttable <- counttable [ rownames(counttable) %in% genesToConsider, ]
"|--> Genes after filtering"
nrow(counttable)
# " Read in GTF file "
# txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", circ_seqs = character())
# ebg <- exonsBy(txdb, by="gene")
# genesToConsider <- names(ebg)
# 


plot.volcano <- function(res, lfc, pval){
    tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
    head(tab)
    par(mar = c(5, 4, 4, 4))
    plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
    signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
    points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
    abline(h = -log10(pval), col = "green3", lty = 2)
    abline(v = c(-lfc, lfc), col = "blue", lty = 2)
    mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
    mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
}


medianNormCounts <- function(dds, colData){
    normCounts = counts(dds, normalize = T)
    conditions = unique(as.character(colData$condition))
    countMatrix = c()
    for(cond in conditions){
        reqCols = rownames(colData[colData$condition == cond,])
        temp = rowMedians(normCounts [, reqCols])
        countMatrix = cbind(countMatrix, temp)
    }
    rownames(countMatrix) = rownames(normCounts)
    colnames(countMatrix) = conditions
    countMatrix = round(countMatrix, 3)
    return(countMatrix)
}

getDESeq2 <- function(counttable, sample.matrix, annot_data, c1, c2){
    
    " Extract data corresponding to the conditions"
    message("reference = ", c1)
    message("test = ", c2)
    colData = subset(sample.matrix, condition %in% c(c1,c2))
    colData$condition = droplevels(colData$condition) # Drop unwanted levels
    colData$condition = relevel(colData$condition, ref = c1) # Keep first condition as the reference level
    
    " Filter count data for low count data, reads/gene > 0"
    countdata <- counttable
   # message(rownames(colData))
    countdata = as.matrix(countdata [ , rownames(colData)])
    colData
   
    rmax = apply(countdata, 1, max)
   # message("[ Info ] Total genes = ", nrow(countdata))
    countdata = countdata [ rmax >= 1, ]
#    message("[ Info ] Genes with reads = ", nrow(countdata))
    
    
    " -- Check if the row order in contrast and columns in count data are in the same order"
    if( ! all(rownames(colData) == colnames(countdata)) ){
        message("[ Error ] Total rows in filtered contrast matrix != Total columns in count table")
        message("[ Error ] Aborting")
        save.image("deseq_error.RData")
        quite(save="no")
    } else {
        message("[ Info ] All seems correct. Proceeding forward ...")
    }
    
    " Create DESeq object"
    condition = colData$condition
    condition = relevel(condition, ref = c1)
    dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~ condition)
    #dds$condition <- factor(dds$condition, levels=c(c1,c2))
   # print(condition)
    " Get gene length"
    " -- Including exon granges in rowRanges. You can also add in the actual gene sizes in mcols(dds)$basepairs"
    temp <- rowRanges(dds)
    rangesToShoveIn <- ebg[names(temp)] 
    rowRanges(dds) <- rangesToShoveIn
    
    " Perform DESeq "
    dds <- DESeq(dds)

    " Get results "
    #res <- results(dds, alpha=0.05, format = "DataFrame")
    res <- results(dds, format = "DataFrame")
    medCounts = medianNormCounts(dds, colData)
    outpt = merge(medCounts, res,  by = "row.names")
    outpt = merge(annot_data, outpt, by.y = "Row.names", by.x = "gene_id")
    # summary(res)
    " Extract results"
    # genesInRes = as.character(rownames(res))
   
    outpt.list = list(
        "dds" = dds,
        "outpt" = outpt
    )
    # outpt.list[["dds"]] = dds
    # outpt.list[["outpt"]] = outpt
    # 
    " -- Data from DESeq "
    
    return(outpt.list)
}

dir.create(outpt)

message("Reading comparisons file ...")
comparisons = read.table(comparisonsFile)
colnames(comparisons) = c("condition1", "condition2")

allComps = list(comps = comparisons)
alldds = list()
for(i in 1:nrow(comparisons)){
    c1 = as.character(comparisons[i,1])
    c2 = as.character(comparisons[i,2])
    message("-- Comparing ", c1, " and ", c2)
    res.results = getDESeq2(counttable, sample.matrix, annot_data, c1, c2)
    n = paste0(c1,"__",c2)
    t = paste0("comp",i)
    allComps[[t]] = res.results$outpt[order(res.results$outpt$padj, decreasing = F), ]
    alldds[[t]] = res.results$dds
    temp = colData(alldds[[t]])
    comment(allComps[[t]]$gene_id) <- paste(temp$condition, rownames(temp), sep=":", collapse=", ")
    write.csv(allComps[[t]], file = paste0(outpt, "/",n,".csv" ))
}

" Write data to excel sheet"
WriteXLS(allComps, ExcelFileName=paste(outpt,"/",title,".deseqAll.xlsx",sep=""), row.names=T, BoldHeaderRow = T, AdjWidth = T)

# pdf(paste(title,".deseq.pdf",sep=""),5,5)
# plotMA(res, frame=F, alpha=0.05)
# plot.volcano(res, 1, 0.05)
# plotDispEsts(dds)
# dev.off()

save.image(paste0(outpt,"/",title,".deseq2.RData"))
q(status = 0)


# > head(counttable)
# BCOR_shRNA1 BCOR_shRNA2 Control_shRNA1 Control_shRNA2 Control2_shRNA1 Control2_shRNA2 Ring1_RNF2_shRNA1
# ENSG00000223972.5           1           0              0              0               0               0                 0
# ENSG00000227232.5           1           0              0              0               1               0                 1
# ENSG00000278267.1           0           0              0              0               0               0                 0
# ENSG00000243485.5           0           0              0              0               2               0                 2
# ENSG00000284332.1           0           0              0              0               0               0                 0
# ENSG00000237613.2           0           0              0              0               0               0                 0
# Ring1_RNF2_shRNA2
# ENSG00000223972.5                 0
# ENSG00000227232.5                 0
# ENSG00000278267.1                 0
# ENSG00000243485.5                 1
# 
# > head(sample.matrix)
# condition    libtype
# BCOR_shRNA1     BCOR_shRNA paired-end
# BCOR_shRNA2     BCOR_shRNA paired-end
# Control_shRNA1     Control paired-end
# Control_shRNA2     Control paired-end
# Control2_shRNA1   Control2 paired-end
# Control2_shRNA2   Control2 paired-end
# 
# > head(annot_data)
# gene_name           gene_id           coords strand  source
# 1     DDX11L1 ENSG00000223972.5 chr1:11869-14409      +  HAVANA
# 2      WASH7P ENSG00000227232.5 chr1:14404-29570      -  HAVANA
# 3   MIR6859-1 ENSG00000278267.1 chr1:17369-17436      - ENSEMBL
# 4 MIR1302-2HG ENSG00000243485.5 chr1:29554-31109      +  HAVANA
# 5   MIR1302-2 ENSG00000284332.1 chr1:30366-30503      + ENSEMBL
# 6     FAM138A ENSG00000237613.2 chr1:34554-36081      -  HAVANA
# 
#res.results = getDESeq2(counttable, sample.matrix, annot_data, "Control", "BCOR_shRNA")


