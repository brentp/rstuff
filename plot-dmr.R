
library(Gviz)
library(GenomicRanges)
library(lattice)
source('~/src/rstuff/rmodels.R')
options(stringsAsFactors=FALSE)
data(cpgIslands)

args = commandArgs(TRUE)

meth_bed = args[1]
clin = args[2]
id_col = args[3]
groups_name = args[4]
genome = args[5]

mat = read.mat(meth_bed)
cov = read.delim(clin)
cov$lbl = cov[,id_col]
cov$grp = cov[,groups_name]

shared = intersect(cov$lbl, colnames(mat)[1:ncol(mat)])

locs = data.frame(chrom = unlist(lapply(strsplit(rownames(mat), ":", fixed=TRUE), function(r) { r[1] })))
locs$start = unlist(lapply(strsplit(rownames(mat), ":", fixed=TRUE), function(r) { as.numeric(r[2]) - 1 }))
locs$end = unlist(lapply(strsplit(rownames(mat), ":", fixed=TRUE), function(r) { as.numeric(r[2]) }))


mat = mat[, shared]
gen = genome(cpgIslands)

rownames(cov) = cov$lbl
cov = cov[shared,]

groups = as.character(cov$grp)

for (dmr in args[6:length(args)]){

    gc()


    png_name = paste(sub(":", "-", dmr), id_col, groups_name, genome, "png", sep=".")
    png(png_name, width=1500, height=1200)
    message(png_name)


    chrom = unlist(strsplit(dmr, ":"))[1]
    se = unlist(strsplit(dmr, ":"))[2]
    dmr = unlist(lapply(strsplit(se, "-"), as.integer))


    m1 = mat[(chrom == locs[,"chrom"]) & (dmr[1] < locs[,"end"]) & (dmr[2] > locs[,"start"]),,drop=FALSE]
    lm1 = locs[(chrom == locs[,"chrom"]) & (dmr[1] < locs[,"end"]) & (dmr[2] > locs[,"start"]),,drop=FALSE]

    dat = matrix(as.numeric(m1), nrow=ifelse(is.null(nrow(m1)), length(m1), nrow(m1)))

    colnames(dat) = colnames(m1)[1:ncol(dat)]
    dat = t(dat) #[,shared])
    rownames(dat) = colnames(m1)[1:ncol(m1)]
    colnames(dat) = rownames(m1)



    cpgtrack = AnnotationTrack(cpgIslands, name="CpG", chromosome=chrom)
    gtrack = GenomeAxisTrack(range=IRanges(start=dmr[1], end=dmr[2]), fill.range="darkblue", add53=TRUE, add35=TRUE, cex=1.4)
    itrack = IdeogramTrack(genome=genome, chromosome=chrom)


    ranges = GRanges(seqnames=rep(chrom, nrow(m1)),
         ranges=IRanges(start=as.numeric(lm1[,"start"]), end=as.numeric(lm1[,"start"] ) + 10))



    colnames(dat) = paste(lm1[,1], lm1[,3], sep=":")


    details = function(identifier, ...){
        d = data.frame(signal=dat[,identifier], group=groups)
        #write(paste(dim(d)), stderr())
        print(densityplot(~signal, group=group, data=d,
                 main=list(label=identifier, cex=0.7),
                           scales=list(draw=FALSE, x=list(draw=TRUE)),
                           ylab="", xlab="",

              col=c("blue", "darkgreen"),
              col.line=c("blue", "darkgreen"),
              col.symbol=c("blue", "darkgreen"),
              jitter.y=TRUE,
              alpha=0.7,
              cex=0.7,
              plot.points="rug",
                           ), newpage=FALSE,
                          prefix="plot")
    }
    aptrackEven = AnnotationTrack(range=ranges, id=colnames(dat), fun=details, 
                                  selectFun=function(index, ...){ 
                                      index %in% seq(2, 100, by=2)
                                  },
                                  details.size=0.96,
                                  )

    aptrackOdd = AnnotationTrack(range=ranges, id=colnames(dat), fun=details, 
                                  selectFun=function(index, ...){ 
                                      index %in% seq(1, 100, by=2)
                                  },
                                  details.size=0.96,
                                  )

    dtrack = DataTrack(range=ranges,
                  dat=dat,
                   groups=groups,
                    genome=genome, name="methylation",
                    type=c("boxplot"),
                    col=c("blue", "darkgreen"),
                   col.frame=c(NA, NA),
                    box.width=14, legend=TRUE, cex.legend=2
                    )


    plotTracks(list(itrack, gtrack, dtrack, cpgtrack, aptrackEven, aptrackOdd),
               from=(dmr[1] - 40), to=(dmr[2] + 40),
               fill=c("blue", "darkgreen"),
               background.title="darkblue", background.panel="#EFEFEF")

    #warnings()
    dev.off()
}
