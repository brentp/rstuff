
library(Gviz)
library(GenomicRanges)
source('rmodels.R')
options(stringsAsFactors=FALSE)
data(cpgIslands)

args = commandArgs(TRUE)

dmr = args[1]
meth_bed = args[2]
clin = args[3]
id_col = args[4]
groups = args[5]
genome = args[6]


png_name = paste(sub(":", "-", dmr), id_col, groups, genome, "png", sep=".")
print(png_name)
png(png_name, width=1500, height=1200)


chrom = unlist(strsplit(dmr, ":"))[1]
se = unlist(strsplit(dmr, ":"))[2]
dmr = unlist(lapply(strsplit(se, "-"), as.integer))

mat = read.mat(meth_bed)
cov = read.delim(clin)
cov$lbl = cov[,id_col]
cov$grp = cov[,groups]

shared = intersect(cov$lbl, colnames(mat)[3:ncol(mat)])
mat = mat[,c("start", "end", shared)]
print(dim(mat))

rownames(cov) = cov$lbl
cov = cov[shared,]
print(dim(cov))

groups = as.character(cov$grp)
print(length(groups))
print(groups)

gen = genome(cpgIslands)

m1 = mat[(dmr[1] < mat[,"end"]) & (dmr[2] > mat[,"start"]),]

dat = matrix(as.numeric(m1[,3:ncol(m1)]), nrow=nrow(m1))
print("OK0")
print(dim(dat))
print(dim(m1))

colnames(dat) = colnames(m1)[3:ncol(m1)]
dat = t(dat) #[,shared])
print("OKa")
rownames(dat) = colnames(m1)[3:ncol(m1)]
print("OKb")
colnames(dat) = rownames(m1)

print("OK1")

cpgtrack = AnnotationTrack(cpgIslands, name="CpG", chromosome=chrom)
gtrack = GenomeAxisTrack(range=IRanges(start=dmr[1], end=dmr[2]), fill.range="darkblue", add53=TRUE, add35=TRUE, cex=1.4)
itrack = IdeogramTrack(genome=genome, chromosome=chrom)

print("OK2")

ranges = GRanges(seqnames=rep(chrom, nrow(m1)),
     ranges=IRanges(start=as.numeric(m1[,"start"]), end=as.numeric(m1[,"start"] ) + 10))


library(lattice)
print("OK3")

colnames(dat) = paste(m1[,1], m1[,2], sep=":")


details = function(identifier, ...){
    d = data.frame(signal=dat[,identifier], group=groups)
    write(paste(dim(d)), stderr())
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
print("OK4")
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
