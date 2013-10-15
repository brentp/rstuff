intensity.450k = function(targets, build="hg19", id_col=1){
    library('minfi')
    rgset = read.450k.exp(targets=targets, extended=TRUE)
    #locs = mapToGenome(rgset, genomeBuild=build)
    locs = getLocations(rgset)
    #rgset.pf = pfilter(rgset)
    mset = preprocessRaw(rgset)
    intensity = getMeth(mset) + getUnmeth(mset)
    colnames(intensity) = targets[,id_col]
    stopifnot(nrow(intensity) == length(locs@seqnames))
    rownames(intensity) = paste(locs@seqnames, locs@ranges@start, sep=":")
    return(intensity)
}

# this is copied from the bioconductor ChAMP package and simplified
# for my understanding. ... and parallellized.
cna.450k = function(targets, prefix="cn.450k", mc.cores=6){
    library(DNAcopy)
    library(preprocessCore)
    library(parallel)

    intensity = intensity.450k(targets)
    intensity = intensity[grep('chr[Y|X]', rownames(intensity),
                               perl=TRUE, invert=TRUE),]

    #intensity = read.mat('intensity.txt')
    #write.table(intensity, row.names=T, sep="\t", file="intensity.txt", quote=F)
    samples = colnames(intensity)
    intsqn = log2(normalize.quantiles(as.matrix(intensity)))
    colnames(intsqn) = samples
    
    chrom = unlist(lapply(strsplit(rownames(intensity), ":", fixed=TRUE),
                           function(r) r[[1]]))
    pos = unlist(lapply(strsplit(rownames(intensity), ":", fixed=TRUE),
                 function(r) as.integer(r[[2]])))

    #refs = rowMeans(intsqn) # use the mean of all samples as the reference
    # median should handle places where we really see a big difference in
    # cases
    refs = apply(intsqn, 1, median)

    tmp = mclapply(1:ncol(intsqn), function(i){
        log.ratio = intsqn[, i, drop=FALSE] - refs
        CNA.obj = CNA(log.ratio, chrom, pos, data.type="logratio",
                            sampleid=samples[i])
        CNA.obj = smooth.CNA(CNA.obj)
        CNA.obj = segment(CNA.obj, verbose=1, alpha=0.005,
                          undo.splits="sdundo", undo.SD=1.96)

        png(paste0(prefix, samples[i], ".png"))
        plot(CNA.obj, plot.type="w")
        dev.off()

        segs = CNA.obj$output
        segs = segs[order(segs$chrom, segs$loc.start),!colnames(segs) == "ID"]
        write.table(segs, sep="\t", col.names=T, row.names=F, quote=F, 
                    file=paste0(prefix, samples[i], ".txt"))
        message(paste0(prefix, samples[i], ".txt"))
        return(TRUE)
    }, mc.cores=mc.cores)
}
