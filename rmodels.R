library(limma)
library(gtools)
library(data.table)

options(scipen=4, stringsAsFactors=FALSE) # stop it from printing 1e6 instead of 1000000

logit = defmacro(p, expr=log(p) - log(1 - p))

genomic_control = function(pvals, mid.fn=median){
    # http://en.wikipedia.org/wiki/Population_stratification
    mid.fn(qchisq(pvals, df=1, lower.tail=F)) / 0.4549
}

spia.ez = function(symbols, values, top){
    if(length(symbols) != length(values)){
        stop("must send in a list of gene symbols with corresponding (logFC) values")
    }
    if(length(top) != length(values)){
        stop("must send in a list booleans with TRUE indicating genes of interest")
    }

    library(SPIA)
    library(org.Hs.eg.db)

    entrez = select(org.Hs.eg.db, keys=symbols, keytype="SYMBOL", cols="ENTREZID")
    # subset to genes that have an entrez mapping
    names(values) = symbols
    names(symbols) = symbols
    names(top) = symbols


    # select to those that were converted
    values = values[entrez[,1]]
    symbols = symbols[entrez[,1]]
    top = top[entrez[,1]]

    save = !duplicated(entrez[,2])
    values = values[save]
    symbols = symbols[save]
    top = top[save]

    entrez.ids = entrez[save, 2]
    names(values) = entrez.ids
    names(symbols) = entrez.ids
    names(top) = entrez.ids
    print(head(top))
    print(head(values[top]))
    print(head(entrez.ids))
    spia(de=values[top], all=entrez.ids, nB=200, organism="hsa") #, plots=TRUE)
}

# return pvalues after the genomic control correction
gc.adjust = function(pvals){
    gc = genomic_control(pvals)
    qchi = qchisq(pvals, df=1, lower.tail=F)
    qchi = qchi / gc
    pchisq(qchi, df=1, lower.tail=F)
}

.write.raw.meth = function(RGset, dp, out_prefix){
    raw = preprocessRaw(RGset)
    colnames(dp) = paste("detection-p-", colnames(dp), sep="")
    methylated = getMeth(raw)
    colnames(methylated) = paste("methylated-", colnames(methylated), sep="")
    unmethylated = getUnmeth(raw)
    colnames(unmethylated) = paste("unmethylated-", colnames(unmethylated), sep="")
    dpout = cbind(unmethylated, methylated, format(dp, digits=3, trim=TRUE))
    write.table(dpout, sep="\t", quote=FALSE, file=paste(out_prefix, "raw-values.txt", sep=""), row.names=T)
    rm(methylated); rm(unmethylated); gc();
}

read.tab = function(fname, sep="\t", header=TRUE, ...){
    if(class(fname) != "character"){
        return(fname)
    }
    clin = read.delim(fname, sep=sep, header=TRUE, comment.char="", ...)
    return(clin)
}

read.mat = function(fname, sep="\t"){ 
    # much faster way to read a matrix.
    #  m = read.delim(fname, row.names=row.names, header=header, sep=sep,
    #                                 comment.char="", ...)
    library(MatrixEQTL);
    s = SlicedData$new();
    s$fileDelimiter = sep;
    s$fileOmitCharacters = 'NA';
    s$fileSkipRows = 1;
    s$fileSkipColumns = 1;
    s$fileSliceSize = 25000;
    s$LoadFile(fname);
    return(as.matrix(s))
}

shared_rows_cols = function(a, b, acolname=NA){
    bnames = colnames(b)
    if(is.na(acolname)){
        anames = rownames(a)
    } else {
        anames = a[,acolname]
    }
    shared = intersect(anames, bnames)
    if(length(shared) == 0){ stop("no shared ids") }
    return(shared)
}


normalize.450K.method = function(targets, method=c("dasen", "swan"), id_col=1, prefix=NULL){
    library('minfi')
    rgset = read.450k.exp(targets=targets, extended=TRUE)

    qcReport(rgset, sampNames=targets$StudyID, sampGroups=targets$Sample_Plate)
    if(!is.null(prefix)){
        clin2 = pData(rgset)
        write.table(clin2, file=paste0(prefix, "clin.txt"), sep="\t", row.names=F, quote=F)
    }
    m.norm = NA
    if(method == "dasen"){
        library('wateRmelon')
        rgset.pf = pfilter(rgset)
        m.norm = Beta2M(dasen(rgset.pf))
    } else {
        detP = detectionP(rgset)
        failed = detP > 0.01
        bad_probes = rowMeans(failed) > 0.10
        rgset.pf = preprocessSWAN(rgset)[!bad_probes,]
        failed = failed[!bad_probes,]
        m.norm = getM(rgset.pf)
        message(sprintf("setting %i probes with detectionP > 0.01 to NA", sum(failed)))
        m.norm[failed] = NA
    }
    locs = getLocations(rgset.pf)

    message(sprintf("seqnames length: %d, dim: %d, %d", length(locs@seqnames),
                nrow(locs@ranges@start), ncol(locs@ranges@start)))

    rownames(m.norm) = paste(locs@seqnames, locs@ranges@start, sep=":")
    colnames(m.norm) = targets[, id_col]
    snames = as.character(locs@seqnames)
    starts = as.integer(locs@ranges@start)
    m.norm = m.norm[order(snames, starts),]
    if(!is.null(prefix)){
        write.matrix(m.norm, file=paste0(prefix, method, ".M.txt"))
    }
    return(m.norm)
}

normalize.450k = function(fclin, out_prefix, base_path, id_col=1){
    library(minfi)
    require('IlluminaHumanMethylation450kmanifest')

    clin = read.tab(fclin)

    idats = list.files(base_path, pattern="*_Grn.idat", recursive=TRUE, full.names=T)
    idats = data.frame(names=unlist(strsplit(basename(idats), "_Grn.idat")), Basename=idats)
    clin = merge(clin, idats, by.x=id_col, by.y="names")

    RGset = read.450k.exp(base = ".", targets = clin)
    pd = pData(RGset)

    qcReport(RGset, sampNames=pd[,id_col], pdf=paste(base_path,
        "qcReport.pdf", sep="."), sampGroups=rep(1, ncol(RGset)))

    Mset.swan <- preprocessSWAN(RGset)

    # 
    dp = detectionP(RGset)
    .write.raw.meth(RGset, dp, out_prefix)

    beta = getBeta(Mset.swan)
    colnames(beta) = pd[,id_col]
    message("writing beta ...")
    write.matrix(beta, file=paste(out_prefix, "beta.txt", sep=""))
    rm(beta)
    gc()

    M = getM(Mset.swan)
    colnames(M) = pd[,id_col]

    message("writing M...")
    write.matrix(M, file=paste(out_prefix, "M.txt", sep=""))
    M[dp > 0.05] = NA
    message("values > 0.05:", sum(dp > 0.05))
    write.matrix(M, file=paste(out_prefix, "M.pgt05.txt", sep=""))
    return(M)
}

normalize.charm = function(sample_description, out_prefix, id_col=1, subject=NA){
    library(charm)
    if(is.na(subject)){
        require(BSgenome.Hsapiens.UCSC.hg18)
        subject=Hsapiens
    }

    pd = read.tab(sample_description)
    rawData =  readCharm(files=pd$filename, sampleKey=pd, path="")

    ################################
    # remove arrays with low quality
    ################################
    qual = qcReport(rawData, file=paste(out_prefix, "qc-report.pdf", sep=""))
    qc.min = 70
    rawData = rawData[,qual$pmSignal>=qc.min]
    qual=qual[qual$pmSignal>=qc.min,]
    pd=pd[pd$sampleID%in%rownames(qual),]
    pData(rawData)$qual=qual$pmSignal
    rm(qual)

    ctrlIdx = getControlIndex(rawData, subject=subject, noCpGWindow=600)
    p  = methp(rawData, controlIndex=ctrlIdx, betweenSampleNorm="quantile")
    message("DIMS:")
    message(paste(dim(p), collapse=" "))

    pmq = pmQuality(rawData)
    okqc = which(rowMeans(pmq) > 75)

    chr = pmChr(rawData)
    pns = probeNames(rawData)
    pos = pmPosition(rawData)
    pd  = pData(rawData)
    rm(rawData)


    ################################
    # remove probes with low quality
    ################################
    Index = setdiff(okqc, ctrlIdx)
    Index = Index[order(chr[Index], pos[Index])]

    p = p[Index,]
    chr = chr[Index]
    pos = pos[Index]
    pns = pns[Index]
    pns = clusterMaker(chr,pos) 

    message("DIMS:")
    message(paste(dim(p), collapse=" "))

    rm(ctrlIdx)
    rownames(p) = paste(chr, pos, sep=":")
    colnames(p) = pd[,id_col]

    message("writing...")
    write.matrix(p, file=paste(out_prefix, ".methp.txt", sep=""))
    write.table(pd, sep="\t", quote=FALSE, file=paste(out_prefix, ".clinical.txt", sep=""), row.names=F)
    return(p)
}

get_null_model = function(model){
    full_formula = as.formula(model)
    null_formula = attr(terms(full_formula), "term.labels")
    # null model is the full model without variable of interest.
    if(length(null_formula) > 1){
        null_formula = as.formula(paste("~", paste(null_formula[2:length(null_formula)], collapse=" + ")))
    } else {
        null_formula = as.formula("~ 1")
    }
    return (null_formula)
}

err.log = function(...){
    write(paste("#>", ..., sep=" "), stderr())
}


# remove the batch effects defined in svs (from peer or sva) from the data.
# from A Jaffe: http://permalink.gmane.org/gmane.science.biology.informatics.conductor/42857
cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

sva.limma.ez = function(data, clin, model,
                prefix="results/sva.limma.",
                contrasts=NULL, probe_len=1,
                batch_correct=TRUE){
    coef = 2

    # the first term is always the variable of interest.
    full_formula = as.formula(model)
    err.log("full model:", full_formula)
    null_formula = get_null_model(model)
    err.log("null model:", null_formula)


    # remove any rows with na
    complete = complete.cases(clin[,attr(terms(full_formula), "term.labels")])
    err.log("removing:", sum(!complete), "because of missing data")
    err.log("leaving:", sum(complete), "rows of data.")
    data_complete = as.matrix(data[,complete])
    rm(data); gc()

    library(Matrix)
    mod = model.matrix(full_formula, data=clin)

    if(batch_correct){
        library(sva)
        mod0 = model.matrix(null_formula, data=clin)

        stopifnot(nrow(mod) == ncol(data_complete))
        err.log("starting sva to find surrogate variables")
        svobj = sva(data_complete, mod, mod0)
        err.log("number of surrogate variables:", svobj$n.sv)
        err.log("running limma with original model and surrogate variables")
        # add the surrogate variables to the model.

        colnames(svobj$sv) = paste('sv_', 1:ncol(svobj$sv), sep="")
        modsv = cbind(mod, svobj$sv)
        write.table(modsv, sep="\t", row.names=F, file=paste(prefix,
                        "mod.sva.txt", sep=""), quote=FALSE)

        data_complete = cleanY(data_complete, mod, svobj$sv)
        write.matrix(data_complete, file=paste(prefix, "cleaned.y.txt", sep=""))
    }

    fit = limma.ez(data_complete, mod, coef, contrasts, prefix, probe_len)
    return(fit)
}

# stolen from genefilter
rowVars = function (x, ...) {
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1))
}

# if ids has repeated values (e.g. genes. then keep only the corresponding
# rows with highest variance.
keep_highest_var = function(ids, mat){
   rownames(mat) = as.character(1:nrow(mat))
   dup_ids = ids[duplicated(ids)] # 266 unique
   is_dup = ids %in% dup_ids # 547 total
   vars = rowVars(mat[is_dup,])

   max.vars = tapply(vars, ids[ids %in% dup_ids], function(f){ names(which.max(f)) })
   keep = ids[!is_dup]
   stopifnot(intersect(names(max.vars), keep) == character(0))

   ret = mat[max.vars,]
   ret = rbind(ret, mat[!is_dup,])
   rownames(ret) = c(names(max.vars), ids[!is_dup])

   as.matrix(ret)
}

remove_low_variance = function(mat, p_drop=0.25){
    if(p_drop > 1){ p_drop = p_drop / 100 }
    rvm = rowVars(mat)
    mat[rvm > sort(rvm)[p_drop * length(rvm)],]
}

read.450k.manifest = function(fname, chrom_prefix="chr", version=c("hg19", "hg18")){
    version=match.arg(version, c("hg19", "hg18"))
    print(version)
    df = read.delim(fname, sep=",", skip=7, quote="", stringsAsFactors=FALSE)
    df$hg18_pos = paste0(chrom_prefix, df$Chromosome_36, ":", df$Coordinate_36)
    df$hg19_pos = paste0(chrom_prefix, df$CHR, ":", df$MAPINFO)
    df$Coordinate_36 = as.integer(df$Coordinate_36)
    df$MAPINFO = as.integer(df$MAPINFO)
    if(version == "hg18"){
        pos = data.frame(chrom=paste0(chrom_prefix, df$Chromosome_36), start=df$Coordinate_36 - 1, end=df$Coordinate_36)
    } else {
        pos = data.frame(chrom=paste0(chrom_prefix, df$CHR), start=df$MAPINFO - 1, end=df$MAPINFO)
    }
    df = data.frame(pos, df)
    df[order(df$chrom, df$start),]
}


write.matrix = function(mat, file, name="probe", quote=FALSE, sep="\t", digits=3, ...){
    mat = cbind(rownames(mat), format(round(mat, digits=digits), digits=digits, trim=TRUE))
    colnames(mat)[1] = name
    write.table(mat, file=file, quote=quote, sep=sep, row.names=FALSE, ...)
    gc();
}    


peer.limma.ez = function(data, clin, model=NULL,
                prefix="results/peer.limma.",
                contrasts=NULL, probe_len=1,
                batch_correct=5){
    #https://github.com/PMBio/peer/wiki/Tutorial
    coef = 2

    if(!is.null(model)){
    # the first term is always the variable of interest.
        full_formula = as.formula(model)
        err.log("full model:", full_formula)

        # remove any rows with na
        complete = complete.cases(clin[,attr(terms(full_formula), "term.labels"), drop=TRUE])
        mod  = model.matrix(full_formula, data=clin)
    } else { # they sent in a matrix, not a formula
        complete = complete.cases(clin)
        mod  = as.matrix(clin[complete,], ncol=ncol(clin))
    }

    err.log("removing:", sum(!complete), "because of missing data")
    err.log("leaving:", sum(complete), "rows of data.")
    data_complete = as.matrix(data[,complete])
    rm(data); gc()



    if(!is.na(batch_correct) && as.logical(batch_correct)){
        n_factors = as.integer(batch_correct)
        data_complete = run.peer(mod, data_complete, prefix, n_factors)
        prefix = paste(prefix, batch_correct, sep=".")
    }
    # including all peer factors, just include those before 1/alpha levels off?
    fit = limma.ez(data_complete, mod, coef, contrasts, prefix, probe_len)
    return(fit)
}


peer.factors = function(mod, data_complete, n_factors=5){
    library(peer)
    peer_obj = PEER()
    stopifnot(nrow(mod) == ncol(data_complete))
    probes = rownames(data_complete)
    data_complete = t(data_complete)
    
    # set sensible defaults
    PEER_setCovariates(peer_obj, mod)
    PEER_setTolerance(peer_obj, 1e-10)
    PEER_setVarTolerance(peer_obj, 1e-12)
    PEER_setNmax_iterations(peer_obj, 1000)
    PEER_setPhenoMean(peer_obj, data_complete)

    nK = min(n_factors, nrow(data_complete) - ncol(mod) - 1)
    PEER_setNk(peer_obj, nK)
    PEER_update(peer_obj)

    X = PEER_getX(peer_obj)
    modpeer = X[,(ncol(mod) + 1):(ncol(mod) + nK)]

    colnames(modpeer) = paste('peer_', 1:ncol(modpeer))
    rownames(modpeer) = rownames(mod)
    return(modpeer)
}

run.peer = function(mod, data_complete, prefix, n_factors=5){
    library(peer)
    peer_obj = PEER()
    stopifnot(nrow(mod) == ncol(data_complete))
    probes = rownames(data_complete)
    data_complete = t(data_complete)

    # set sensible defaults
    PEER_setCovariates(peer_obj, mod)
    PEER_setTolerance(peer_obj, 1e-10)
    PEER_setVarTolerance(peer_obj, 1e-10)
    PEER_setNmax_iterations(peer_obj, 1000)

    #PEER_setTolerance(peer_obj, 1e-4)
    #PEER_setVarTolerance(peer_obj, 1e-12)
    #PEER_setNmax_iterations(peer_obj, 195)
    PEER_setPhenoMean(peer_obj, data_complete)

    nK = min(n_factors, nrow(data_complete) - ncol(mod) - 1)
    PEER_setNk(peer_obj, nK)
    err.log("set number of factors to infer as", nK)

    PEER_update(peer_obj)
    modpeer = PEER_getX(peer_obj)
    wgt_all = PEER_getW(peer_obj)

    f = paste(prefix, "peer.alpha.txt", sep="")
    write.table(PEER_getAlpha(peer_obj), file=f, row.names=F, sep="\t")

    colnames(modpeer) = c(colnames(mod), paste('peer_', 1:(ncol(modpeer) - ncol(mod)), sep=""))
    rownames(modpeer) = rownames(mod)
    write.matrix(modpeer, file=paste0(prefix, "peer.", n_factors, ".mod.txt"))

    colnames(wgt_all) = colnames(modpeer)
    rownames(wgt_all) = probes

    # regres out only the peer columns.
    peer_cols = (ncol(mod) + 1):(ncol(modpeer))
    peer_factors = modpeer[,peer_cols]
    peer_wgts = wgt_all[,peer_cols]


    write.matrix(wgt_all, file=paste0(prefix, "peer.", n_factors, ".weights.txt"))

    cleanY = t(data_complete - (peer_factors %*% t(peer_wgts)))
    if(any(is.na(cleanY))){
        message("NAs in cleanY")
    }
    rownames(cleanY) = probes
    colnames(cleanY) = rownames(mod)
    write.matrix(cleanY, file=paste(prefix, "cleaned.", n_factors, ".y.txt", sep=""))
    cleanY
}



limma.ez = function(data, mod, coef, contrasts, prefix, probe_len=1, genome_control=FALSE){
    fit = lmFit(data, mod)
    if(is.null(contrasts)){
        fit = eBayes(fit, trend=FALSE)
    } else {
        if(!class(contrasts) == "matrix"){
            err.log("levels", colnames(mod))
            contr.matrix = makeContrasts(contrasts=contrasts, levels=mod)
        } else {
            contr.matrix = contrasts
        }
        fit = eBayes(contrasts.fit(fit, contr.matrix), trend=FALSE)
        coef = colnames(fit)
    }
    err.log("variables in fit", paste(colnames(fit), collapse=" "))
    for(c in coef){
        if(is.numeric(c)){
            c = colnames(fit)[c]
        }
        err.log("getting p-values for:", c)
        tt = topTable(fit, coef=c, n=Inf, sort.by="none")

        # calculate unmoderated t-statistic
        ta = fit$coef[,c] / fit$stdev.unscaled[,c] / fit$sigma
        tt$unmoderated.t = ta
        tt$unmoderated.p = 2 * pt(-abs(ta), df=fit$df.residual[1])

        lambda = genomic_control(tt$P.Value)
        message(sprintf("genomic control value: %.2f", lambda))
        if(genome_control){
            if(lambda > 1.01){
                tt$P.Value = gc.adjust(tt$P.Value)
            }
            tt$adj.P.Val = p.adjust(tt$P.Value, "fdr")
        }
        chroms = unlist(lapply(strsplit(as.character(rownames(tt)), ":", fixed=TRUE), function(r){ r[1] }))
        starts = unlist(lapply(strsplit(as.character(rownames(tt)), ":", fixed=TRUE), function(r){ as.numeric(r[2]) - 1 }))
        if(!any(is.na(starts))){
            ends = starts + probe_len
            err.log(length(chroms))
            err.log(length(starts))
            err.log(colnames(tt))

            tt = data.frame(chrom=chroms, start=starts, end=ends, tt[,colnames(tt) != "ID"])
            tt = tt[order(chroms, starts),]
        }

        write.table(tt, row.names=F, sep="\t", quote=F, file=paste(prefix, ".", c, ".pvals.bed", sep=""))
    }
    return(fit)
}


.fix_anno = function(anno){
    probe = anno$probe
    info = anno$info
    anno = format(anno[,!colnames(anno) %in% c("probe", "info")], digits=4, trim=TRUE)
    anno = cbind(anno, probe, info)
    return(anno)
}

annotate_top_table = function(tt, probe_info="probe_lookups.txt"){
    lookup = read.table(probe_info, header=T,  sep="\t", quote="")
    anno = .fix_anno(merge(lookup, tt, by.y="ID", by.x="probe"))
    return(anno)
}

.adjust_prefix = function(prefix){
    if(substr(prefix, nchar(prefix), nchar(prefix)) %in% c(".", "-", "/")){
        return(prefix)
    }
    return(paste(prefix, ".", sep=""))
}

.shuffle_clinical = function(mod, seed, prefix){
    if(is.na(seed)){
        return(c(mod, prefix))
    } 
    rnm = rownames(mod)
    seed = as.integer(seed)
    set.seed(seed)
    perm = sample(nrow(mod))  
    if(seed > 0){
        err.log("shuffling clinical data, only last column")
        mod[, ncol(mod)] = mod[perm, ncol(mod)]
        prefix = paste(prefix, "shuffle.", colnames(mod)[ncol(mod)], ".", seed, ".", sep="")
    }
    else {
        err.log("shuffling clinical data, all columns")
        mod[, ] = mod[perm, ]
        prefix = paste(prefix, "shuffle.all", seed, ".", sep="")
    }
    rownames(mod) = rnm
    return(c(mod, prefix))
}

shuffle_matrix = function(mat, seed, dim=c("col", "row")){
    seed = as.integer(seed)
    set.seed(seed)
    if(dim == "col"){
        cn = colnames(mat)
        mat = mat[,sample(ncol(mat))]
        colnames(mat) = cn
        return(mat)
    }
    mat = mat[sample(nrow(mat)),]
    return(mat)
}

matrix.eQTL.post = function(prefix, expr_locs, marker_locs=NULL, anno=NULL){
    tra = read.tab(paste(prefix, 'eQTL_tra.txt', sep=""))
    cis = read.tab(paste(prefix, 'eQTL_cis.txt', sep=""))
    library(ChIPpeakAnno)
    
    if(is.null(marker_locs)){
        marker_locs = get_marker_locs(union(tra$SNP, cis$SNP))
        marker_locs$chromStart = marker_locs$pos - 1
        marker_locs$chromEnd = marker_locs$pos
        marker_locs$name = marker_locs$snp
        marker_locs = marker_locs[,c("chrom", "chromStart", "chromEnd", "name")]
    }
    # output: snp_chrom snp_start snp_end expr_chrom expr_start expr_end
    # cis_trans dist t-stat pval FDR expr_gene snp_stuff_from_chippeakanno.

    rg = BED2RangedData(marker_locs)
    anno = annotatePeakInBatch(rg, 
                  AnnotationData=anno, select="all")
    library(org.Hs.eg.db)
    # TODO: make this changeable.
    anno = addGeneIDs(anno, "org.Hs.eg.db", c("symbol"))

    eg = getEnrichedGO(anno, orgAnn="org.Hs.eg.db", maxP=0.01, multiAdj=TRUE,
                       multiAdjMethod="BH" )
    go_f = paste0(prefix, "snp.go.enrich.txt")
    write.table(rbind(eg$mf, eg$bp, eg$cc), row.names=FALSE, sep="\t", file=go_f)

    anno = as.data.frame(anno)
    colnames(anno)[1] == "chrom"
    anno[,1] = paste("chr", anno[,1], sep="")
    write.table(anno, file=paste0(prefix, "anno.txt"), sep="\t",
                quote=FALSE, row.names=FALSE)
    anno

}


get_marker_locs = function(names){
        chrm_snp = unlist(lapply(strsplit(as.character(names), ":", fixed=TRUE), 
                          function(r){ r[1] }))
        pos = unlist(lapply(strsplit(as.character(names), ":", fixed=TRUE),
                          function(r){ as.numeric(r[2]) }))

        snpspos = data.frame(snp=names, chrom=chrm_snp, pos=pos)
        #write.table(snpspos, row.names=T, sep="\t", quote=F)
        return(snpspos)
}


#library(ChIPpeakAnno)
#data(TSS.human.NCBI36)

matrix.eQTL.ez = function(expr_data, marker_data, clinical, model, prefix,
                            expr_locs, marker_locs=NULL,
                            cis_dist=1e6, seed=NA, gseed=NA,
                            linear_cross=FALSE,
                            p_thresh=1.0,
                            snp_size=1,
                            anno=TSS.human.NCBI36){
   
    prefix = .adjust_prefix(prefix)

    stopifnot(all(colnames(marker_data) == colnames(expr_data)))
    stopifnot(all(rownames(clinical) == colnames(expr_data)))


    full_formula = as.formula(model) 
    mod = as.matrix(model.matrix(full_formula, data=clinical))
    cnames = colnames(mod)
    rnames = rownames(mod)
    n = sum(!cnames %in% "(Intercept)")
    mod = matrix(mod[,!cnames %in% "(Intercept)"], nrow=nrow(mod), ncol=n)
    colnames(mod) = cnames[!cnames %in% "(Intercept)"]
    rownames(mod) = rnames

    complete = complete.cases(clinical[,attr(terms(full_formula), "term.labels")])
    err.log("removing:", sum(!complete), "because of missing data")
    err.log("leaving:", sum(complete), "rows of data.")

    marker_complete = as.matrix(marker_data[,complete])
    rm(marker_data); gc()

    library(MatrixEQTL)
    # TODO: just shuffle the locs? or just the marker locs?
    if(!is.na(gseed)){
        # shuffle genotype columns.
        gseed = as.integer(gseed)
        # avoid copy by not calling function.

        set.seed(gseed)
        cn = colnames(marker_complete)
        marker_complete = marker_complete[,sample(ncol(marker_complete))]
        colnames(marker_complete) = cn

        if (gseed < 0){
            err.log("shuffling genotype with clinical")
            mod = shuffle_matrix(mod, gseed, dim="col");
            prefix = paste(prefix, "shuffle.genotype_and_clin", ".", gseed, ".", sep="")
        } else {
            err.log("shuffling genotype columns")
            prefix = paste(prefix, "shuffle.genotype", ".", gseed, ".", sep="")
        }
    }

    expr_complete = SlicedData$new(as.matrix(expr_data[,complete]));
    expr_complete$ResliceCombined()
    rm(expr_data); gc(TRUE)

    marker_complete = SlicedData$new(marker_complete);
    marker_complete$ResliceCombined()

    err.log("expr using:", paste(dim(expr_complete), collapse=", ", sep=", "))
    err.log("snps using:", paste(dim(marker_complete), collapse=", ", sep=", "))
    if(linear_cross){
        err.log("using column:", colnames(mod)[ncol(mod)])
    }

    if(!is.na(seed)){
        mod_prefix = .shuffle_clinical(mod, seed, prefix)
        mod = as.matrix(mod_prefix[1])
        prefix = mod_prefix[2]
    } 
    if(ncol(mod) != 0){
        write.matrix(mod, name="ID", file=paste(prefix, "model.txt", sep=""))
    }
    clin = SlicedData$new(t(mod));

    # get the location of the SNPs.
    if(is.null(marker_locs)){
        snpspos = get_marker_locs(rownames(marker_complete));
    } else {
        snpspos = marker_locs
        rm(marker_locs);
    }
    gc()
    if(nrow(clin) != 0){
        stopifnot(all(colnames(marker_complete) == colnames(clin)))
    }
    stopifnot(all(colnames(expr_complete) == colnames(marker_complete)))
    stopifnot(nrow(snpspos) == nrow(marker_complete))

    genepos = data.frame(geneid=expr_locs$probe, chrm_probe=expr_locs$chrom,
                         start=expr_locs$start, end=expr_locs$end)
    rm(expr_locs)
    gc(TRUE)

    output_file_name_tra = paste(prefix, 'eQTL_tra.txt', sep="")
    output_file_name_cis = paste(prefix, 'eQTL_cis.txt', sep="")

    while(p_thresh * marker_complete$nRows() * expr_complete$nRows() > 1000000){
        p_thresh = p_thresh / 10.0
        # set to the maximum allowable p_threshold by MatrixeQTL must be a
        # multipe of 10.
    }
    err.log("p-value p_threshold:", p_thresh)

    me = Matrix_eQTL_main(
        snps = marker_complete,
        gene = expr_complete,
        cvrt = clin,
        output_file_name  = output_file_name_tra,
        pvOutputThreshold = p_thresh,
        #pvOutputThreshold = 0,
        # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/manual.html#models
        useModel = ifelse(linear_cross, modelLINEAR_CROSS, modelLINEAR),
        errorCovariance = numeric(),
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = p_thresh,
        #pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cis_dist,
        pvalue.hist = "qqplot");

    jpeg(paste(prefix, "qq.plot.jpeg", sep=""))
    plot(me)
    dev.off()
    err.log("cis tests:", me$cis$ntests)
    err.log("trans tests:", me$trans$ntests)
    #err.log(summary(me))
    me$prefix = prefix

    out_file = .add_dist(output_file_name_cis, output_file_name_tra, snpspos, genepos, prefix, snp_size)
    #plot_enrichment(output_file_name_cis, output_file_name_tra, me$cis$ntests,
    #           me$trans$ntests, paste0(prefix, "cis-enrichment.png"))
    if(!all(is.na(anno))){
        anno = matrix.eQTL.post(prefix, expr_locs, anno=anno)
    }
    return(me)
}

plot_enrichment = function(fcis, ftra, ncis, ntra, out_name,
   column="FDR"){

    png(out_name)
    cis = read.delim(fcis)
    tra = read.delim(ftra)

    p.max = max(cis[,column], tra[,column])
    p.min = min(cis[,column], tra[,column])

    p.rng = 10^-seq(-log10(p.max), -log10(p.min), length.out=51)

    cis_counts = rep(NA, 20)
    tra_counts = rep(NA, 20)
    chi.ps = rep(NA, 20)
    i = 0
    for(pcutoff in p.rng){
        i = i + 1
        cis_counts[i] = sum(cis[,column] < pcutoff)
        tra_counts[i] = sum(tra[,column] < pcutoff)
        chi.ps[i] = max(1e-300, chisq.test(matrix(c(cis_counts[i], ncis,
                                     tra_counts[i], ntra), nrow=2))$p.value)

        if(tra_counts[i] == 0 || cis_counts[i] == 0){ break }
    }
    cis_counts = cis_counts[1:i]
    tra_counts = tra_counts[1:i]
    if(length(cis_counts) == 1){ return }
    chi.ps = chi.ps[1:i]
    cis_enrichment = (cis_counts / ncis) / (tra_counts / ntra)
    ymax = max(-log10(chi.ps), cis_enrichment) * 1.15
    plot(-log10(p.rng[1:i]), 
        cis_enrichment,
        type='b',
        pch=5,
        ylim=c(0, ymax),
        ylab="", xlab=sprintf("-log10(%s-cutoff)", column),
        main="cis enrichment of QTLs", col="black")
    points(-log10(p.rng[1:i]), -log10(chi.ps), pch=19, col="blue", type='b')

    legend(-log10(p.rng[i]), ymax, c("cis enrichment", 
        "-log10(chisq-p) of enrichment"), col=c("black", "blue"),
        pch=c(5, 19), xjust=1, bg="transparent")

    dev.off()
}

mart_anno = function(peakList, dataset="hsapiens_gene_ensembl", 
                                      host="may2009.archive.ensembl.org"){
    mart = useMart(biomart="ensembl", dataset=dataset, host=host)
    TSS = getAnnotation(mart, featureType="TSS")
    utr5 = getAnnotation(mart, featureType="5utr")
    utr3 = getAnnotation(mart, featureType="3utr")
    exon = getAnnotation(mart, featureType="Exon")
    assignChromosomeRegion(myPeakList, exon, TSS, utr5, utr3)
}

.add_dist = function(fcis, ftrans, snpspos, genepos, prefix, snp_size=0){
    # snp_size accounts for when the snp is actually a probe...
    gene_chrom = paste(colnames(genepos)[2], "_gene", sep="")
    gene_start = paste(colnames(genepos)[3], "_gene", sep="")
    gene_end = paste(colnames(genepos)[4], "_gene", sep="")
    snp_chrom = paste(colnames(snpspos)[2], "_snp", sep="")

    qtls = read.tab(fcis)
    out_names = colnames(qtls)

    qtls$cis_tra = "cis"

    trans = read.tab(ftrans)
    trans$cis_tra = "tra"

    qtls = rbind(qtls, trans)
    rm(trans); gc()

    output_file = paste(prefix, 'eQTL_ct_dist.txt', sep="")
    # SNP gene t-stat p-value FDR
    colnames(snpspos) = paste(colnames(snpspos), "_snp", sep="")
    colnames(genepos) = paste(colnames(genepos), "_gene", sep="")

    all_dat = merge(qtls, snpspos, by.x=1, by.y=1)
    all_dat = merge(all_dat, genepos, by.x=2, by.y=1)

    dist_s = all_dat$pos_snp - all_dat[,gene_start]
    dist_e = all_dat$pos_snp - all_dat[,gene_end]
    dist_s1 = (all_dat$pos_snp + snp_size) - all_dat[,gene_start]
    dist_e1 = (all_dat$pos_snp + snp_size) - all_dat[,gene_end]

    all_dat$dist = apply(cbind(abs(dist_s), abs(dist_e), abs(dist_s1), abs(dist_e1)), 1, min)
    all_dat$dist[(dist_s > 0 | dist_s1 > 0) & (dist_e < 0 | dist_e1 < 0)] = 0 # account for snp inside gene. set dist to 0
    all_dat$dist[all_dat[,gene_chrom] != all_dat[,snp_chrom]] = NA
    all_dat = all_dat[,c("chrom_snp", "pos_snp", "pos_snp", out_names, "cis_tra", "dist")]

    all_dat[,3] = all_dat[,3] + snp_size
    colnames(all_dat)[1:3] = c("chrom", "start", "end")
    write.table(all_dat, file=output_file, row.names=FALSE, sep="\t", quote=FALSE)
    return(output_file)
}

read.agilent = function(targets, names=NULL, path=NULL){
  read.maimages(targets, source="agilent", green.only=TRUE, 
                    names=names, path=path)
}


normalize.agilent = function(x, offset=10){

  y = backgroundCorrect(x, method="normexp", offset=offset)
  y = normalizeBetweenArrays(y, method="quantile")
  
  # from the limma manual
  # the 95th percentile of the negative control probes on the array
  neg95 = apply(y$E[y$genes$ControlType==-1,], 2, function(x) quantile(x, p=0.95))
  
  cutoff = matrix(1.1 * neg95, nrow(y), ncol(y), byrow=TRUE)
  isexpr = rowSums(y$E > cutoff) >= (ncol(x) / 3)
  
  expr = y[y$genes$ControlType==0 & isexpr,]
  expr.ave = avereps(expr, ID=expr$genes[,"ProbeName"])
  expr.ave
}

agilent.limma = function(targets, model, names=NULL, coef=2,
                         offset=10, path=NULL){

  mm = model.matrix(model, targets)
  x = read.agilent(targets, names, path=path)

  y = normalize.agilent(x, offset=offset)
  
  topTable(eBayes(lmFit(y, mm), trend=TRUE), coef=coef, n=Inf, adjust.method="fdr")
  
}

glht.fit.ez = function(dat, clin, model_str, comparison, mc.cores=4){
  # this is used for fitting lme4 functions, where a model is, e.g.
  # ~ 0 + disease + age + (1|family)
  # with comparison of 'diseaseCOPD - diseaseIPF = 0' as would be
  # specified to multcomp::glht
  # this function parallelizes that and returns the pvalue, coefficient, and
  # t-statistic
  library(lme4)
  library(multcomp)
  library(parallel)
  model = gsub("^\\s+", "", as.character(model_str)) # remove initial whitespace
  if(substring(model, 1, 1) == "~"){
    model = paste0("y ", model)
  } else { 
    if(!substring(model, 1, 1) == "y"){
        stop(paste("model should start with '~':", model))
    }
  }

  res = mclapply(1:nrow(dat), function(i){
    if(i %% 10000 == 0){ message(paste("at record", i)) }
    y = dat[i,]
    mod = lmer(as.formula(model), clin)    
    r = glht.fit.one(y, mod, comparison)
    r$probe = rownames(dat)[i]
    r$cmp = rownames(r)
    rownames(r) = NULL
    r
  }, mc.cores=mc.cores)
  res = rbindlist(res)
  res$qvalue = p.adjust(res$pvalue, "fdr")
  res 
}

glht.fit.one = function(y, mod, comparison){
  s = summary(glht(mod, linfct=comparison))
  data.frame(pvalue=s$test$pvalues,
             t=s$test$tstat,
             coefficient=s$test$coefficients)
}

