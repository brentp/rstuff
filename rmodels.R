library(limma)
library(gtools)

options(scipen=15, stringsAsFactors=FALSE) # stop it from printing 1e6 instead of 1000000

logit = defmacro(p, expr=log(p) - log(1 - p))

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
    s$fileSliceSize = 18000;
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

normalize.charm = function(sample_description, out_prefix, id_col=1){
    library(charm)
    require(BSgenome.Hsapiens.UCSC.hg18)
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

    ctrlIdx = getControlIndex(rawData, subject=Hsapiens, noCpGWindow=600)
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

remove_low_variance = function(mat, p_drop=0.25){
    if(p_drop > 1){ p_drop = p_drop / 100 }
    rvm = rowVars(mat)
    mat[rvm > sort(rvm)[p_drop * length(rvm)],]
}


write.matrix = function(mat, file, name="probe", quote=FALSE, sep="\t", digits=3, ...){
    mat = cbind(rownames(mat), format(round(mat, digits=digits), digits=digits, trim=TRUE))
    colnames(mat)[1] = name
    write.table(mat, file=file, quote=quote, sep=sep, row.names=FALSE, ...)
    gc();
}    


peer.limma.ez = function(data, clin, model,
                prefix="results/peer.limma.",
                contrasts=NULL, probe_len=1,
                batch_correct=5){
    #https://github.com/PMBio/peer/wiki/Tutorial
    coef = 2

    # the first term is always the variable of interest.
    full_formula = as.formula(model)
    err.log("full model:", full_formula)

    # remove any rows with na
    complete = complete.cases(clin[,attr(terms(full_formula), "term.labels")])
    err.log("removing:", sum(!complete), "because of missing data")
    err.log("leaving:", sum(complete), "rows of data.")
    data_complete = as.matrix(data[,complete])
    rm(data); gc()

    mod  = model.matrix(full_formula, data=clin)


    if(!is.na(batch_correct) && as.logical(batch_correct)){
        n_factors = as.integer(batch_correct)
        peer_factors = run.peer(mod, t(data_complete), prefix, n_factors)
        modpeer = cbind(mod, peer_factors)
        write.table(modpeer, sep="\t", row.names=F, file=paste(prefix, "mod.peer.", n_factors, ".txt", sep=""), quote=F)

        # we remove the batch effects from the data, rather than including the
        # peers as covariates.
        data_complete = cleanY(data_complete, mod, peer_factors)
        write.matrix(data_complete, file=paste(prefix, "cleaned.", n_factors, ".y.txt", sep=""))
    }
    # including all peer factors, just include those before 1/alpha levels off?
    fit = limma.ez(data_complete, mod, coef, contrasts, prefix, probe_len)
    return(fit)
}


run.peer = function(mod, data_complete, prefix=NULL, n_factors=5){
    library(peer)
    peer_obj = PEER()
    PEER_setCovariates(peer_obj, mod)
    PEER_setTolerance(peer_obj, 1e-10)
    PEER_setVarTolerance(peer_obj, 1e-10)
    PEER_setNmax_iterations(peer_obj, 1000)
    PEER_setPhenoMean(peer_obj, data_complete)
    nK = min(n_factors, nrow(data_complete) - ncol(mod) - 1)
    PEER_setNk(peer_obj, nK)
    err.log("set number of factors to infer as", nK)
    PEER_update(peer_obj)
    modpeer = PEER_getX(peer_obj)
    if(!is.null(prefix)){
        f = paste(prefix, "peer.alpha.txt", sep="")
        write.table(PEER_getAlpha(peer_obj), file=f, row.names=F, sep="\t")
    }
    colnames(modpeer) = c(colnames(mod), paste('peer_', 1:(ncol(modpeer) - ncol(mod)), sep=""))
    # return only the peer columns.
    return(modpeer[,(ncol(mod) + 1):(ncol(modpeer))])
}



limma.ez = function(data, mod, coef, contrasts, prefix, probe_len=1){
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

        chroms = unlist(lapply(strsplit(as.character(tt$ID), ":", fixed=TRUE), function(r){ r[1] }))
        starts = unlist(lapply(strsplit(as.character(tt$ID), ":", fixed=TRUE), function(r){ as.numeric(r[2]) - 1 }))
        ends = starts + probe_len
        err.log(length(chroms))
        err.log(length(starts))
        err.log(colnames(tt))

        tt = data.frame(chrom=chroms, start=starts, end=ends, tt[,colnames(tt) != "ID"])

        write.table(tt, row.names=F, sep="\t", quote=F, file=paste(prefix, c, ".pvals.bed", sep=""))
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

matrix.eQTL.post = function(prefix, expr_locs, marker_locs=NULL, anno_lib="BSgenome.Hsapiens.UCSC.hg18"){
    tra = read.tab(paste(prefix, 'eQTL_tra.txt', sep=""))
    cis = read.tab(paste(prefix, 'eQTL_cis.txt', sep=""))
    
    if(is.null(marker_locs)){
        marker_locs = get_marker_locs(intersect(tra$SNP, cis$SNP))
        marker_locs$chromStart = marker_locs$pos - 1
        marker_locs$chromEnd = marker_locs$pos
        marker_locs$name = marker_locs$SNP
    }
    # output: snp_chrom snp_start snp_end expr_chrom expr_start expr_end
    # cis_trans dist t-stat pval FDR expr_gene snp_stuff_from_chippeakanno.
    library(ChIPpeakAnno)
    library(anno_lib)
    library(rtracklayer)

    rg = BED2RangedData(marker_locs)
    anno = as.data.frame(annotatePeakInBatch(rg, Annotation=anno_lib, select="all"))





# TODO: write separate annotate fn.

    # http://www.stat.berkeley.edu/share/biolab/Courses/CIPF10/Data/lab3_solutions.r


    #test.bed = data.frame(cbind(chrom = c("4", "6"),
    #chromStart=c("100", "1000"),chromEnd=c("200", "1100"),
    #name=c("peak1", "peak2")))
    #test.rangedData = BED2RangedData(test.bed)
    #as.data.frame(annotatePeakInBatch(test.rangedData,
    #BED2RangedData() 

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



matrix.eQTL.ez = function(expr_data, marker_data, clinical, model, prefix,
                            expr_locs, marker_locs=NULL,
                            cis_dist=1e6, seed=NA, gseed=NA,
                            linear_cross=FALSE,
                            p_thresh=1.0,
                            snp_size=1){
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
        mod = mod_prefix[1]
        prefix = mod_prefix[2]
    } 
    if(ncol(mod) != 0){
        write.matrix(mod, name="ID", file=paste(prefix, "model.txt", sep=""))
    }
    clin = SlicedData$new(t(mod));

    # get the location of the SNPs.
    if(is.null(marker_locs)){
        snpspos = get_marker_locs(rownames(marker_complete));
    }
    else {
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

    .add_dist(output_file_name_cis, output_file_name_tra, snpspos, genepos, prefix, snp_size)
    return(me)
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

    all_dat = all_dat[,c(out_names, "cis_tra", "dist")]
    write.table(all_dat, file=output_file, row.names=FALSE, sep="\t", quote=FALSE)
    return(output_file)
}


freedman_lane_permute = function(y, model_matrix, cols){
   # cols are the columns of interest in the design matrix.
   # if only a single column is specified, the t-stat is used. otherwise
   # the F-statistic is used.
   # Returns a fit object with the simulated p-values using
   #  the Freedman-Lane method
   # of shuffling residuals from the reduced model.
   # TODO:
   # if use_beta = True, then it utilizes irrizary et al's bump hunting
   proportion = 0.02
   design = model_matrix
   fit = eBayes(lmFit(y, design), proportion=proportion)

   reduced_design = design[,!colnames(design) %in% cols]
   reduced_fit = eBayes(lmFit(y, reduced_design))
   reduced_resid = residuals(reduced_fit, y)
   reduced_fitted = fitted(reduced_fit, y)
   rm(reduced_fit); gc()

   stat_col = ifelse(length(cols) > 1, "F", "t")
   tt = topTable(fit, coef=cols, n=Inf, sort.by="none")
   rm(fit); gc()

   stat_orig = abs(tt[,stat_col])
   
   n_cols = ncol(reduced_resid)
   n_greater = rep(0, nrow(y))
   n_perms = rep(0, nrow(y))
   g_subset = rep(TRUE, nrow(y))
   cutoff = 0.20
   # THIS sections calls the simulation on shuffled data. after each loop.
   # it takes only the subset that has a perm_p below some less stringent cutoff
   # so it does not waste time retesting probes that have a high p-value after 25
   # sims.
   for (n_perm in c(25, 80, 240, 720, 1650, 5000)){
           print(paste(cutoff, n_perm))
           n_greater[g_subset] = .freedman_lane_sim(reduced_fitted[g_subset,],
                                                    reduced_resid[g_subset,],
                                                    design,
                                                    cols,
                                                    n_greater[g_subset],
                                                    n_perm,
                                                    stat_orig[g_subset],
                                                    proportion)
           n_perms[g_subset] = n_perms[g_subset] + n_perm
           if((sum(n_greater[g_subset] < cutoff)) == 0){ break }
	   g_subset = g_subset & (n_greater < cutoff)
           proportion = proportion * 2.0
           cutoff = cutoff / 2
         
   }
   sim_p = as.matrix((1 + n_greater) / (1 + n_perms), ncol=1)
   return(cbind(tt, sim_p))
}

.freedman_lane_sim = function(reduced_fitted, reduced_resid, design, cols, n_greater, n_perms, stat_orig, proportion){
   # number of simulations with a stat greater than the observed.
   nc = ncol(reduced_resid)
   stat_col = ifelse(length(cols) > 1, "F", "t")
   print(dim(reduced_resid)) 
   
   for(i in 1:n_perms){
      ystar = reduced_fitted + reduced_resid[, sample(1:nc)]
      fit_sim = eBayes(lmFit(ystar, design), proportion=proportion)
      tt_sim = topTable(fit_sim, coef=cols, n=Inf, sort.by="none")
      n_greater = n_greater + (abs(tt_sim[,stat_col]) > stat_orig)
   }
   return(n_greater)
}
