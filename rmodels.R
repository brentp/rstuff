library(limma)
options(scipen=100) # stop it from printing 1e6 instead of 1000000


.write.raw = function(RGset, dp, out_prefix){
    raw = preprocessRaw(RGset)
    colnames(dp) = paste("detection-p-", colnames(dp), sep="")
    methylated = getMeth(raw)
    colnames(methylated) = paste("methylated-", colnames(methylated), sep="")
    unmethylated = getUnmeth(raw)
    colnames(unmethylated) = paste("unmethylated-", colnames(unmethylated), sep="")
    dpout = cbind(unmethylated, methylated, format(dp, digits=3))
    write.table(dpout, sep="\t", quote=FALSE, file=paste(out_prefix, "raw-values.txt", sep=""), row.names=T)
    rm(methylated); rm(unmethylated); gc();
}

read.tab = function(fname, sep="\t", header=TRUE, ...){
    if(class(fname) != "character"){
        return(fname)
    }
    clin = read.delim(fname, sep=sep, header=TRUE, ...)
    return(clin)
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
    .write.raw(RGset, dp, out_prefix)

    beta = getBeta(Mset.swan)
    colnames(beta) = pd[,id_col]
    message("writing beta ...")
    write.table(format(beta, digits=3), sep="\t", quote=FALSE, file=paste(out_prefix, "beta.txt", sep=""), row.names=T)
    rm(beta)
    gc()

    M = getM(Mset.swan)
    colnames(M) = pd[,id_col]

    message("writing M...")
    write.table(format(M, digits=3), sep="\t", quote=FALSE, file=paste(out_prefix, "M.txt", sep=""), row.names=T)
    M[dp > 0.05] = NA
    message("values > 0.05:", sum(dp > 0.05))
    write.table(format(M, digits=3), sep="\t", quote=FALSE, file=paste(out_prefix, "M.pgt05.txt", sep=""), row.names=T)
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
    qual = qcReport(rawData, file="qc-report.pdf")
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
    write.table(format(p, digits=3), sep="\t", quote=FALSE, file=paste(out_prefix, ".methp.txt", sep=""), row.names=T)
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
    library(sva)

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

    library(Matrix)
    mod  = model.matrix(full_formula, data=clin)


    if(batch_correct){
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
        tmp = cbind(rownames(data_complete), data_complete)
        colnames(tmp)[1] = "probe"
        write.table(tmp, row.names=F, sep="\t", quote=F, file=paste(prefix, ".cleaned.y.txt", sep=""))
        rm(tmp); gc();
    }

    fit = limma.ez(data_complete, mod, coef, contrasts, prefix, probe_len)
    return(fit)
}


peer.limma.ez = function(data, clin, model,
                prefix="results/peer.limma.",
                contrasts=NULL, probe_len=1,
                batch_correct=TRUE){
    #https://github.com/PMBio/peer/wiki/Tutorial
    library(peer)
    coef = 2

    # the first term is always the variable of interest.
    full_formula = as.formula(model)
    err.log("full model:", full_formula)

    # remove any rows with na
    complete = complete.cases(clin[,attr(terms(full_formula), "term.labels")])
    err.log("removing:", sum(!complete), "because of missing data")
    err.log("leaving:", sum(complete), "rows of data.")
    data_complete = as.matrix(data[,complete])

    mod  = model.matrix(full_formula, data=clin)

    if(batch_correct){
        peer_factors = run.peer(mod, t(data_complete), prefix)
        modpeer = cbind(mod, peer_factors)
        write.table(modpeer, sep="\t", row.names=F, file=paste(prefix, "mod.peer.txt", sep=""), quote=F)

        # we remove the batch effects from the data, rather than including the
        # peers as covariates.
        data_complete = cleanY(data_complete, mod, peer_factors)
        tmp = cbind(rownames(data_complete), data_complete)
        colnames(tmp)[1] = "probe"
        write.table(tmp, row.names=F, sep="\t", quote=F, file=paste(prefix, ".cleaned.y.txt", sep=""))
    }

    # including all peer factors, just include those before 1/alpha levels off?
    fit = limma.ez(data_complete, mod, coef, contrasts, prefix, probe_len)
    return(fit)
}


run.peer = function(mod, data_complete, prefix=NULL){
    peer_obj = PEER()
    PEER_setCovariates(peer_obj, mod)
    PEER_setTolerance(peer_obj, 1e-9)
    PEER_setVarTolerance(peer_obj, 1e-9)
    PEER_setNmax_iterations(peer_obj, 1000)
    PEER_setPhenoMean(peer_obj, data_complete)
    PEER_setNk(peer_obj, min(5, nrow(data_complete) - ncol(mod) - 1))
    PEER_update(peer_obj)
    modpeer = PEER_getX(peer_obj)
    if(!is.null(prefix)){
        f = paste(prefix, "peer.alpha.txt", sep=".")
        write.table(PEER_getAlpha(peer_obj), file=f, row.names=F, sep="\t")
    }
    print(dim(modpeer))
    print(dim(mod))
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