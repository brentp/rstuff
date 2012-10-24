ff = function(r0, r1, df0, df1, p=FALSE){
    n = ncol(r0)
    stopifnot(n == ncol(r1))
    rss0 = rowSums(r0 * r0) 
    rss1 = rowSums(r1 * r1)
    fstats = ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
    if(!p){ 
        return(fstats)
    }
    p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
    return(data.frame(fstats=fstats, p=p))
}


# This is a modified implementation of the freedman-lane algorithm described in:
# http://www.ncbi.nlm.nih.gov/pubmed/17630650 by Wagner et al.
# where one permutes the residuals of a reduced model in order to generate p-values.
# It also implements modification of the adusted p-value described in that paper.
# The additions implemented here are:
# 1. significance by comparing beta values (in addition to F stats as above)
# 2. p-value adjustment is actually a q-value, not a fixed level.
# 3. optimized matrix-based calculation for speed
# 4. very low memory-use. itermediate values are not stored.
# 5. because it is very computationally intensive to perform, e.g. 100K
#    simulations on 40K probes, this algorithm is performed iteratively.
#    For example after only 40 iterations on the full data set, the function
#    will then select only those probes with a current simulated p-value less
#    than 0.2. At that time it will perform 80 iterations on an expected 8K
#    probes. Then 160 iterations on 4K. And so on, so that only the probes
#    in which we are most interested (those with the lowest p-values) undergo
#    the highest number of iterations. In this way, the probes with the lowest
#    p-values may undergo 100K iterations so that their simulated p and q values
#    are both precise and accurate.
# The output is such that this can be parallelized and combined. However, on
# a modest laptop. This function runs a set of 12,042 genes through as many
# as 163K iterations in just over 6 minutes.
freedman_lane_permute = function(y, mod, mod0, use_beta=TRUE, max_iter=100000){
   # if use_beta = True, then it runs compares the simulated parameter
   # estimate to the observed. oterhwise, it compares f-statistics.
   design = mod
   reduced_design = mod0

   fit0 = lm.fit(mod0, t(y))
   fit = lm.fit(mod, t(y))
   reduced_fitted = t(fit0$fitted.values)
   reduced_resid = t(fit0$residuals)
   asym = ff(reduced_resid, t(fit$residuals), ncol(mod0), ncol(mod), p=TRUE)

   if(use_beta){
           use_beta = setdiff(colnames(mod), colnames(mod0))
           if(length(use_beta) > 1){
               stop("can only have a single covariate if use_beta is True")
           }
           stat_orig = abs(fit$coef[use_beta,])
           # actually not reduced residuals, but leaving name for now...
           # use residuals from full model and fit from null model.
           asym$beta = stat_orig
   } else {
       stat_orig = asym$fstats
   }
   rm(fit0, fit); gc()
   
   n_greater = rep(0, nrow(y))
   qvals = rep(0, nrow(y))
   n_perms = rep(0, nrow(y))
   g_subset = rep(TRUE, nrow(y))
   cutoff = 0.2
 
   n_perm = 40
   # THIS sections calls the simulation on shuffled data. after each loop.
   # it takes only the subset that has a perm_p below some less stringent cutoff
   # so it does not waste time retesting probes that have a high p-value after 25
   # sims.
   n_samples_full = nrow(y)

   print_subset = order(asym$p)[1:20]

   loops = rep(0, nrow(y))
   while(sum(g_subset) > 0 && cutoff > 1 / nrow(reduced_resid) && max(n_perms) < max_iter){
       write(sprintf("performing %i shufflings of %i rows then limiting to < %.4g",n_perm, sum(g_subset), cutoff), stderr())
       t = Sys.time()

       li = .freedman_lane_sim(reduced_fitted[g_subset,],
                                                reduced_resid[g_subset,],
                                                design,
                                                reduced_design,
                                                n_greater[g_subset],
                                                n_perm,
                                                stat_orig[g_subset],
                                                use_beta,
                                                n_samples_full)

       # don't add n_greater because it's added in sim func.
       n_greater[g_subset] = li[["n_greater"]]
       qvals[g_subset] = qvals[g_subset] + li[["qvals"]]
       loops[g_subset] = loops[g_subset] + 1
       n_perms[g_subset] = n_perms[g_subset] + n_perm

       # this line will show the evolution of the 20 probes with the lowest
       # asymptotic p-value
       #print(qvals[print_subset] / loops[print_subset])

       # n_great / n_perms is the simulated p-value
       # for the next iteration, we only care about probes that
       # have a current simulated p-value less than the cutoff
       g_subset = g_subset & (((n_greater / n_perms) < cutoff) ) # | (asym$p[g_subset] < cutoff))
       # the cutoff becomes more stringent each time.
       cutoff = cutoff / 2
       n_perm = n_perm * 2
       write(Sys.time() - t, stderr())
   }
   # qvals is the average across all loops.
   qvals = pmin(1, qvals / loops)

   sim_p = n_greater / n_perms
   ret = data.frame(
          probes=rownames(dat),
          sim_p,
          n_greater=n_greater,
          n_perms=n_perms,
          asym_p=format(asym$p, digits=4, nsmall=3, trim=TRUE),
          F=format(asym$fstats, digits=4, nsmall=2, trim=TRUE),
          beta=format(asym$beta, digits=4, nsmall=3, trim=TRUE),
          sim_q=qvals
   )
   rownames(ret) = rownames(dat)
   return(ret)
}



.freedman_lane_sim = function(reduced_fitted, reduced_resid, design, 
                              reduced_design, n_greater, n_perms,
                              stat_orig, use_beta, n_samples_full){
   # number of simulations with a stat greater than the observed.
   nc = ncol(reduced_resid)
   qvals = rep(0, nrow(reduced_resid))
   qdiv =  n_samples_full / (1:length(qvals))

   reduced_resid = t(reduced_resid)
   reduced_fitted = t(reduced_fitted)

   pb = txtProgressBar(min = 0, max = n_perms, style = 3)

   for(i in 1:n_perms){
      ystar = reduced_fitted + reduced_resid[sample(1:nc),]
      if(use_beta != "FALSE"){
          stat = abs(lm_fit_fast(ystar, design, use_beta))
          # we have changed use_beta to the name of the variable to extract.
      } else {
          #fit = lm.fit(design, ystar)
          #fit0 = lm.fit(reduced_design, ystar)
          #f = ff(t(fit0$residuals), t(fit$residuals),
          #       ncol(reduced_design), ncol(design), p=FALSE)
          fit0_res = lm_fit_fast(ystar, reduced_design, "res")
          fit_res = lm_fit_fast(ystar, design, "res")
          stat = abs(ff(t(fit0_res), t(fit_res),
                 ncol(reduced_design), ncol(design), p=FALSE))
      }
      n_greater = n_greater + (stat > stat_orig)
      qv = ecdf(-(stat))(-stat_orig)
      o = order(qv)
      ro = order(o)
      #qv = ((qv[o] * length(qv)) / (1:length(qv)))[ro]
      # correct by the number in the full sample...
      # n_samples_full is like the n argument to p.adjust()
      #qv = ((qv[o] * n_samples_full) / (1:length(qv)))[ro]
      qv = (qv[o] * qdiv)[ro]
      qvals = qvals + qv
      setTxtProgressBar(pb, i)

   }
   close(pb)
   return(list("n_greater"=n_greater, "qvals"=qvals / n_perms))
}

# from Douglas Bates' artical in R News! June 2004
# naive.sol <- solve(t(X) %*%  X) %*% t(X) %*% dat
# cpod.sol <- solve(crossprod(X), crossprod(X, dat))
# ch <- chol(crossprod(X))
#    chol.sol <- backsolve(ch, forwardsolve(ch, crossprod(X, dat), upper = TRUE, trans = TRUE))

lm_fit_fast <- function(dat, X, out) {
    # originally from R charm package.
    #sol = solve(t(X)%*%X)%*%t(X) %*% t(dat)
    beta = solve(crossprod(X), t(crossprod(dat, X)))
    if(out %in% rownames(beta)) return(beta[out,])

    if(out=="fit") {
        return(t(X%*% beta))
    } else if(out=="res") {
        return(dat - (X%*% beta))
    } else stop("invalid out arg")
}

