rstuff
======

R stuff I wrote

perm-p
======

This is a modified implementation of the freedman-lane algorithm described in:
http://www.ncbi.nlm.nih.gov/pubmed/17630650 by Wagner et al.
where one permutes the residuals of a reduced model in order to generate p-values.
It also implements modification of the adusted p-value described in that paper.
The additions implemented here are:

1. significance by comparing beta values (in addition to F stats as above)
2. p-value adjustment is actually a q-value, not a fixed level.
3. optimized matrix-based calculation for speed
4. very low memory-use. itermediate values are not stored.
5. because it is very computationally intensive to perform, e.g. 100K
   simulations on 40K probes, this algorithm is performed iteratively.
   For example after only 40 iterations on the full data set, the function
   will then select only those probes with a current simulated p-value less
   than 0.2. At that time it will perform 80 iterations on an expected 8K
   probes. Then 160 iterations on 4K. And so on, so that only the probes
   in which we are most interested (those with the lowest p-values) undergo
   the highest number of iterations. In this way, the probes with the lowest
   p-values may undergo 100K iterations so that their simulated p and q values
   are both precise and accurate.

The output is such that this can be parallelized and combined. However, on
a modest laptop. This function runs a set of 12,042 genes through as many
as 163K iterations in just over 6 minutes.

