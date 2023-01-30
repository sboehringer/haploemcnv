Npairs = function(N)( (N^2 + N)/2 )
# N^2 + N - 2*Np = 0
Npairs2N = function(Npairs)floor( -1/2 + sqrt(1/4 + 2*Npairs) )
# alleles: 0, 1, ...
gtFromPairS = function(Nalleles, a)return(Npairs(a[2]) + a[1]);
gtFromPair = function(Nalleles, a)gtFromPairS(Nalleles, sort(a));
# gts: 0, 1, 2, ...
pairFromGt = function(Nalleles, gt) {
  mx = Npairs2N(gt);
  return(c(gt - Npairs(mx), mx));
}

dtapply = function(dt, f, ...)apply(matrix(dt, ncol = 2, byrow = T), 1, f, ...)
dtsapply = function(dts, f, ..., iter = lapply)iter(dts, dtapply, f = f, ...)
gtAf = function(g)(mean(g)/2);
gtsAfs = function(gts)apply(gts, 2, gtAf);

#
#	<p> EM algorithm
#

vector.assign1p = function(v, idcs, ...)vector.assign(v, idcs + 1, ...)
# to which haplotypes does a given diplotype contribute?
dtContribMat = function(Nhts) {
  Ndts = Npairs(Nhts);
  templ = rep(0, Nhts);
  dt2ht = apply(sapply(0:(Ndts - 1), pairFromGt, Nalleles = Nhts), 2, vector.assign1p, e = 1, v = templ);
  dt2ht2 = t(dt2ht) * 2 / apply(dt2ht, 2, sum);	# normalized to two alleles
  return(dt2ht2);
}
# map diplotype indices to indeces of Nhts x Nhts diplotype frequency matrix
dtFreqExtraction = function(Nhts) {
  combs = merge.multi.list(list(list(1:Nhts), list(1:Nhts)));
  combs = combs[combs[, 1] <= combs[, 2], ];
  gts = apply(combs, 1, gtFromPair, Nalleles = Nhts);	
  return(combs[order(gts), ]);
}
# express diplotype frequncies in a Nhts x Nhts matrix, assuming HWE
htfs2dtfs = function(htfs, dtfI) {
  dtfsRaw = htfs %*% t(htfs);
  dtfs = array.extract(2*dtfsRaw - diag(diag(dtfsRaw)), dtfI[, 1], dtfI[, 2]);
  return(dtfs);
}
# standardize matrix by column
colStd = function(m)t(t(m) / apply(m, 2, sum))
rowStd = function(m)t(colStd(t(m)))

vector.assign = function(v, idcs, e, na.rm = 0, N) {
  if (!missing(N)) v = rep(v, N);
  v[idcs] = e;
  if (!is.na(na.rm)) v[is.na(v)] = na.rm;
  v
}

# indicator matrix of possible diplotypes (genotypes), N x |{gts}| matrix
alleleFreqEstEM = function(gts, eps = 1e-5, Nitmax = 1e2) {
  Nhts = Npairs2N(ncol(gts));
  Mdt2ht = dtContribMat(Nhts);	# contribution of diplotypes to haplotypes
  dtfI = dtFreqExtraction(Nhts);	# indices in frequency matrix
  afs = cbind(1 - gtsAfs(gts), gtsAfs(gts));
  # initialize frequency vector
  htfs = vn(rep(1, Nhts));
  
  for (i in 1:Nitmax) {
    dtfs = htfs2dtfs(htfs, dtfI);
    dtCond = colStd(t(gts) * dtfs);	# conditional probabilities
    htfsN = vn(apply(t(dtCond) %*% Mdt2ht, 2, sum));
    LogS(5, "%{i}03d [%{dev}.3f]: %{freqs}s",
         dev = max(abs(htfsN - htfs)), freqs = join(round(htfsN, 4), sep = ' '));
    if (max(abs(htfsN - htfs)) <= eps) break;
    htfs = htfsN;
  }
  return(list(htfs = htfsN, converged = i < Nitmax, eps = eps, iterations = i));	
}

#
#	<p> multi-locus reconstruction
#

wal = function(v)which(as.logical(v))
dfNameIdx = function(v, i, prefix = 'a')SetNames(Df(v), names = paste0(prefix, i))

# recursively form combinations based on herozygosity
# loci by column
# noParentOfOrigin: ignore parent-of-origin status
alleleCombinations = function(allelePairs, noParentOfOrigin = T) {
  als = allelePairs[, 1];
  doFlip = als[1] != als[2] && !noParentOfOrigin;
  combsRest = if (ncol(allelePairs) > 1)
    alleleCombinations(allelePairs[, -1, drop = F], noParentOfOrigin = F) else c();
  # rely on recyling
  all = if (doFlip) rbind(cbind(als, combsRest), cbind(rev(als), combsRest)) else cbind(als, combsRest);
  return(all);
}

# combine single pair
combineGtsP = function(g, N) {
  # split genotype into allele pairs
  pairs = mapply(pairFromGt, N, g);
  # allele combinations -> new haplotypes
  combs = alleleCombinations(pairs);
  # map allele tuple to haplotype enumeration
  Nmap = pop(c(1, N));
  hts = matrix(apply(Nmap * t(combs), 2, sum), ncol = 2, byrow = T);
  # psirs of haplotypes to diplotype
  gts = apply(hts, 1, gtFromPair, Nalleles = prod(N));
  # gt numbers to indicator vector
  gtIndicators = vector.assign(0, gts + 1, 1, N = Npairs(prod(N)));
  #print(list(g = g, hts = hts, gts = gts, gtIndicators = gtIndicators, N = length(gtIndicators)));
  #print(gtIndicators);
  return(gtIndicators);
}

combineGtsV = function(g1, g2) {
  N = c(Npairs2N(length(g1)), Npairs2N(length(g2)));
  # genotype combinations
  combs = merge(data.frame(g1 = wal(g1) - 1), data.frame(g2 = wal(g2) - 1));
  # diplotype combinations
  dts = apply(apply(combs, 1, combineGtsP, N = N), 1, sum);
  #dts = apply(combs, 1, combineGtsP, N = N);
  #print(dts);
  return(dts);
}

combineGts = function(m1, m2) {
  rs = lapply(Seq(1, nrow(m1)), function(i)combineGtsV(m1[i, ], m2[i, ]));
  return(do.call(rbind, rs));
}


merge.multi.list = function(l, .col.names = NULL, .col.names.prefix = "X",
                            .return.lists = F, .first.constant = T, stringsAsFactors = F, .cols.asAre = F, .constraint = NULL, ...) {
  # <p> determine column names of final data frame
  .col.names.generic = paste(.col.names.prefix, 1:length(l), sep = "");
  if (is.null(.col.names)) .col.names = names(l);
  if (is.null(.col.names)) .col.names = .col.names.generic;
  .col.names[.col.names == ""] = .col.names.generic[.col.names == ""];
  names(l) = .col.names;		# overwrite names
  # <p> construct combinations
  if (.first.constant) l = rev(l);
  df0 = data.frame();
  if (length(l) >= 1) for (i in 1:length(l)) {
    newNames = if (.cols.asAre) names(l[[i]]) else names(l)[i];
    # <p> prepare data.frame: handle lists as well as data.frames
    # <!> changed 22.3.2016
    #dfi = if (is.list(l[[i]])) unlist(l[[i]]) else l[[i]];
    dfi = if (!is.data.frame(l[[i]])) unlist(l[[i]]) else l[[i]];
    df1 = data.frame.types(dfi, names = newNames, stringsAsFactors = stringsAsFactors);
    # <p> perform merge
    df0 = if (i > 1) merge(df0, df1, ...) else df1;
  }
  if (.first.constant) df0 = df0[, rev(names(df0)), drop = F];
  if (.return.lists) df0 = apply(df0, 1, as.list);
  if (!is.null(.constraint)) {
    df0 = df0[apply(df0, 1, function(r).do.call(.constraint, as.list(r))), ];
  }
  df0
}



if (0) {
  print(dtContribMat(1));
  print(dtContribMat(2));
}

# enumeration of genotypes for allele pairs
if (0) {
  Nalleles = 3;
  gts = merge.multi(0:(Nalleles - 1), 0:(Nalleles - 1));
  gtEnum = apply(gts, 1, gtFromPair, Nalleles = Nalleles);
  gtsRt = do.call(rbind, lapply(gtEnum, function(gtE)pairFromGt(Nalleles, gtE)));
  print(cbind(gtEnum, gts, gtsRt));
  
}

if (T) {
  # two alleles: 0, 1 -> three genotypes: 0, 1, 2 -> columns 1, 2, 3 below
  # order of genotypes defined by gtFromPair
  gts = matrix(
    c(0, 1, 1, # 1st individual with genotypes 2 or 3
      1, 0, 0,
      0, 1, 1,
      0, 1, 0), ncol = 3, byrow = T);
  gts1 = matrix(
    c(0, 1, 1, # 1st individual with genotypes 2 or 3
      1, 0, 0,
      1, 1, 0,
      0, 0, 1), ncol = 3, byrow = T);
}
if (0) {
  afs = alleleFreqEstEM(gts);
  print(afs);
}

# Example2: should be .5
if (0) {
  afs = alleleFreqEstEM(gts1);
  print(afs);
}

if (0) {
  gts2 = combineGts(gts, gts1);
  print(gts2);
  afs = alleleFreqEstEM(gts2);
  print(afs);
}
