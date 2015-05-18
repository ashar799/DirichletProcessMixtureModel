library(hgu133a.db)
library(limma)
library(affy)
library(graph)
library(ipred)
library(gdata)
library(penalized)
library(org.Hs.eg.db)
library(penalizedSVM)
library(CoxBoost)
library(ggplot2)
library(BioNet)

options(error=recover)

# The working directory
WDIR="."

test.signature = function(signature, DIR=WDIR, DIROUT=WDIR, signame="", cutoff=0.05,
                          tests=c("DT", "DO"), map2entrez=NULL){        
  require(GOstats)
  require(gdata)        
  sig = setdiff(map2entrez[signature], NA)
  if(length(sig) > 0)
    write.csv(sig, file.path(DIR, paste("Signature_", signame, ".csv", sep="")))
  else
    write.csv(signature, file.path(DIR, paste("Signature_", signame, ".csv", sep="")),
              sep="\t")        
  print(head(signature))
  if(length(signature) == 0)
    tests = c()
  if("DT" %in% tests){                
    if(length(signature) == 0)
      p.value.dt = 1
    else{                        
      if(length(sig) > 0)
        signature = sig                        
      allgenes = unique(map2entrez)
      drugs = read.csv(file.path(DIR, "DrugTargets.csv"), skip=1, sep="\t")                        
      consensus.drugTargets = table(signature %in%
                                      as.character(setdiff(unique(drugs$Chosen.Gene.IDs.EntrezGene.IDs.), NA)))                        
      others.drugTargets = table(setdiff(allgenes, signature) %in%
                                   as.character(setdiff(unique(drugs$Chosen.Gene.IDs.EntrezGene.IDs.), NA)))
      conf.table = rbind(consensus.drugTargets[c(2,1)], others.drugTargets[c(2,1)])
      conf.table[is.na(conf.table)] = 0                                                
      F = fisher.test(conf.table, alternative="g")
      cat("Enrichment of drug targets for signature ", signame, ":\n")
      print(F)
      if(signame != "")
        save(F, file=file.path(DIROUT, paste(signame, "_DrugTargetEnrichment.rda")))
      p.value.dt = F$p.value
      names(p.value.dt) = "Enrichment of drug targets"
    }
  }
  else
    p.value.dt = c()
  if("DO" %in% tests){                                
    DOlite = read.xls(file.path(DIR, "do_lite.xls"))                                
    DOterms = c("Cancer", "Advanced cancer", "Breast cancer")                
    allgenes = unique(as.character(DOlite$Entrez))                
    if(length(signature) == 0)
      pvalues = rep(1, length(DOterms))
    else
      if(length(sig) > 0)
        signature = sig
    signature = intersect(signature, allgenes)
    pvalues = as.numeric(sapply(DOterms, function(DO){                                
      DO.match = table(signature %in% as.character(DOlite[DOlite$Disease == DO,
                                                          "Entrez"]))
      others = table(setdiff(allgenes, signature) %in%
                       as.character(DOlite[DOlite$Disease == DO, "Entrez"]))
      conf.table = rbind(DO.match[c(2,1)], others[c(2,1)])
      conf.table[is.na(conf.table)] = 0
      print(conf.table)
      F = fisher.test(conf.table, alternative="g")$p.value
      print(F)
      F
    }))
    names(F) = "Enrichment of disease suspectible genes"
    names(pvalues) = DOterms
    pvalues = p.adjust(pvalues, "bonf")
  }
  else
    pvalues = c()        
  list(p.value.dt=p.value.dt, p.value.DO=pvalues, p.value.KEGG=p.value.KEGG,
       p.value.GO=p.value.GO)
}


fit.survival = function(XX, Y, lambda1, lambda2=0, standardize=FALSE,
                        coef.penalty=NULL){        
  cv.loglik. = function(lambdas, ...){
    cat("lambdas = ", lambdas, "\n")
    if(length(lambdas) > 1){
      if(is.null(coef.penalty))
        Q = tryCatch(cvl(Y, XX, fold=5, lambda1=2^lambdas[1], lambda2=2^lambdas[2],
                         standardize=standardize, maxiter=200), error=function(e) list(cvl=-Inf,
                                                                                       model=NULL))
      else
        Q = tryCatch(cvl(Y, XX, fold=5, lambda1=2^lambdas[1]*coef.penalty,
                         lambda2=2^lambdas[2]*coef.penalty, standardize=standardize, maxiter=200),
                     error=function(e) list(cvl=-Inf, model=NULL))
    }
    else{
      if(is.null(coef.penalty))
        Q = tryCatch(cvl(Y, XX, fold=5, lambda1=2^lambdas, standardize=standardize,
                         maxiter=200), error=function(e) list(cvl=-Inf, model=NULL))
      else
        Q = tryCatch(cvl(Y, XX, fold=5, lambda1=2^lambdas*coef.penalty,
                         standardize=standardize, maxiter=200), error=function(e) list(cvl=-Inf,
                                                                                       model=NULL))
    }
    if(Q$cvl == -Inf)
      Q$cvl = -1e16
    cat("CV likelihood = ", Q$cvl, "\n")
    
    ret <- list(q.val = -Q$cvl, model = Q$fullfit)
    class(ret) <- "penalized"
    return(ret)
  }        
  # determine optimal hyperparameters via EPSGO algorithm
  options(error=dump.frames)
  if(length(lambda2) > 1){
    bounds = t(data.frame(log2lambda1=lambda1, log2lambda2=lambda2))
    colnames(bounds)<-c("lower", "upper")
    print(bounds)
    ep = try(EPSGO(cv.loglik., bounds=bounds, parms.coding=c("log2", "log2"),
                   seed=12345), silent=TRUE)
  }
  else{
    bounds = t(data.frame(log2lambda1=lambda1))
    colnames(bounds)<-c("lower", "upper")
    print(bounds)
    ep = try(EPSGO(cv.loglik., bounds=bounds, parms.coding=c("log2"), seed=12345),
             silent=TRUE)                
  }                        
  options(error=recover)        
  # train penalized object
  if(!is.null(ep)){                        
    if(length(lambda2) > 1){
      if(is.null(coef.penalty))
        fit = penalized(Y, XX, lambda1=2^ep$xmin[1], lambda2=2^ep$xmin[2],
                        standardize=standardize)
      else
        fit = penalized(Y, XX, lambda1=2^ep$xmin[1]*coef.penalty,
                        lambda2=2^ep$xmin[2]*coef.penalty, standardize=standardize)
    }
    else{
      if(is.null(coef.penalty))
        fit = penalized(Y, XX, lambda1=2^ep$xmin[1], standardize=standardize)
      else
        fit = penalized(Y, XX, lambda1=2^ep$xmin[1]*coef.penalty, standardize=standardize)
    }
  }                
  else{
    stop("No model could be fitted!!!")                                
  }
  fit
}


fastGeneRank = function(ex, Aq, d=0.5){        
  b = (1 - d)*ex                        
  r = solve(Aq, b)
  as.numeric(r)
}

compute.geneRanks = function(ex, net.genes, Aq){                
  net.genes = intersect(net.genes, names(ex))
  ex = ex[net.genes]                        
  gr = fastGeneRank(ex, Aq)
  names(gr) = net.genes
  gr
}

perform.cox = function(x, surv){
  p = c()
  for(i in 1:ncol(x)){                        
    cat(".")
    cox = tryCatch(coxph(surv ~ x[,i]), error=function(e) NULL)
    if(is.null(cox))
      p = c(p, 1)
    else
      p = c(p, summary(cox)$logtest["pvalue"])
  }
  p[is.na(p)] = 1                
  names(p) = colnames(x)
  p
}

# DIR: directory where to find datasets
init = function(KEGG=FALSE, DIR=WDIR){
  cat("initializing networks and precomputing QR decomposition for GeneRank
      calculation ...")
  if(!KEGG){
    load(file=file.path(DIR, "PathwayCommons_graph.rda"))
    KEGGstring = ""
  }
  else{
    load(file=file.path(DIR, "KEGGgraphTotalUndirected.rda"))
    KEGGstring = "_KEGG"
    nodes(G) = sub("hsa:", "", nodes(G))
  }
  require(RBGL)
  require(Matrix)
  G = ugraph(G)
  Gmat = as(G, "matrix")        
  mapping = data.frame(probesetID=names(map2entrez), graphID=map2entrez)
  mapping = mapping[!is.na(mapping$graphID), ]
  common.nodes = intersect(mapping$graphID, nodes(G))
  W = Gmat[common.nodes, common.nodes]
  cs = colSums(W)
  W = W[cs > 0, cs > 0]
  # precompute several things for speed up reasons
  D = diag(1/colSums(W))
  A =  Matrix(diag(ncol(W))) - 0.5*Matrix(W)%*%Matrix(D)
  Aq = qr(A)                        
  save(W, Aq, file=file.path(DIR, paste("GeneRanks", KEGGstring, ".rda", sep="")))
  cat("done\n")
}


# x = expression / gene rank matrix (samples x genes)
# times = survival times
# events = 1:event; 0:no event
crossval = function(x, times, events, nfolds=10, repeats=10, method="CoxBoost",
                    W=NULL, Aq=NULL, map2entrez=NULL, DIR=WDIR){
  if(method %in% c("BioNet", "BioNet (KEGG)"))
    G = as(W, "graphNEL")
  n = nrow(x)        
  ibrier = matrix(0, ncol=repeats, nrow=nfolds)
  genes = list()        
  orig.genes = list()
  fits = list()
  kk = 1        
  coef.penalty = NULL
  set.seed(1234567)                
  if(method %in% c("penGenFilter", "penGenFilter (KEGG)", "penGenFilterDG",
                   "penGenRank", "penGenRank (KEGG)", "CoxBoost", "CoxBoost (KEGG)", "filter (net.
                   genes)", "filter (net. genes KEGG)", "CoxBoost (net. genes)", "CoxBoost (net. genes
                   KEGG)", "penalized (net. genes)", "penalized (net. genes KEGG)", "BioNet", "BioNet
                   (KEGG)")){
    x = x[, !is.na(map2entrez[colnames(x)])]
    x = t(avereps(t(x), ID=map2entrez[colnames(x)]))
    if(method %in% c("penGenFilter", "penGenFilter (KEGG)", "penGenRank", "penGenRank
                     (KEGG)", "filter (net. genes)", "filter (net. genes KEGG)", "CoxBoost (net.
                     genes)", "CoxBoost (net. genes KEGG)","penalized (net. genes)", "penalized (net.
                     genes KEGG)", "BioNet", "BioNet (KEGG)")){
      net.genes = intersect(rownames(W), colnames(x))
      x = x[,net.genes]
}
}                
for(r in 1:repeats){                
  perm = sample(1:n)                                
  for(k in 1:nfolds){
    cat("CV fold ", k," (repeat ", r, ")\n\n")
    if(nfolds > 1){
      tst = perm[seq(k, n, by=nfolds)]                        
      trn = setdiff(1:n, tst)                                        
    }
    else{
      tst = 1:n
      trn = 1:n
    }
    xtrn = x[trn,,drop=F]
    xtst = x[tst,,drop=F]                
    if(method %in% c("CoxBoost", "CoxBoost (KEGG)")){
      cidx = match(colnames(W), colnames(x))
      optim.res = optimCoxBoostPenalty(times[trn], events[trn], x[trn,], K=5,
                                       maxstepno=1000, pendistmat=W, connected.index=cidx)
      cat("optimal penalty =", optim.res$penalty, "\n")
      fit <- CoxBoost(times[trn], events[trn], x[trn,], stepno=max(1,
                                                                   optim.res$cv.res$optimal.step), penalty=optim.res$penalty, pendistmat=W,
                      connected.index=cidx)
      if(NROW(xtst) > 0)
        pred = predict(fit, xtst, type="risk", times=times[tst])
      genes.idx = which(fit$coefficients[nrow(fit$coefficients),] != 0)
      genes[[kk]] = colnames(x)[genes.idx]         
      orig.genes[[kk]] = colnames(xtrn)
    }                                                
    else if(method %in% c("CoxBoost (w.o. netw.)", "CoxBoost (net. genes)", "CoxBoost
                          (net. genes KEGG)")){
      optim.res = optimCoxBoostPenalty(times[trn], events[trn], x[trn,], K=5,
                                       maxstepno=1000)
      fit <- CoxBoost(times[trn], events[trn], x[trn,], stepno=max(1,
                                                                   optim.res$cv.res$optimal.step), penalty=optim.res$penalty)
      if(NROW(xtst) > 0)
        pred = predict(fit, xtst, type="risk", times=times[tst])
      genes.idx = which(fit$coefficients[nrow(fit$coefficients),] != 0)
      genes[[kk]] = colnames(x)[genes.idx]                        
      orig.genes[[kk]] = colnames(xtrn)
  }
  else if(method %in% c("penGenFilter", "penGenFilter (KEGG)", "filter", "filter
                        (net. genes)", "filter (net. genes KEGG)")){ # GeneRank als Filter                                
    p = perform.cox(xtrn, Surv(times[trn], event=events[trn]))
    qs = c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025)
    if(method %in% c("penGenFilter", "penGenFilter (KEGG)")){
      raw.ranks = -log(p + 1e-100)
      gr = compute.geneRanks(raw.ranks, rownames(W), Aq=Aq)                                        
    }
    else{
      qs = round(qs*ncol(xtrn))
      qs = qs[qs < nrow(xtrn) & qs > 0]
    }                                                                                
    best = Inf
    for(q in qs){
      if(method %in% c("filter", "filter (net. genes)", "filter (net. genes KEGG)")){
        sel.tmp = names(sort(p)[1:q])
      }
      else{
        sel.tmp = names(sort(gr, decreasing=T)[1:round(min(c(q*length(gr),
                                                             NROW(xtrn)-1)))])                                                
      }                                        
      xtrn.tmp = scale(xtrn[, sel.tmp, drop=F])                                                                                
      cox.total.tmp = try(coxph(Surv(times[trn], event=events[trn]) ~ xtrn.tmp),
                          silent=T)
      if(!is(cox.total.tmp, "try-error")){
        bic = -cox.total.tmp$loglik[2] + length(sel.tmp)*NROW(xtrn)/(NROW(xtrn) -
                                                                       length(sel.tmp) - 1) # AICc criterion                                                
        print(bic)
        if(bic < best){
          best = bic
          sel = sel.tmp                                                                                                
          mytrn = xtrn.tmp
        }                                                        
      }
    }                                                                                                        
    if(NROW(xtst) > 0){                                                
      xtst = scale(xtst[, sel, drop=F], center=attr(mytrn, "scaled:center"),
                   scale=attr(mytrn, "scaled:scale"))
    }                                                                                                                                                
    cox.total = coxph(Surv(times[trn], event=events[trn]) ~ mytrn,
                      control=coxph.control(iter.max=100))
    fit = survfit(cox.total, as.data.frame(xtst))
    class(fit) = "survfit"
    genes[[kk]] = colnames(mytrn)
    orig.genes[[kk]] = colnames(xtrn)
    print(head(genes[[kk]]))
  }
  else if(method %in% c("penalized", "penalized (net. genes)", "penalized (net.
                        genes KEGG)", "penGenRank", "penGenRank (KEGG)", "penPath", "penAvg", "BioNet",
                        "BioNet (KEGG)", "sigPath", "penCons")){                                                                                                                                                                                                                                                                                                        
    if(method %in% c("penGenRank", "penGenRank (KEGG)")){                        # GeneRank als embedded
      Methode                                                                
      p = perform.cox(xtrn, Surv(times[trn], event=events[trn]))
      gr = compute.geneRanks(-log(p + 1e-100), rownames(W),
                             Aq=Aq)                                                                                                                                                                                                
      coef.penalty = abs(1/gr)                                
      xtrn = xtrn[,names(gr)]
      xtst = xtst[,names(gr)]
    }                                
    else if(method %in% c("BioNet", "BioNet (KEGG)")){                                        
      p = perform.cox(xtrn, Surv(times[trn], event=events[trn]))                                                                        
      bum = fitBumModel(p, plot=TRUE)
      sc = scoreFunction(bum, fdr=0.1)                                                                                        
      if(!all(sc < 0) && !all(sc == Inf)){                                                                                                                        
        G = as(W, "graphNEL")
        Gsub = runFastHeinz(G, sc)
        print(Gsub)
        sel = nodes(Gsub)
        xtrn = xtrn[, sel, drop=F]
        if(NROW(xtst) > 0)
          xtst = xtst[, sel, drop=F]
      }
      else{ # leeres Modell
        xtrn = as.matrix(rep(1, length(trn)))
        xtst = as.matrix(rep(1, length(tst)))
      }                                                                                                                        
    }
    else if(method == "penPath"){        # perform global test to pre-select GENES                                                                                        
      mytest = result(gtKEGG(Surv(times[trn], events[trn]), x[trn,],
                             id=Rkeys(hgu133aPATH), annotation="hgu133a.db", multtest="BH"))
      print(mytest)
      mypaths = setdiff(rownames(mytest)[mytest$BH < 0.1], NA)
      if(length(mypaths) > 0){                                        
        sel = unlist(mget(mypaths, revmap(hgu133aPATH),ifnotfound=NA))
        xtrn = xtrn[,sel, drop=F]
        if(NROW(xtst) > 0)
          xtst = xtst[,sel, drop=F]
      }
      else{
        xtrn = as.matrix(rep(1, length(trn)))
        xtst = as.matrix(rep(1, length(tst)))
      }
    }                        
    else if(method %in% c("sigPath")){ # perform global test to pre-select PATHWAYS                                
      mytest = result(gtKEGG(Surv(times[trn], events[trn]), x[trn,],
                             id=Rkeys(hgu133aPATH), annotation="hgu133a.db", multtest="BH"))
      print(mytest)
      mypaths = setdiff(rownames(mytest)[mytest$BH < 0.1], NA)
      if(length(mypaths) > 0){                                                                                        
        sel = mget(mypaths, revmap(hgu133aPATH),ifnotfound=NA)
        xtrn = sapply(sel, function(s) rowMeans(xtrn[,s,drop=F]))                                        
        if(NROW(xtst) > 0)
          xtst = sapply(sel, function(s) rowMeans(xtst[,s,drop=F]))
      }
      else{
        xtrn = as.matrix(rep(1, length(trn)))
        xtst = as.matrix(rep(1, length(tst)))
      }
    }
    else if(method == "penAvg"){ # average PATHWAY expression
      sel = mget(Rkeys(hgu133aPATH), revmap(hgu133aPATH),ifnotfound=NA)
      xtrn = sapply(sel, function(s) rowMeans(xtrn[,s,drop=F]))
      if(NROW(xtst) > 0)
        xtst = sapply(sel, function(s) rowMeans(xtst[,s,drop=F]))
    }
    else if(method == "penCons"){ # penalized consensus signature
      consensus = read.csv(file.path(DIR, "Consensus_annotated.csv"))$entrez
      sel = setdiff(unlist(mget(as.character(consensus), revmap(hgu133aENTREZID),
                                ifnotfound=NA)), NA)
      xtrn = xtrn[,sel,drop=F]
      if(NROW(xtst) > 0)
        xtst = xtst[,sel, drop=F]        
    }
    if(!all(xtrn == 1)) # empty model
      xtrn = scale(xtrn)
    if(NROW(xtst) > 0 & !all(xtst == 1))
      xtst = scale(xtst, center=attr(xtrn, "scaled:center"), scale=attr(xtrn,
                                                                        "scaled:scale"))                                                                
    if(NCOL(xtrn) > 1)
      fit = fit.survival(xtrn, Surv(times[trn], events[trn]), lambda1=c(-10,10),
                         lambda2=c(-10,0), standardize=FALSE, coef.penalty=coef.penalty)
    else
      fit = penalized(Surv(times[trn], events[trn]), xtrn, lambda1=1e-100, lambda2=0,
                      standardize=FALSE)
    if(NROW(xtst) > 0){
      pred1 = predict(fit, xtst)
      pred = sapply(times[tst], function(t) survival(pred1, t))
    }
    if(method %in% c("penAvg", "sigPath") && nfolds == 1){                                        
      genes[[kk]] = unlist(mget(names(coef(fit)), revmap(hgu133aPATH),
                                ifnotfound=NA))                                        
    }
    else{                                                                                                                                        
      genes[[kk]] = names(coef(fit))                                        
    }
    orig.genes[[kk]] = colnames(xtrn)
  }                                                                
  if(NROW(xtst) > 0){
    if(method %in% c("penGenFilter", "penGenFilter (KEGG)", "filter", "filter (net.
                     genes)", "filter (net. genes KEGG)")){                                        
      ibrier[k,r] = sbrier(Surv(times[tst], event=events[tst]), fit)
  }
  else
    ibrier[k,r] = sbrier(Surv(times[tst], events[tst]), t(pred))                        
  cat("Integrated Brier score (test set) = ", ibrier[k,r], "\n\n")
  }
  fits[[kk]] = fit
  kk = kk + 1
  }                                
}                        
list(ibrier=ibrier, genes=genes, fit=fits, orig.genes=orig.genes)
}

prediction.analysis = function(dataset, ex, times, events, map2entrez, DIR=WDIR){        
  print(dataset)
  methods = c("penalized", "penalized (net. genes)", "penalized (net. genes KEGG)",
              "filter", "filter (net. genes)", "filter (net. genes KEGG)", "sigPath", "penPath",
              "penAvg", "penGenRank", "penGenRank (KEGG)", "penGenFilter", "penGenFilter (KEGG)",
              "BioNet", "BioNet (KEGG)", "CoxBoost", "CoxBoost (KEGG)", "CoxBoost (w.o. netw.)",
              "CoxBoost (net. genes)", "CoxBoost (net. genes KEGG)", "penCons")                
  results = data.frame()                
  for(method in methods){
    print(method)
    if(file.exists(file.path(DIR, paste(dataset, "_", method, "_cvresults.rda",
                                        sep="")))){                        
      load(file.path(DIR, paste(dataset, "_", method, "_cvresults.rda", sep="")))                        
    }                        
    else{                        
      if(length(grep("KEGG", method)) > 0)
        load(file=file.path(DIR, "GeneRanks_KEGG.rda"), envir=globalenv())
      else
        load(file=file.path(DIR, "GeneRanks.rda"), envir=globalenv())                                        
      cv.res = crossval(ex, times, events, method=method, W=W, Aq=Aq,
                        map2entrez=map2entrez)                        
      save(cv.res, file=file.path(DIR, paste(dataset, "_", method, "_cvresults.rda",
                                             sep="")))
    }                
    results = rbind(results, data.frame(integrated.brier=colMeans(cv.res$ibrier),
                                        method=method, dataset=dataset))
  }                
  qplot(method, integrated.brier, data=results, fill=method, position="dodge",
        geom="boxplot", xlab="", ylim=c(0,1), ylab="integraded Brier score",
        main=paste(dataset, "- Prediction Performance (10 x 10-fold CV)")) + theme_bw() +
    coord_flip() + opts(legend.position="none", strip.text.y = theme_text(angle=0)) 
  ggsave(file=file.path(DIR, paste(dataset, "_cvresults.pdf")))
  save(results, file=file.path(DIR, paste(dataset, "_cvresultsAll.rda")))
}

prediction.analysis.crossdataset = function(mydataset, ex, times, events,
                                            datasets=setdiff(c("Desmedt et al.", "Ivshina et al.", "Pawitan et al.", "Schmidt et
                                                               al.", "Wang et al."), mydataset), DIR=WDIR){        
  methods = c("penalized", "penAvg", "penGenFilter", "CoxBoost", "penCons")
  methods2 =  c("ela. net w.o. netw.", "AvgPath", "GeneRankFilter", "CoxBoost",
                "ConsSig")
  ibrier = matrix(0, ncol=length(datasets), nrow=length(methods))
  dimnames(ibrier) = list(methods, datasets)
  brier.scores = data.frame()
  for(method in methods){                                
    if(file.exists(file.path(DIR, paste(mydataset, "_", method, "_totalfit.rda",
                                        sep="")))){
      load(file.path(DIR, paste(mydataset, "_", method, "_totalfit.rda", sep="")))                        
    }                
    if(method == "penAvg"){                                
      pathways2genes = mget(unlist(fit$orig.genes), revmap(hgu133aPATH), ifnotfound=NA)
      xtrn = sapply(pathways2genes, function(s) rowMeans(ex[,s,drop=F]))
    }
    else
      xtrn = ex
    xtrn = scale(xtrn)
    for(dataset in datasets){
      cat("method = ", method, "dataset = ", dataset, "\n")
      if(dataset == "Schmidt et al."){
        load(file.path(DIR, "Mainz_FARMS.rda"))
        xtst = t(exprs(an))
        ttimes = pData(an)$t.dmfs
        tevents = pData(an)$e.dmfs                                
      }
      if(dataset == "Pawitan et al."){
        load(file.path(DIR, "PawitanGSE1456_FARMS.rda"))
        xtst = t(exprs(an))
        ttimes = pData(an)$SURV_RELAPSE..
        tevents = pData(an)$RELAPSE..
      }
      if(dataset == "Ivshina et al."){
        load(file.path(DIR, "GSE4922_U_FARMS.rda"))
        xtst = t(exprs(an))
        ttimes = pData(an)$DFS.TIME..yrs. 
        tevents =
          pData(an)$DFS.EVENT..0.censored..1.event.defined.as.any.type.of.recurrence..local..regional.or.distant..or.death.from.breast.cancer
      }
      if(dataset == "Desmedt et al."){
        load(file.path(DIR, "GSE7390_Desmedt_FARMS.rda"))
        xtst = t(exprs(an))
        ttimes = pData(an)$t.rfs
        tevents = pData(an)$e.rfs
      }
      if(dataset == "Wang et al."){
        load(file.path(DIR, "Wang_FARMS.rda"))
        xtst = t(exprs(an))
        ttimes = pData(an)$time.to.relapse.or.last.follow.up..months.
        tevents = pData(an)$relapse..1.True.
      }                        
      if(method == "penAvg"){
        xtst = sapply(pathways2genes, function(s) rowMeans(xtst[,s,drop=F]))
      }                        
      xtst = scale(xtst, center=attr(xtrn, "scaled:center"), scale=attr(xtrn,
                                                                        "scaled:scale"))
      if(method == "penGenFilter"){
        ibrier[method,dataset] = sbrier(Surv(ttimes, tevents), fit$fit)
        bs = sapply(1:length(ttimes), function(t) sbrier(Surv(ttimes, tevents), fit$fit,
                                                         btime=ttimes[t]))
      }
      else{
        if(method %in% c("penCons")){
          pred1 = predict(fit$fit, xtst[, unlist(fit$orig.genes), drop=F])
          pred = sapply(ttimes, function(t) survival(pred1, t))
        }
        else if(method %in% c("penalized", "penAvg")){                                        
          pred1 = tryCatch(predict(fit$fit, xtst), error=function(e) predict(fit$fit,
                                                                             xtst[,names(coef(fit$fit)),drop=F])) 
          pred = sapply(ttimes, function(t) survival(pred1, t))
        }
        else if(method == "CoxBoost")
          pred = predict(fit$fit, xtst, type="risk", times=ttimes)
        ibrier[method,dataset] = sbrier(Surv(ttimes, tevents), t(pred))
        bs = sapply(1:length(ttimes), function(t) sbrier(Surv(ttimes, tevents),
                                                         t(pred)[t,,drop=F], btime=ttimes[t]))
      }
      bs = bs + rep(runif(1, min=-0.01, max=0.01), length(bs)) # small jitter, just to
      see curves separated
      cat("Integrated Brier score  = ", ibrier[method,dataset], "\n\n")                                        
      df = unique(data.frame(method=method, dataset=dataset, time=ttimes, BS=bs))
      brier.scores = rbind(brier.scores, df=df[order(df$time),])
    }                
  }
  rownames(ibrier) = methods2
  levels(brier.scores$method) = methods2
  brier.scores$time[brier.scores$dataset == "Schmidt et al."] =
    brier.scores$time[brier.scores$dataset == "Schmidt et al."]/365
  brier.scores$time[brier.scores$dataset == "Desmedt et al."] =
    brier.scores$time[brier.scores$dataset == "Desmedt et al."]/365
  brier.scores$time[brier.scores$dataset == "Wang et al."] =
    brier.scores$time[brier.scores$dataset == "Wang et al."]/12
  qplot(time, BS, data=brier.scores, col=method, facets=~dataset, geom="line",
        xlab="time (years)", ylab="Brier score", main=paste("Prediction performance -
                                                            models trained on ", mydataset)) + theme_bw() + theme(strip.text.y =
                                                                                                                    element_text(angle=0)) + scale_color_grey() + geom_point(aes(shape=method))        
  ggsave(file=file.path(DIR, paste("PerfCrossdataset_TrainedOn", mydataset, ".pdf",
                                   sep="")))
  save(ibrier, brier.scores, file=file.path(DIR, paste("PerfCrossdataset_TrainedOn",
                                                       mydataset, ".rda", sep="")))        
  }

perf.comparison = function(resultsDO, resultsDT, datasets=c("Desmedt et al.",
                                                            "Ivshina et al.", "Pawitan et al.", "Schmidt et al.", "Wang et al."), DIR=WDIR){
  methods = c("penalized", "penalized (net. genes)", "penalized (net. genes KEGG)",
              "filter", "filter (net. genes)", "filter (net. genes KEGG)", "CoxBoost (w.o.
              netw.)", "sigPath", "penPath", "penAvg", "penGenRank", "penGenRank (KEGG)",
              "penGenFilter", "penGenFilter (KEGG)", "BioNet", "BioNet (KEGG)", "CoxBoost",
              "CoxBoost (KEGG)", "CoxBoost (net. genes)", "CoxBoost (net. genes KEGG)",
              "penCons")
  methods2 = c("ela. net w.o. netw.", "ela. net. net. genes", "ela. net. KEGG genes",
               "CoxFilter w.o. netw.", "CoxFilter net. genes", "CoxFilter KEGG genes", "CoxBoost
               w.o. netw.", "sigPathAvg", "sigPathGene", "AvgPath", "GeneRankWeight",
               "GeneRankWeight (KEGG)", "GeneRankFilter", "GeneRankFilter (KEGG)", "BioNet",
               "BioNet (KEGG)", "CoxBoost", "CoxBoost (KEGG)", "CoxBoost (net. genes)", "CoxBoost
               (KEGG genes)", "ConsSig")        
  results = data.frame()
  results.stab = data.frame()
  results.stab2 = data.frame()
  results.ngenes = data.frame()
  for(method in methods){
    for(dataset in datasets){
      cat("method = ", method, "dataset = ", dataset, "\n")
      load(file.path(DIR, paste(dataset, "_", method, "_cvresults.rda", sep="")))
      mymethod = methods2[which(methods == method)]
      results = rbind(results, data.frame(integrated.brier=colMeans(cv.res$ibrier),
                                          method=mymethod, dataset=dataset))
      #                        
      ngenes = sapply(cv.res$genes, length)
      cv.res$genes=sapply(cv.res$genes, unique)
      freq = sort(table(unlist(cv.res$genes)), decreasing=T)[1:10]                        
      H = hist(freq, breaks=seq(0,100, by=10), plot=F)$counts / length(freq)                        #
      convert counts into fractions!                                                
        SI = sum(freq, na.rm=T) / (sum(freq != 0, na.rm=T) * 100)                        
      results.stab = rbind(results.stab, data.frame(method=mymethod,
                                                    orig.method=method, SI=SI, dataset=dataset))                
      results.stab2 = rbind(results.stab2, data.frame(method=mymethod,
                                                      orig.method=method, appearance.freq=seq(10,100,by=10), fraction=H,
                                                      dataset=dataset))
      results.ngenes = rbind(results.ngenes, data.frame(method=mymethod,
                                                        orig.method=method, ngenes=ngenes, dataset=dataset))
    }
  }        
  qplot(method, SI, data=results.stab, fill=method, facets=~dataset, geom="bar",
        position="dodge", ylim=c(0, 1), xlab="", ylab="stability index (SI)",
        main="Signature stability (top 10 features)") + theme_bw() + coord_flip() +
    opts(legend.position="none", strip.text.y = theme_text(angle=0))
  ggsave(file=file.path(DIR, "SI.pdf", height=10))        
  #
  qplot(appearance.freq, fraction, data=results.stab2, colour=method,
        facets=~dataset, geom="point", xlab="times selected", ylim=c(0,1), ylab="relative
        frequency", main="Probeset selection stability (top 10)") + theme_bw() +
    geom_line(aes(colour=method))
  ggsave(file=file.path(DIR, "GeneSelectionStabilityTotal.pdf"), height=10)
  #        
  results.ngenes2 = results.ngenes[results.ngenes$ngenes != 0, ]
  qplot(method, log10(ngenes), data=results.ngenes2, fill=method, facets=~dataset,
        geom="boxplot", position="dodge", xlab="", ylab="log10 (no. genes)",
        main="Signature size") + theme_bw() + coord_flip() + opts(legend.position="none",
                                                                  strip.text.y = theme_text(angle=0))
  ggsave(file=file.path(DIR, paste("SignatureSizesTotal.pdf", sep="")), height=10)
  #
  empty.models = aggregate(ngenes~dataset + method, data=results.ngenes,
                           FUN=function(x) mean(x == 0))
  qplot(method, ngenes, data=empty.models, fill=method, facets=~dataset, geom="bar",
        position="dodge", xlab="", ylab="fraction of empty models", main="Fraction of empty
        models") + theme_bw() + coord_flip() + opts(legend.position="none", strip.text.y =
                                                      theme_text(angle=0))
  ggsave(file=file.path(DIR, paste("EmptyModelsTotal.pdf", sep="")), height=10)
  #
  qplot(method, integrated.brier, data=results, fill=method, facets=~dataset,
        position="dodge", geom="boxplot", xlab="", ylim=c(0,0.5), ylab="integrated Brier
        score", main= "Prediction Performance (10 x 10-fold CV)") + theme_bw() +
    coord_flip() + opts(legend.position="none", strip.text.y = theme_text(angle=0)) 
  ggsave(file=file.path(DIR, "CVTotal.pdf"), height=10)        
  #                
  resultsDO$method = as.character(resultsDO$method)
  for(m in 1:length(methods))
    resultsDO$method[resultsDO$method == methods[m]] = methods2[m]        
  resultsDO$method = factor(resultsDO$method, levels=unique(resultsDO$method))
  resultsDO$p.value.DO[resultsDO$p.value.DO == Inf] = 150
  qplot(method, p.value.DO, data=resultsDO, geom="bar", facets=~dataset, fill=Term,
        position="dodge", xlab="", ylab="-log10(p-value)", main="Enrichment of suspected
        disesase genes") + geom_hline(yintercept=-log10(0.05)) + coord_flip() + theme_bw()
  + opts(legend.position="top", strip.text.y = theme_text(angle=0))
  ggsave(file=file.path(DIR, "DiseaseGeneEnrichment.pdf"), height=10)
  #
  resultsDT$method = as.character(resultsDT$method)
  for(m in 1:length(methods))
    resultsDT$method[resultsDT$method == methods[m]] = methods2[m]        
  resultsDT$method = factor(resultsDT$method, levels=unique(resultsDT$method))
  resultsDT$p.value.DT[resultsDT$p.value.DT == Inf] = 150
  qplot(method, p.value.DT, data=resultsDT, geom="bar", facets=~dataset,
        position="dodge", xlab="", ylab="-log10(p-value)", main="Enrichment of known drug
        targets") + geom_hline(yintercept=-log10(0.05)) + coord_flip() + theme_bw() +
    opts(legend.position="none", strip.text.y = theme_text(angle=0))
  ggsave(file=file.path(DIR, "DTEnrichment.pdf"), height=10)                                
  
  sig = lapply(levels(results$dataset), function(d){
    pwt = matrix(NA, ncol=nlevels(results$method), nrow=nlevels(results$method))
    dimnames(pwt) = list(levels(results$method), levels(results$method))
    for(m1 in levels(results$method)){
      for(m2 in levels(results$method)){
        if(m1 != m2){
          pwt[m1, m2] = wilcox.test(results$integrated.brier[results$dataset ==d &
                                                               results$method == m1], results$integrated.brier[results$dataset ==d &
                                                                                                                 results$method == m2], paired=T, alternative="l")$p.value
        }
      }
    }                                
    pwt = matrix(p.adjust(pwt, method="BH"), ncol=ncol(pwt))
    dimnames(pwt) = list(levels(results$method), levels(results$method))                                                                                                                        
    write.csv(pwt, file=file.path(DIR, paste("PerformanceComparison", d, "csv",
                                             sep="")))
    pwt
  })
  
  require(MADAM)
  pwt2 = matrix(NA, ncol=ncol(sig[[1]]), nrow=nrow(sig[[1]]))
  dimnames(pwt2) = dimnames(sig[[1]])
  for(i in 1:ncol(sig[[1]])){                
    for(j in 1:nrow(sig[[1]])){
      if(i != j){
        pwt2[i, j] = fisher.method(t(as.matrix(sapply(sig, function(x) x[i,j]))),
                                   p.corr="none")$p.adj
      }                        
    }
  }                
  pwt2 = matrix(p.adjust(pwt2, method="BH"), ncol=ncol(pwt2))
  dimnames(pwt2) = dimnames(sig[[1]])
  write.csv(pwt2, file=file.path(DIR, "PerformanceComparisonAverage.csv"))
  
  
  pdf(file="ANOVAVerificationSI.pdf")
  plot(lm(SI~method + dataset, data=results.stab))
  dev.off()
  model2 = aov(SI~method + dataset, data=results.stab)        
  print(summary(model2)) 
  tuk = TukeyHSD(model2, which="method", ordered=TRUE) 
  tuk$method = tuk$method[order(tuk$method[,"diff"]),]
  tuk = data.frame(tuk$method, comparison=rownames(tuk$method))
  tuk$comparison = factor(tuk$comparison, levels=tuk$comparison)
  ggplot (tuk, aes (x=comparison, y=diff, ymin=lwr, ymax=upr)) + geom_pointrange () +
    ylab ("difference in mean SI") + theme_bw() + coord_flip ()
  ggsave(file=file.path(DIR, "TukeyMethodStability.pdf"), height=10)
  write.csv(tuk, file=file.path(DIR, "TukeyMethodComparisonStability.csv"))
}

interpretability.analysis = function(dataset, ex, times, events, map2entrez,
                                     DIR=WDIR){        
  resultsDO = data.frame()
  resultsDT = data.frame()        
  methods = c("penalized", "penalized (net. genes)", "penalized (net. genes KEGG)",
              "filter", "filter (net. genes)", "filter (net. genes KEGG)", "CoxBoost (w.o.
              netw.)", "sigPath", "penPath", "penAvg", "penGenRank", "penGenRank (KEGG)",
              "penGenFilter", "penGenFilter (KEGG)", "BioNet", "BioNet (KEGG)", "CoxBoost",
              "CoxBoost (KEGG)", "CoxBoost (net. genes)", "CoxBoost (net. genes KEGG)",
              "penCons")        
  for(method in methods){
    print(method)                
    if(file.exists(file.path(file.path(DIR, paste(dataset, "_", method,
                                                  "_totalfit.rda", sep=""))))){
      load(file.path(DIR, paste(dataset, "_", method, "_totalfit.rda", sep="")))                        
    }                        
    else{
      if(length(grep("KEGG", method)) > 0)
        load(file=file.path(DIR, "GeneRanks_KEGG.rda"))
      else
        load(file=file.path(DIR, "GeneRanks.rda"))                        
      drugs = read.csv(file.path(DIR, "DrugTargets.csv"), skip=1, sep="\t")                
      DOlite = read.xls(file.path(DIR, "do_lite.xls"))
      dgenes =
        union(as.character(setdiff(unique(drugs$Chosen.Gene.IDs.EntrezGene.IDs.), NA)),
              DOlite$Entrez[DOlite$Disease %in% c("Cancer", "Advanced Cancer", "Breast
                                                  Cancer")])
      fit = crossval(ex, times, events, method=method, nfolds=1, repeats=1, W=W, Aq=Aq,
                     map2entrez=map2entrez, L=L, dgenes=dgenes) # fit model to WHOLE dataset
      print(length(unlist(fit$genes)))                        
      mytest = test.signature(unlist(fit$genes), signame=paste(dataset, "_", method,
                                                               sep=""), map2entrez=map2entrez)  # enrichment analysis                        
      save(fit, mytest, file=file.path(DIR, paste(dataset, "_", method,
                                                  "_totalfit.rda", sep="")))
      print(fit)                        
    }                        
    print(length(unlist(fit$genes)))                        
    mytest = test.signature(unlist(fit$genes), signame=paste(dataset, "_", method,
                                                             sep=""), map2entrez=map2entrez)  # enrichment                        
    if(length(unlist(fit$genes)) > 0){
      
      if(length(mytest$p.value.DO) > 0)
        resultsDO = rbind(resultsDO, data.frame(p.value.DO=-log10(mytest$p.value.DO),
                                                method=method, Term=names(mytest$p.value.DO), dataset=dataset))                
      if(length(mytest$p.value.dt) > 0)
        resultsDT = rbind(resultsDT, data.frame(p.value.DT=-log10(mytest$p.value.dt),
                                                method=method, dataset=dataset))
    }
  }                        
  if(NROW(resultsDO) > 0){
    qplot(Term, p.value.DO, facets=method~., data=resultsDO, geom="bar",
          position="dodge", xlab="", ylab="-log10(p-value)", main=paste("DO Enrichment (",
                                                                        dataset, ")", sep="")) + geom_hline(yintercept=-log10(0.05)) + coord_flip() +
      theme_bw() + opts(legend.position="none", strip.text.y = theme_text(angle=0))
    ggsave(file=file.path(DIR, paste(dataset, "DOEnrichment.pdf", sep="")))                
  }
  
  if(NROW(resultsDT) > 0){
    qplot(method, p.value.DT, data=resultsDT, geom="bar", position="dodge", xlab="",
          ylab="-log10(p-value)", main=paste("Drug Target Enrichment (", dataset, ")",
                                             sep="")) + geom_hline(yintercept=-log10(0.05)) + coord_flip()+ theme_bw() +
      opts(legend.position="none")
    ggsave(file=file.path(DIR, paste(dataset, "DTEnrichment.pdf", sep="")))                
  }        
  list(resultsDO=resultsDO, resultsDT=resultsDT)        
}

#init(TRUE) # initialization for KEGG pathways => DONE (file GeneRanks_KEGG.rda)!
#init(FALSE) # initialization for PPI network => DONE (file GeneRanks.rda)!
#datasets = c("Desmedt et al.", "Ivshina et al.", "Pawitan et al.", "Schmidt et
al.", "Wang et al.") 
datasets = "Schmidt et al." # just one example dataset
resultsDO = data.frame()
resultsDT = data.frame()
for(dataset in datasets){
print(dataset)

if(dataset == "Schmidt et al."){
load(file.path(WDIR, "Mainz_FARMS.rda"))                
map2entrez = unlist(mget(featureNames(an), hgu133aENTREZID, ifnotfound=NA))
prediction.analysis("Schmidt et al.", t(exprs(an)), pData(an)$t.dmfs,
pData(an)$e.dmfs, map2entrez=map2entrez)                
res = interpretability.analysis("Schmidt et al.", t(exprs(an)), pData(an)$t.dmfs,
pData(an)$e.dmfs, map2entrez=map2entrez)                

}

if(dataset == "Pawitan et al."){
load(file.path(WDIR, "PawitanGSE1456_FARMS.rda"))
map2entrez = unlist(mget(featureNames(an), hgu133aENTREZID, ifnotfound=NA))
prediction.analysis("Pawitan et al.", t(exprs(an)), pData(an)$SURV_RELAPSE..,
pData(an)$RELAPSE.., map2entrez=map2entrez)
res = interpretability.analysis("Pawitan et al.", t(exprs(an)),
pData(an)$SURV_RELAPSE.., pData(an)$RELAPSE.., map2entrez=map2entrez)                
}

if(dataset == "Ivshina et al."){
load(file.path(WDIR, "GSE4922_U_FARMS.rda"))
map2entrez = unlist(mget(featureNames(an), hgu133aENTREZID, ifnotfound=NA))
prediction.analysis("Ivshina et al.", t(exprs(an)), pData(an)$DFS.TIME..yrs.,
pData(an)$DFS.EVENT..0.censored..1.event.defined.as.any.type.of.recurrence..local..regional.or.distant..or.death.from.breast.cancer,
map2entrez=map2entrez)
res = interpretability.analysis("Ivshina et al.", t(exprs(an)),
pData(an)$DFS.TIME..yrs.,
pData(an)$DFS.EVENT..0.censored..1.event.defined.as.any.type.of.recurrence..local..regional.or.distant..or.death.from.breast.cancer,
map2entrez=map2entrez)                
}

if(dataset == "Desmedt et al."){
load(file.path(WDIR, "GSE7390_Desmedt_FARMS.rda"))
map2entrez = unlist(mget(featureNames(an), hgu133aENTREZID, ifnotfound=NA))
prediction.analysis("Desmedt et al.",  t(exprs(an)), pData(an)$t.rfs,
pData(an)$e.rfs, map2entrez=map2entrez)
res = interpretability.analysis("Desmedt et al.", t(exprs(an)), pData(an)$t.rfs,
pData(an)$e.rfs, map2entrez=map2entrez)                
}

if(dataset == "Wang et al."){
load(file.path(WDIR, "Wang_FARMS.rda"))
map2entrez = unlist(mget(featureNames(an), hgu133aENTREZID, ifnotfound=NA))
prediction.analysis("Wang et al.", t(exprs(an)),
pData(an)$time.to.relapse.or.last.follow.up..months., pData(an)$relapse..1.True.,
map2entrez=map2entrez)
res = interpretability.analysis("Wang et al.", t(exprs(an)),
pData(an)$time.to.relapse.or.last.follow.up..months., pData(an)$relapse..1.True.,
map2entrez=map2entrez)
}

resultsDO = rbind(resultsDO, res$resultsDO)
resultsDT = rbind(resultsDT, res$resultsDT)
list(resultsDO=resultsDO, resultsDT=resultsDT)
}

perf.comparison(resultsDO, resultsDT, datasets=datasets)

for(dataset in datasets){
print(dataset)
if(dataset == "Schmidt et al."){
load(file.path(WDIR, "Mainz_FARMS.rda"))                                
prediction.analysis.crossdataset("Schmidt et al.", t(exprs(an)), pData(an)$t.dmfs,
pData(an)$e.dmfs)                
}

if(dataset == "Pawitan et al."){
load(file.path(WDIR, "PawitanGSE1456_FARMS.rda"))                
prediction.analysis.crossdataset("Pawitan et al.", t(exprs(an)),
pData(an)$SURV_RELAPSE.., pData(an)$RELAPSE..)                        
}

if(dataset == "Ivshina et al."){
load(file.path(WDIR, "GSE4922_U_FARMS.rda"))                
prediction.analysis.crossdataset("Ivshina et al.", t(exprs(an)),
pData(an)$DFS.TIME..yrs.,
pData(an)$DFS.EVENT..0.censored..1.event.defined.as.any.type.of.recurrence..local..regional.or.distant..or.death.from.breast.cancer)        
}

if(dataset == "Desmedt et al."){
load(file.path(WDIR, "GSE7390_Desmedt_FARMS.rda"))                
prediction.analysis.crossdataset("Desmedt et al.",  t(exprs(an)), pData(an)$t.rfs,
pData(an)$e.rfs)                                
}

if(dataset == "Wang et al."){
load(file.path(WDIR, "Wang_FARMS.rda"))                
prediction.analysis.crossdataset("Wang et al.", t(exprs(an)),
pData(an)$time.to.relapse.or.last.follow.up..months.,
pData(an)$relapse..1.True.)                
}
}
