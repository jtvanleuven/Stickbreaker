#' Simulate data from priors then use to calcualte posterior probability of models given data
#'
#' @param fit.smry X
#' @param data X
#' @param analysis.name X
#' @param coes.prior X
#' @param sig.prior X
#' @param d.range X
#' @param d.adj.max X
#' @param wts X
#' @param n.samps.per.mod X
#' @return List:\cr
#' \code{$posteriors} Posterior probability of each epistasis model. \cr
#' \code{$multinomial.model} Multinomial model fit to the simulated data used to assign posterior probabilities.
#' @export

simulate.data.calculate.posteriors <- function(fit.smry, data, analysis.name, coes.prior, sig.prior, d.range, d.adj.max, wts, n.samps.per.mod){

  n.genos <- length(data[,1])
  n.muts <- length(data[1,])-1
  geno.matrix <- data[,seq(1, n.muts)]
  fit.matrix <- as.matrix(data[,(n.muts+1)])

  outdir <- system.file("extdata", package="Stickbreaker")
  simdata.file.name <- paste("simulated_data_from_priors_", analysis.name, ".txt", sep="")
  simdata.outpath <- paste(outdir, simdata.file.name, sep="/")

  simulate.partial.data.from.priors.for.mod.selection(geno.matrix,
                                                      coes.prior=coes.prior,
                                                      sigs.prior=sig.prior,
                                                      mods.to.sim=c("stick", "mult", "add"),
                                                      d.true=1,
                                                      d.range=d.range,
                                                      d.adj.max=d.adj.max,
                                                      w.wt=1,
                                                      wts=wts,
                                                      outpath=simdata.outpath,
                                                      n.samps.per.mod=n.samps.per.mod,
                                                      coe.sim.model="identical",
                                                      coe.dist.par=NA)

  print("Fitting multinomial model to simulated data")
  simdata <- read.table(file=simdata.outpath, header=TRUE)
  data$R2.stick[which(data$R2.stick<0)] <- -5
  data$R2.mult[which(data$R2.mult<0)] <- -5
  data$R2.add[which(data$R2.add<0)] <- -5
  mod.formula <- as.formula(model ~ R2.stick + R2.mult + R2.add + P.stick + P.mult + P.add)
  mod <- fit.nnet.multinomial.regression(simdata, mod.formula)
  posteriors <- predict(mod, newdata=fit.smry, type="probs")
  return(list(posteriors=posteriors, multinomial.model=mod))
}


#' Simulate data at specified parameter values for doing model selection
#'
#' @param n.muts Number of mutations to simulate. Vector of the form \code{c(3,4,5)}.
#' @param coes.to.sim Coefficient values to simulate. Vector of the form \code{c(0.1, 0.3, 0.5)}.
#' @param sigs.to.sim Sigma values to simulate. Vector of the form \code{c(0.02, 0.05, 0.08)}.
#' @param mods.to.sim Models to simulate under. Vector with model names of the form c("stick", "mult", "add").
#' @param d.true True value of d.
#' @param d.range Range of possible values for d. If estimate is outside this, estimate is not considered valid.
#' @param w.wt Wild type fitness.
#' @param wts Weights when estimating coefficients and d.
#' @param outpath Full path including file name to write results
#' @param n.reps.ea Number of simulated-fit datasets per parametric condition
#' @param coe.sim.model Coefficient simulation model. Possible values: "identical" means all coefficients take same value-- from coe.v.
#' "uniform" indicates to  sample individual coefficients from a uniform with mean given by coe.v and distributed U(coe.v-coe.dist.par, coe.v+coe.dist.parm).
#' "normal" means sample coefficients from normal distribuiton with mean from coe.v and sigma given by coe.dist.par.
#' Default = "identical".
#' @param coe.dist.par Coefficient distribution parameter. If coe.sim.model=="uniform", then uniform is U(-coe.dist.par, +coe.dist.parm).
#' If coe.sim.model=="normal", then distributed normal with mean given by coe.v and sigma give by coe.dist.parm.
#' @return Nothing. Instead results are written to \code{outpath} file for later analysis.
#' @details Function loops over all parametric combinations and simulates datasets. For each dataset it
#' fits to each of the models and outputs a row of metrics that summarize the fits.
#' @export

simulate.data.for.mod.selection <- function(n.muts, coes.to.sim, sigs.to.sim, mods.to.sim, d.true, d.range, w.wt, wts, outpath, n.reps.ea, coe.sim.model="identical", coe.dist.par=NA){
  geno.matrix <- generate.geno.matrix(n.muts)
  first.results <- TRUE
  for (mod.i in 1:length(mods.to.sim)){
    model <- mods.to.sim[mod.i]
    for (coe.i in 1:length(coes.to.sim)){
      if (coe.sim.model=="identical"){
        coes <- rep(coes.to.sim[coe.i], n.muts)
      } else if (coe.sim.model == "uniform"){
        coes <- runif(n.muts, coes.to.sim[coe.i]-coe.dist.par, coes.to.sim[coe.i]+coe.dist.par)
      } else if (coe.sim.model == "normal"){
        coes <- rnorm(n=n.muts, mean=coes.to.sim[coe.i], sd=coe.dist.par)
      }

      for (sig.i in 1:length(sigs.to.sim)){
        sig.val <- sigs.to.sim[sig.i]

        for (rep.i in 1:n.reps.ea){
          if (model=="stick"){
            fit.matrix <- simulate.stick.data(n.muts=n.muts, coes=coes, sigma=sig.val, d.true=d.true, w.wt=w.wt)$fit.matrix
          } else if (model=="mult"){
            fit.matrix <- simulate.mult.data(n.muts=n.muts, selcoes=coes, sigma=sig.val, w.wt=w.wt)$fit.matrix
          } else if (model=="add"){
            fit.matrix <- simulate.add.data(n.muts=n.muts, addcoes=coes, sigma=sig.val, w.wt=w.wt)$fit.matrix
          }

          # --- Fit data to each model ---
            # -- stickbreaking --
              d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
              d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
              d.hat.max <- max(fit.matrix, na.rm=TRUE) - w.wt
              d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range)
              fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, wts=wts, run.regression=TRUE)
              u.hat.mean <- mean(fit.stick$u.hats)
              fit.stick <- c(u.mean=u.hat.mean, d.hat=d.hat.seq, round(unlist(fit.stick[1:4]),5), P=round(fit.stick$regression.results$P,5))
              oldnames <- c("R2", "sig.hat", "logLike", "P")    # rename to include model
              newnames <- c("R2.stick", "sig.stick", "lnL.stick", "P.stick")
              names(fit.stick)[sapply(oldnames, function(x) which(names(fit.stick)==x))] <- newnames
            # -- multiplicative --
              fit.mult <-fit.mult.model(geno.matrix, fit.matrix, wts)
              s.hat.mean <- mean(fit.mult$s.hats)
              fit.mult <- c(s.mean=s.hat.mean, round(unlist(fit.mult[1:4]),5), P=round(fit.mult$regression.results$P,5))
              newnames <- c("R2.mult", "sig.mult", "lnL.mult", "P.mult")
              names(fit.mult)[sapply(oldnames, function(x) which(names(fit.mult)==x))] <- newnames
            # -- additive --
              fit.add <- fit.add.model(geno.matrix, fit.matrix, wts)
              w.hat.mean <- mean(fit.add$w.hats)
              fit.add <- c(w.mean=w.hat.mean, round(unlist(fit.add[1:4]),5), P=round(fit.add$regression.results$P,5))
              newnames <- c("R2.add", "sig.add", "lnL.add", "P.add")
              names(fit.add)[sapply(oldnames, function(x) which(names(fit.add)==x))] <- newnames

              parm.row <- data.frame(model=model, n.muts=n.muts, coes=coes.to.sim[coe.i], sigma=sig.val)
              fit.row <- data.frame(cbind(t(fit.stick), t(fit.mult)), t(fit.add))

              one.row <- data.frame(cbind(parm.row, fit.row))

              if (first.results == TRUE){
                write.table(x=one.row, file=outpath, append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
                first.results <- FALSE
              } else{
                write.table(x=one.row, file=outpath, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
              }
        } #next rep.i
      }  #next sig.i
    }  #next coe.i
  } #next model.i
}





#' Simulate data from priors for doing model selection
#'
#' @param n.muts Number of mutations to simulate. Vector of the form \code{c(3,4,5)}.
#' @param coes.prior Range for uniform prior on coefficients
#' @param sigs.prior Range for uniform prior on sigma
#' @param mods.to.sim Models to simulate under. Vector with model names of the form c("stick", "mult", "add").
#' @param d.true True value of d.
#' @param d.range Range of possible values for d. If estimate is outside this, estimate is not considered valid.
#' @param d.adj.max Factor to increase observed distanct to max fitness by for ad hoc d estimate (when other estimators fail)
#' @param w.wt Wild type fitness.
#' @param wts Weights when estimating coefficients and d.
#' @param outpath Full path including file name to write results
#' @param n.samps.per.mod Number of datasets to simulate per model
#' @param coe.sim.model Coefficient simulation model. See details.
#' @param coe.dist.par Coefficient distribution parameter. If coe.sim.model=="uniform", then uniform is U(-coe.dist.par, +coe.dist.parm).
#' If coe.sim.model=="normal", then distributed normal with mean given by coe.v and sigma give by coe.dist.parm.
#' @param print.interval Every this many replicates, prints out replicate number. If NA (default) no printing is done.
#' @return Nothing. Instead results are written to \code{outpath} file for later analysis
#' @details This function generates datasets by drawing from priors. It generates \code{n.samps.per.mod}
#' per model. It then analyzes each dataset under all three models writes one row of summary statistics
#' to the output file (defined by \code{outpath}).
#' \code{coe.sim.model}: The expected coefficient is sampled from uniform prior (coes.prior): E[coe]
#' The coe.sim.model determines how the individual coefficients are generated.
#'  Possible values: "identical" means all coefficients take same value--E[coe].
#' "uniform" indicates to  sample individual coefficients from a uniform distribution: U(E[coe]-coe.dist.par, E[coe]+coe.dist.parm).
#' "normal" means sample coefficients from normal distribuiton with mean E[coe] and sigma given by coe.dist.par.
#' Default = "identical".
#' @export

simulate.data.from.priors.for.mod.selection <- function(n.muts, coes.prior, sigs.prior, mods.to.sim, d.true, d.range, d.adj.max, w.wt, wts, outpath, n.samps.per.mod, coe.sim.model="identical", coe.dist.par=NA, print.interval=NA){
  geno.matrix <- generate.geno.matrix(n.muts)
  if (is.na(print.interval)==FALSE){
    print.pts <- c(1, seq(print.interval, n.samps.per.mod, by=print.interval))
  }
  first.results <- TRUE
  for (mod.i in 1:length(mods.to.sim)){
    model <- mods.to.sim[mod.i]
    if (is.na(print.interval)==FALSE){
      print(model)
    }
    for (rep.i in 1:n.samps.per.mod){

      if (is.na(print.interval)==FALSE){
        if (rep.i %in% print.pts){
          print(rep.i)
        }
      }
      if (coe.sim.model=="identical"){
        coes.prior <- sort(coes.prior)
        coe.center <- runif(n=1,coes.prior[1], coes.prior[2])
        coes <- rep(coe.center, n.muts)
      } else if (coe.sim.model == "uniform"){
        coes <- runif(n.muts, coe.center-coe.dist.par, coe.center+coe.dist.par)
      } else if (coe.sim.model == "normal"){
        coes <- rnorm(n=n.muts, mean=coe.center, sd=coe.dist.par)
      }

      sig.val <- runif(n=1, sigs.prior[1], sigs.prior[2])

      if (model=="stick"){
        fit.matrix <- simulate.stick.data(n.muts=n.muts, coes=coes, sigma=sig.val, d.true=d.true, w.wt=w.wt, geno.matrix)$fit.matrix
      } else if (model=="mult"){
        fit.matrix <- simulate.mult.data(n.muts=n.muts, selcoes=coes, sigma=sig.val, w.wt=w.wt, geno.matrix)$fit.matrix
      } else if (model=="add"){
        fit.matrix <- simulate.add.data(n.muts=n.muts, addcoes=coes, sigma=sig.val, w.wt=w.wt, geno.matrix)$fit.matrix
      }

      # --- Fit data to each model ---
      # -- stickbreaking --
        d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
        d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
        d.hat.max <- max(fit.matrix, na.rm=TRUE) - w.wt
        d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range, d.adj.max)
        fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, wts=wts,run.regression=TRUE)
        u.hat.mean <- mean(fit.stick$u.hats)
        fit.stick <- c(u.mean=u.hat.mean, d.hat=d.hat.seq, round(unlist(fit.stick[1:4]),5), P=round(fit.stick$regression.results$P,5))
        oldnames <- c("R2", "sig.hat", "logLike", "P")    # rename to include model
        newnames <- c("R2.stick", "sig.stick", "lnL.stick", "P.stick")
        names(fit.stick)[sapply(oldnames, function(x) which(names(fit.stick)==x))] <- newnames
        # -- multiplicative --
        fit.mult <-fit.mult.model(geno.matrix, fit.matrix, wts)
        s.hat.mean <- mean(fit.mult$s.hats)
        fit.mult <- c(s.mean=s.hat.mean, round(unlist(fit.mult[1:4]),5), P=round(fit.mult$regression.results$P,5))
        newnames <- c("R2.mult", "sig.mult", "lnL.mult", "P.mult")
        names(fit.mult)[sapply(oldnames, function(x) which(names(fit.mult)==x))] <- newnames
        # -- additive --
        fit.add <- fit.add.model(geno.matrix, fit.matrix, wts)
        w.hat.mean <- mean(fit.add$w.hats)
        fit.add <- c(w.mean=w.hat.mean, round(unlist(fit.add[1:4]),5), P=round(fit.add$regression.results$P,5))
        newnames <- c("R2.add", "sig.add", "lnL.add", "P.add")
        names(fit.add)[sapply(oldnames, function(x) which(names(fit.add)==x))] <- newnames

        parm.row <- data.frame(model=model, n.muts=n.muts, coes=coe.center, sigma=sig.val)
        fit.row <- data.frame(cbind(t(fit.stick), t(fit.mult)), t(fit.add))

        one.row <- data.frame(cbind(parm.row, fit.row))

        if (first.results == TRUE){
          write.table(x=one.row, file=outpath, append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
          first.results <- FALSE
        } else{
          write.table(x=one.row, file=outpath, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
        }
      } #next rep.i
  } #next model.i
}


#' Fit training data to multinomial regression using \code{nnet} package
#'
#'  @param data Dataset to fit
#'  @export
#'
fit.nnet.multinomial.regression <- function(data, mod.formula){
  mnet.fit <- multinom(data=data, formula = mod.formula)
  return(mnet.fit)
}



#' Calculates posterior probabilities for each row of dataset given model
#'
#' @param data Dataset as generated by \code{\link{simulate.data.from.priors.for.mod.selection}}
#' or \code{\link{simulate.data.for.mod.selection}}
#' @param model Multinomial regression model from nnet package fit using the \code{\link{fit.nnet.multinomial.regression}} function.
#' @return dataset with the posterior probabilities of add, mult and stick models bound (using cbind)
#' to right edge of the supplied dataset.
#' @export

calculate.posteriors.for.datasets <- function(data, mod){

  posteriors <- as.data.frame(matrix(nrow=length(data[,1]), ncol=3))
  colnames(posteriors)=c("add", "mult", "stick")

  for (j in 1:length(data[,1])){
    posteriors[j,] <- predict(mod, newdata=data[j,], type="probs")
  }
  dataset <- cbind(data, posteriors)
  return(dataset)
}







#'  Calculate classification performance on simulated data
#'
#' @param data Dataset as generated by \code{\link{simulate.data.from.priors.for.mod.selection}}
#' or \code{\link{simulate.data.for.mod.selection}}
#' @param rej.cut Rejection cutoff point specifies at what posterior probability is a model rejected. Default is 0.05.
#' @return Dataframe with false rejection rate, correct and unique classification rate, and
#'  mean posterior probability for each parameter combination in the data.
#' @details If a model has posterior probabiltiy < \code{rej.cut}, it is rejected. If a model has posterior
#' probabability \eqn{\ge}  1-\code{rej.cut}, it is considered uniquely accepted.
#' @export

summarize.posteriors.on.simulated.dataset <- function(data, rej.cut=0.05){
  mods.to.lyze <- unique(data$model)
     coes.to.lyze <- unique(data$coes)
     sigs.to.lyze <- unique(data$sigma)
     parm.table <- expand.grid(model=mods.to.lyze, coes=coes.to.lyze, sigma=sigs.to.lyze)
     parm.table <- cbind(model=parm.table[,1], n.muts=rep(n.muts, length(parm.table[,1])), parm.table[,2:3])
     parm.strings <- apply(parm.table, 1, function(x) paste(x, collapse="_"))
     end.in.zero <- which(sapply(parm.strings, function(x) substring(x, nchar(x))=="0"))
     if (length(end.in.zero)>0){
        parm.strings[end.in.zero] <- sapply(end.in.zero, function(x) substring(parm.strings[x], 1, nchar(parm.strings[x])-1))
     }

     result.cols <- c("n", "false.rej", "true.unq", "mean.post")
     results.table <- as.data.frame(matrix(nrow=length(parm.table[,1]), ncol=length(result.cols)))
     colnames(results.table) <- result.cols
     one.table <- cbind(parm.table, results.table)
     mod.cols <- sapply(c("add", "mult", "stick"), function(x) which(colnames(data)==x))
     for (mod.i in 1:length(mods.to.lyze)){
       for (coe.i in 1:length(coes.to.lyze)){
          for (sig.i in 1:length(sigs.to.lyze)){
              rows <- with(data, which(model==mods.to.lyze[mod.i] & n.muts==muts.to.lyze[mut.i] &
                                         coes==coes.to.lyze[coe.i] & sigma==sigs.to.lyze[sig.i]))
              pattern <- paste(as.character(data[rows[1],1:4]), collapse="_")
              nc <- nchar(pattern)
              if (substring(pattern, nc-1, nc)=="_0"){pattern <- paste(substring(pattern, 1, nc-2), "_0.0", sep="")}
              parm.row <- which(parm.strings==pattern)
              n <- length(rows)
              one.table$n[parm.row] <- n
              true.mod.col <- mod.cols[which(names(mod.cols)==mods.to.lyze[mod.i])]
              true.mod.probs <- data[rows, true.mod.col]
              one.table$false.rej[parm.row] <- length(which(true.mod.probs <= rej.cut))/n
              one.table$true.unq[parm.row] <- length(which(true.mod.probs >= 1-rej.cut))/n
              one.table$mean.post[parm.row] <- mean(true.mod.probs)
          } #next sig.i
        }  #next coe.I
     }  #next mod.i
     return(one.table)
}







#' Simulate partial data from priors for doing model selection
#'
#' @param geno.matrix Matrix specifying which genotypes to simulate
#' @param coes.prior Range for uniform prior on coefficients
#' @param sigs.prior Range for uniform prior on sigma
#' @param mods.to.sim Models to simulate under. Vector with model names of the form c("stick", "mult", "add").
#' @param d.true True value of d.
#' @param d.range Range of possible values for d. If estimate is outside this, estimate is not considered valid.
#' @param d.adj.max Factor to increase observed distanct to max fitness by for ad hoc d estimate (when other estimators fail)
#' @param w.wt Wild type fitness.
#' @param wts Weights when estimating coefficients and d.
#' @param outpath Full path including file name to write results
#' @param n.samps.per.mod Number of datasets to simulate per model
#' @param coe.sim.model Coefficient simulation model. See details.
#' @param coe.dist.par Coefficient distribution parameter. If coe.sim.model=="uniform", then uniform is U(-coe.dist.par, +coe.dist.parm).
#' If coe.sim.model=="normal", then distributed normal with mean given by coe.v and sigma give by coe.dist.parm.
#' @return Nothing. Instead results are written to \code{outpath} file for later analysis
#' @details This function generates datasets by drawing from priors. It generates \code{n.samps.per.mod}
#' per model. It then analyzes each dataset under all three models writes one row of summary statistics
#' to the output file (defined by \code{outpath}).
#' \code{coe.sim.model}: The expected coefficient is sampled from uniform prior (coes.prior): E[coe]
#' The coe.sim.model determines how the individual coefficients are generated.
#'  Possible values: "identical" means all coefficients take same value--E[coe].
#' "uniform" indicates to  sample individual coefficients from a uniform distribution: U(E[coe]-coe.dist.par, E[coe]+coe.dist.parm).
#' "normal" means sample coefficients from normal distribuiton with mean E[coe] and sigma given by coe.dist.par.
#' Default = "identical".
#' @export

simulate.partial.data.from.priors.for.mod.selection <- function(geno.matrix, coes.prior, sigs.prior, mods.to.sim, d.true, d.range, d.adj.max, w.wt, wts, outpath, n.samps.per.mod, coe.sim.model="identical", coe.dist.par=NA){
  n.muts <- length(geno.matrix[1,])
  first.results <- TRUE
  print.interval <- 25
  print(paste("Simulating data from priors. ", n.samps.per.mod, " reps per mod."))
  print("Model, replicate #:")
  if (n.samps.per.mod >= print.interval){
    print.v <- c(1,seq(print.interval, n.samps.per.mod, by=print.interval))
  } else{
    print.v <- 1
  }
  for (mod.i in 1:length(mods.to.sim)){
    model <- mods.to.sim[mod.i]
    print(model)
    for (rep.i in 1:n.samps.per.mod){
        if (rep.i %in% print.v){
          print(rep.i)
        }
      if (coe.sim.model=="identical"){
        coes.prior <- sort(coes.prior)
        coe.center <- runif(n=1,coes.prior[1], coes.prior[2])
        coes <- rep(coe.center, n.muts)
      } else if (coe.sim.model == "uniform"){
        coes <- runif(n.muts, coe.center-coe.dist.par, coe.center+coe.dist.par)
      } else if (coe.sim.model == "normal"){
        coes <- rnorm(n=n.muts, mean=coe.center, sd=coe.dist.par)
      }

      sig.val <- runif(n=1, sigs.prior[1], sigs.prior[2])

      if (model=="stick"){
        fit.matrix <- simulate.stick.data(n.muts=n.muts, coes=coes, sigma=sig.val, d.true=d.true, w.wt=w.wt, geno.matrix)$fit.matrix
      } else if (model=="mult"){
        fit.matrix <- simulate.mult.data(n.muts=n.muts, selcoes=coes, sigma=sig.val, w.wt=w.wt, geno.matrix)$fit.matrix
      } else if (model=="add"){
        fit.matrix <- simulate.add.data(n.muts=n.muts, addcoes=coes, sigma=sig.val, w.wt=w.wt, geno.matrix)$fit.matrix
      }

      # --- Fit data to each model ---
      # -- stickbreaking --
      d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
      d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
      d.hat.max <- max(fit.matrix, na.rm=TRUE) - w.wt
      d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range, d.adj.max)
      fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, wts=wts,run.regression=TRUE)
      u.hat.mean <- mean(fit.stick$u.hats)
      fit.stick <- c(u.mean=u.hat.mean, d.hat=d.hat.seq, round(unlist(fit.stick[1:4]),5), P=round(fit.stick$regression.results$P,5))
      oldnames <- c("R2", "sig.hat", "logLike", "P")    # rename to include model
      newnames <- c("R2.stick", "sig.stick", "lnL.stick", "P.stick")
      names(fit.stick)[sapply(oldnames, function(x) which(names(fit.stick)==x))] <- newnames
      # -- multiplicative --
      fit.mult <-fit.mult.model(geno.matrix, fit.matrix, wts)
      s.hat.mean <- mean(fit.mult$s.hats)
      fit.mult <- c(s.mean=s.hat.mean, round(unlist(fit.mult[1:4]),5), P=round(fit.mult$regression.results$P,5))
      newnames <- c("R2.mult", "sig.mult", "lnL.mult", "P.mult")
      names(fit.mult)[sapply(oldnames, function(x) which(names(fit.mult)==x))] <- newnames
      # -- additive --
      fit.add <- fit.add.model(geno.matrix, fit.matrix, wts)
      w.hat.mean <- mean(fit.add$w.hats)
      fit.add <- c(w.mean=w.hat.mean, round(unlist(fit.add[1:4]),5), P=round(fit.add$regression.results$P,5))
      newnames <- c("R2.add", "sig.add", "lnL.add", "P.add")
      names(fit.add)[sapply(oldnames, function(x) which(names(fit.add)==x))] <- newnames

      parm.row <- data.frame(model=model, n.muts=n.muts, coes=coe.center, sigma=sig.val)
      fit.row <- data.frame(cbind(t(fit.stick), t(fit.mult)), t(fit.add))

      one.row <- data.frame(cbind(parm.row, fit.row))

      if (first.results == TRUE){
        write.table(x=one.row, file=outpath, append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
        first.results <- FALSE
      } else{
        write.table(x=one.row, file=outpath, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
      }
    } #next rep.i
  } #next model.i
}
