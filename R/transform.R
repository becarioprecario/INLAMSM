#' @name inla.MCAR.transform
#' @aliases inla.Mmodel.transform
#' @rdname transform
#'
#' @title Transform hyperparameters in multivarite spatial models.
#'
#' @description Multivariate spatial models fit will report hyperparameters
#' in the internal scale. These functions will transform the hyperparameters
#' to a different scale. Using this function requires setting 
#' \code{control.compute = list(config = TRUE)} when fitting the model with
#' INLA.
#' 
#' @param obj An 'inla' object with an MCAR, IMCAR or M-model latent effect.
#' @param k Number of variables in the multivariate model.
#' @param model Either "INDIMCAR", "INDPMCAR", "IMCAR" or "PMCAR". Not used for M-models.
#' @param alpha.min Lower bound of the autocorrelation parameter alpha.
#' @param alpha.max Upper bound of the autocorrelation parameter alpha.
#' 
#' @return This function returns a list with the following elements:
#' \itemize{
#'
#'   \item \code{marginals.hyperpar} List with the posterior marginals of 
#'   transformed hyperparameters.
#'
#'   \item \emph{summary.hyperpar} Summary of the posterior marginals.
#'
#'   \item \emph{VAR.p} Variance matrix of between-variables variability 
#'   computed using point estimates from the posterior marginals.
#'
#'  \item \emph{VAR.m} Posterior mean of variance matrix of the between-variables
#'  variability. This is computed using the internal representation of the
#'  posterior joint distribution of the hyperparameters.
#'
#'  \item \emph{confs} Configurations of the hyperparameters used to
#'  compute \emph{VAR.m}. This is obtained from \code{obj$misc$configs$config}.
#'
#'  \item \emph{M.p} M matrix (only in the M-model) obtained using point
#'  estimates of the parameters.
#'  
#'  \item \emph{M.m} M matrix (only in the M-model) obtained using the
#'  the internal representation of the posterior joint distribution of
#'  the hyperparameters.
#'
#' }
#'
#' @export
## @importFrom INLA inla.tmarginal
## @importFrom INLA inla.zmarginal

inla.MCAR.transform <- function(obj, k, model = "IMCAR", alpha.min, alpha.max) {

  # Check
  if(! model %in% c("INDIMCAR", "INDPMCAR", "IMCAR", "PMCAR")) {
    stop("Parameter 'model' must be one of 'INDIMCAR', 'INDPMCAR', 'IMCAR' or 'PMCAR'.") 
  }

  # Is there a spat. autocorrelation parameter?
  spat.cor.param <- ifelse(model %in% c("INDIMCAR", "IMCAR"), 0, 1)

  # Marginals of diagonal elements in the precision matrix
  # Log-precission to VARIANCES
  margs1 <- lapply(obj$marginals.hyperpar[spat.cor.param + 1:k], function(X) {
    INLA::inla.tmarginal(function(x) exp(-x), X)
    }) 

  # Marginals of between-diseases correlations
  if(model %in% c("IMCAR", "PMCAR")) {
  margs2 <- lapply(obj$marginals.hyperpar[-c(1:(k + spat.cor.param))], function(X) {
    INLA::inla.tmarginal(function(x) {((2 * exp(x))/(1 + exp(x)) - 1)}, X)
    }) 
  }

  if(model %in% c("INDPMCAR", "PMCAR")) {
    margs0 <- INLA::inla.tmarginal(function(x) {
        alpha.min + (alpha.max - alpha.min)/(1 + exp(-x))
      },
      obj$marginals.hyperpar[[1]])

    if(model == "INDPMCAR") {
      margs <- c(list(margs0), margs1)
    } else  {
      margs <- c(list(margs0), margs1, margs2)
    }
  } else {
    if(model == "INDIMCAR") {
      margs <- margs1
    } else {
      margs <- c(margs1, margs2)
    }
  }

  # Summary statistics
  zmargs <- lapply(margs, INLA::inla.zmarginal, silent = TRUE)
  zmargs <- lapply(zmargs, unlist)
  zmargs <- do.call(rbind, zmargs)
  # Add names (this fixes a missing name in the PMCAR model)
  row.names(zmargs) <- names(obj$marginals.hyperpar)

  # Variance covariance matrix (from point estimates)
  if(model %in% c("INDIMCAR", "INDPMCAR")) {
    # No need to coompute this
    VAR.p <- Diagonal(k, x = 1)
    VAR.m <- Diagonal(k, x = 1)
    confs <- NULL
  } else {
    n <- (k - 1) * k / 2
    M <- diag(1, k)
    M[lower.tri(M)] <- zmargs[k + 1:n, "mean"]
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    st.dev <- sqrt(zmargs[1:k, "mean"])
    st.dev.mat <- matrix(st.dev, ncol = 1) %*% matrix(st.dev, 
      nrow = 1)
    VAR.p <- M * st.dev.mat

  
    # Posterior mean of matrix elements using representation of 
    # the latent field

    # Hyperparameters and log.posterior.density
    confs  <- lapply(obj$misc$configs$config, function(X) {
      c(X$theta, X$log.posterior)
    })
    confs <- as.data.frame(do.call(rbind, confs))

    if(model == "PMCAR") 
      confs <- confs[, -c(1:spat.cor.param)] #Remove spatial autocorr. parameter

    # Compute weights
    confs$weight <- confs[, ncol(confs)]
    confs$weight <- exp(confs$weight - max(confs$weight))
    confs$weight <- confs$weight / sum(confs$weight)

    # Transform parameters
    #Log-precisions to variances
    confs[, 1:3] <- exp(-confs[, 1:k])
    # Correlations
    confs[, k + 1:n] <- 2 / (1 + exp(-confs[, k + 1:n])) - 1

    # Posterior mean of matrix entries

    # Compute variacne matrix for all configurations
    aux <- apply(confs[, 1:(k + n)], 1, function(param) {
      M <- diag(1, k)
      M[lower.tri(M)] <- param[-c(1:k)]
      M[upper.tri(M)] <- t(M)[upper.tri(M)]
      st.dev <- sqrt(param[1:k])
      st.dev.mat <- matrix(st.dev, ncol = 1) %*% matrix(st.dev, nrow = 1)
      M <- M * st.dev.mat
      return(M)
    })

    aux <- sapply(1:ncol(aux), function(X) {aux[, X] * confs$weight[X]})
    VAR.m <- matrix(apply(aux, 1, sum), ncol = k)
  }

  return(list(marginals.hyperpar = margs, summary.hyperpar = zmargs, 
    VAR.p = VAR.p, VAR.m = VAR.m, confs = confs))
}

#' @export


inla.Mmodel.transform<- function(obj, k, alpha.min, alpha.max) {
  # Transform autocorrelation parameters to their 'model scale'
  margs1 <- lapply(obj$marginals.hyperpar[1:k], function(m) {
    INLA::inla.tmarginal( function(x) { 
        alpha.min + (alpha.max - alpha.min)/(1 + exp(-x))
      }, m)
  })   

  margs <- c(margs1, obj$marginals.hyperpar[-c(1:k)])

  # Summary statistics
  zmargs <- lapply(margs, INLA::inla.zmarginal, silent = TRUE)
  zmargs <- lapply(zmargs, unlist)
  zmargs <- do.call(rbind, zmargs)

  # M^T M matrix from posterior means (which is so, so...)
  M.p <- matrix(zmargs[-c(1:k), "mean"], ncol = k)
  MtM.p <- crossprod(M.p)


  # M^T M from the diferent configurations
  confs  <- lapply(obj$misc$configs$config, function(X) {
    c(X$theta, X$log.posterior)
  })
  confs <- as.data.frame(do.call(rbind, confs))

  # Compute weights
  confs$weight <- confs[, ncol(confs)]
  confs$weight <- exp(confs$weight - max(confs$weight))
  confs$weight <- confs$weight / sum(confs$weight)

  # Remove autocorr. parameters
  aux <- confs[, -c(1:k)]

  M.aux <- sapply(1:nrow(aux), function(X) {unlist(aux[X, 1:(k * k)] * confs$weight[X])})
  M.m <- matrix(apply(as.matrix(M.aux), 1, sum), ncol = k)

  MtM.aux <- sapply(1:nrow(aux), function(X) {
    crossprod(matrix(unlist(aux[X, 1:(k * k)]), ncol = 3)) * confs$weight[X]
  })
  MtM.m <- matrix(apply(as.matrix(MtM.aux), 1, sum), ncol = k)

  return(list(marginals.hyperpar = margs, summary.hyperpar = zmargs,
    M.p = M.p, VAR.p = MtM.p, M.m = M.m, VAR.m = MtM.m, confs = confs))
}
