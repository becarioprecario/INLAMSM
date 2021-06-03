#' @name inla.rgeneric.Mmodel.model
#' @rdname Mmodel
#'
#' @title M-model: Proper multivariate CAR latent effect with a different
#' spatial autorcorrelation parameter for each disease.
#'
#' @description Multivariate generalization of the proper conditional
#' autorregresive model with one common correlation parameter. This model
#' is performed using the M-model aproximation of Rocamora et. al. (2015).
#'
#' @param cmd Arguments used by latent effects defined using the 'rgeneric'
#' latent effect.
#'
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the 'INLA::inla()'.
#'
#' @details This function is used to define a latent effect that is a
#' multivariate spatial effect based on the M-model aproximation of Rocamora
#' et. al. (2015) in which \eqn{\theta} is modelled as a product of a
#' \eqn{\Phi \cdot M} where the colums of \eqn{\Phi} are modeled independently
#' with a proper conditional autorregresive distribution with a different
#' spatial autocorrelation parameter for each disease and M is a square matrix
#' which introduce de dependence between the diseases. Due to this effect is a
#' multivariate spatial latent effect this function requires the following
#' arguments when defining the latent effect:
#' \itemize{
#'
#'   \item \emph{W} Adjacency SPARSE matrix for spatial effect in the basic
#'   binary code.
#'
#'   \item \emph{k} Number of diseases of the multivariate study.
#'
#'   \item \emph{alpha.min} Minimum value of the spatial autocorrelation
#'   parameter.
#'
#'   \item \emph{alpha.max} Maximum value of the spatial autocorrelation
#'   parameter.
#'
#'}
#'
#' This model is defined using the 'f()' function and an index in order to
#' identify the spatial areas. See the example.
#'
#' @references Palmí-Perales F, Gómez-Rubio V, Martinez-Beneito MA (2021). “Bayesian
#' Multivariate Spatial Models for Lattice Data with INLA.” _Journal of
#' Statistical Software_, *98*(2), 1-29. doi: 10.18637/jss.v098.i02 (URL:
#' https://doi.org/10.18637/jss.v098.i02).
#'
#' @section Prior distributions of the hyperparameters:
#' The hyperparamenters of this lattent effect are the common spatial
#' autocorrelation parameters (one for each disease) and the entries of the
#' M matrix (considered all as a random effects).
#'
#' @examples
#'
#' \donttest{
#' if (require("INLA", quietly = TRUE)) {
#' require(spdep)
#' require(spData)
#' require(rgdal)
#'
#' #Load SIDS data
#' nc.sids <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
#' proj4string(nc.sids) <- CRS("+proj=longlat +ellps=clrk66")
#'
#' #Compute adjacency matrix, as nb object 'adj' and sparse matrix 'W'
#' adj <- poly2nb(nc.sids)
#' W <- as(nb2mat(adj, style = "B"), "Matrix")
#'
#' #Compute expected cases
#' r74 <- sum(nc.sids$SID74) / sum(nc.sids$BIR74)
#' nc.sids$EXP74 <- r74 * nc.sids$BIR74
#' nc.sids$SMR74 <- nc.sids$SID74 / nc.sids$EXP74
#' nc.sids$NWPROP74 <- nc.sids$NWBIR74 / nc.sids$BIR74
#'
#' r79 <- sum(nc.sids$SID79) / sum(nc.sids$BIR79)
#' nc.sids$EXP79 <- r79 * nc.sids$BIR79
#' nc.sids$SMR79 <- nc.sids$SID79 / nc.sids$EXP79
#' nc.sids$NWPROP79 <- nc.sids$NWBIR79 / nc.sids$BIR79
#'
#' # Data (replicated to assess scalability)
#'
#' #Real data
#' n.rep <- 1
#' d <- list(OBS = c(nc.sids$SID74, nc.sids$SID79),
#'           NWPROP = c(nc.sids$NWPROP74, nc.sids$NWPROP79),
#'           EXP = c(nc.sids$EXP74, nc.sids$EXP79))
#' d <- lapply(d, function(X) { rep(X, n.rep)})
#' d$idx <- 1:length(d$OBS)
#'
#' #Parameters of the Mmodel
#' k <- 2
#' alpha.min <- 0
#' alpha.max <- 1
#'
#'
#'
#' model <- inla.rgeneric.define(inla.rgeneric.Mmodel.model, debug = FALSE,
#'                               k = k, W = W, alpha.min = alpha.min,
#'                               alpha.max = alpha.max)
#'
#' r.Mmodel <- inla(OBS ~ -1 + f(idx, model = model), data = d, E = EXP,
#'   family = "poisson", control.predictor = list(compute = TRUE))
#'
#' nc.sids$Model1 <- r.Mmodel$summary.random$idx[1:100, "mean"]
#' nc.sids$Model2 <- r.Mmodel$summary.random$idx[100 + 1:100, "mean"]
#'
#' spplot(nc.sids, c("Model1", "Model2"))
#'
#' nc.sids$Fit1 <- r.Mmodel$summary.fitted[1:100, "mean"]
#' nc.sids$Fit2 <- r.Mmodel$summary.fitted[100 + 1:100, "mean"]
#'
#' spplot(nc.sids, c("Fit1", "SMR74", "Fit2", "SMR79"))
#'
#'
#' ## Running UNIVARIATE MODEL
#'
#' #Real data
#' n.rep <- 1
#' d <- list(OBS = nc.sids$SID74,
#'           NWPROP = nc.sids$NWPROP74,
#'           EXP = nc.sids$EXP74)
#' d <- lapply(d, function(X) { rep(X, n.rep)})
#' d$idx <- 1:length(d$OBS)
#'
#' #Fit model
#' r.uni <- inla(OBS ~ 1 + f(idx, model = "besag", graph = W), # + NWPROP,
#'               data = d, E = EXP, family = "poisson",
#'               control.predictor = list(compute = TRUE))
#'
#' summary(r.uni)
#'
#' nc.sids$FITTED74.uni <- r.uni$summary.fitted.values[ , "mean"]
#'
#' #Display univariate VS multivariate  fitted relative risks.
#' dev.new()
#' spplot(nc.sids, c("SMR74", "Fit1", "FITTED74.uni"))
#' spplot(nc.sids, c("Fit1", "FITTED74.uni"),
#'        main=list(label="Relative risk estimation",cex=2))
#' dev.new()
#' plot(nc.sids$FITTED74.uni, nc.sids$Fit1, main="Relative Risk estimations",
#'      xlab="Univariate RR estimations"
#'      , ylab="Multivariate RR estimations")#, xlim=c(0.5, 2.5), ylim=c(0, 2))
#' abline(h=0, col="grey")
#' abline(v=0, col="grey")
#' abline(a=0, b=1, col="red")
#'
#' }
#' }
#'
#' @export
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix bdiag
#' @importFrom MCMCpack dwish
#' @usage inla.rgeneric.Mmodel.model(cmd, theta)

# Define previous variables as global to avoid warnings()
utils::globalVariables(c("k", "W", "alpha.min", "alpha.max"))

'inla.rgeneric.Mmodel.model' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
{
  ## M-model implementation using a proper CAR with different parameters
  ## k: number of diseases/blocks
  ## W: adjacency matrix
  ## alpha.min: minimum value for alpha
  ## alpha.max: maximum value for alpha


  #theta: k spatial correlation parameter alpha, k * k entries of the M matrix,
  #  by columns.
  # NOTE: The CAR distributions do NOT have a precision parameter.
  interpret.theta <- function()
  {
    #Function for changing from internal scale to external scale
    # also, build the inverse of the matrix used to model in the external scale
    # the between-disease variability and then this matrix is inverted.

    # First k parameters are the autocorrelation parameter
    alpha <- alpha.min + (alpha.max - alpha.min) /
      (1 + exp(-theta[as.integer(1:k)]))

    # The next k * k parameters are the entries in the M matrix, by cols.
    M <- matrix(theta[-as.integer(1:k)], ncol = k)

    return (list(alpha = alpha, M = M))
  }


  #Graph of precision function; i.e., a 0/1 representation of precision matrix
  graph <- function()
  {

    # M \kronecker I
    MI <- kronecker(Matrix::Matrix(1, ncol = k, nrow = k),
      Matrix::Diagonal(nrow(W), 1))

    # I + W
    IW <- Matrix::Diagonal(nrow(W), 1) + W

    # Block diagonal
    BlockIW <- Matrix::bdiag(replicate(k, IW, simplify = FALSE))

    G <- (MI %*% BlockIW) %*% MI
    return (G)
  }

  #Precision matrix
  Q <- function()
  {
    #Parameters in model scale
    param <- interpret.theta()

    # M^{-1} \kronecker I
    M.inv <- solve(param$M)
    MI <- kronecker(M.inv, Diagonal(nrow(W), 1))

    # Number of neighbours
    D <- as.vector(apply(W, 1, sum))

    BlockIW <- Matrix::bdiag(lapply(1:k, function(i) {
      Matrix::Diagonal(x = D) - param$alpha[i] * W
    }))

    Q <- (MI %*% BlockIW) %*% kronecker(t(M.inv), Matrix::Diagonal(nrow(W), 1))

    return (Q)
  }

  #Mean of model
  mu <- function() {
    return(numeric(0))
  }

  log.norm.const <- function() {
    ## return the log(normalising constant) for the model
    #param = interpret.theta()
    #
    #val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) +
    # 0.5 * log(1.0 - param$alpha^2)

    val <- numeric(0)
    return (val)
  }

  log.prior <- function() {
    ## return the log-prior for the hyperparameters.
    ## Uniform prior in (alpha.min, alpha.max) on model scale
    param <- interpret.theta()

    # log-Prior for the autocorrelation parameter
    val <-  sum(-theta[as.integer(1:k)]
               - 2 * log(1 + exp(-theta[as.integer(1:k)]))
               )

    # Whishart prior for M^T * M
    sigma2 <- 1000 #Wishart parameter
    val = val + log(MCMCpack::dwish(W = crossprod(param$M),
      v = k, S = diag(rep(sigma2, k))))

    return (val)
  }

  initial <- function() {
    ## return initial values: k spat. autoc. and M-matrix columns

    return ( c(rep(0, k), as.vector(diag(rep(1, k)))))

  }

  quit <- function() {
    return (invisible())
  }

  # FIX for rgeneric to work on R >= 4
  # Provided by E. T. Krainski
  if (as.integer(R.version$major) > 3) {
    if (!length(theta))
      theta = initial()
  } else {
    if (is.null(theta)) {
      theta <- initial()
    }
  }

  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

##' @rdname Mmodel
##' @usage inla.Mmodel.model(...)
##' @param ...  Arguments to be passed to 'inla.rgeneric.define'.
##' @export
inla.Mmodel.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.Mmodel.model, ...)
}

