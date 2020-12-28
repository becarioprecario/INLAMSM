#' @name inla.rgeneric.IMCAR.model
#' @rdname imcar
#'
#' @title \eqn{IMCAR(\Lambda)}: Intrinsic multivariate CAR latent effect.
#'
#' @description Multivariate generalization of the intrinsic conditional
#' autorregresive model. The matrix which models the variability between
#' diseases is a symmetric matrix with the inverse of the marginal
#' precisions on the diagonal elements and the correlation parameters
#' divided by the square root of the precisions on the off-diagonal elements.
#'
#' @param cmd Arguments used by latent effects defined using the 'rgeneric'
#' latent effect.
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the 'INLA::inla()'.
#'
#' @details This function is used to define a latent effect that is a
#' multivariate spatial effect with a intrinsic conditional autorregresive
#' distribution and a symmetric matrix in order to model the whitin-disease and
#' the between-diseases variability, respectively. Due to this effect is a
#' multivariate spatial latent effect this function requires the following
#' arguments when defining the latent effect:
#' \itemize{
#'
#'   \item \emph{W} Adjacency SPARSE matrix for spatial effect in the basic
#'   binary code.
#'
#'   \item \emph{k} Number of diseases of the multivariate study.
#'
#'}
#'
#' This model is defined using the 'f()' function and an index in order to
#' identify the spatial areas. See the example.
#'
#' @section Prior distributions of the hyperparameters:
#' The hyperparamenters of this lattent effect are the marginal precisions of
#' each disease which are equal to the number of diseases and the correlation
#' parameters for the whole pair of diseases.
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
#' # Define sum-to-zero constraints
#' A <- kronecker(Diagonal(2, 1), Matrix(1, ncol = nrow(W), nrow = 1))
#' e  = rep(0, 2)
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
#' # Model parameters
#' k <- 2 * n.rep #Number of diseases
#'
#'
#' #Define model IMCAR
#' model <- inla.rgeneric.define(inla.rgeneric.IMCAR.model, debug = FALSE,
#'   k = k, W = W)
#'
#' #Fit model
#' r <- inla(OBS ~ 1 + f(idx, model = model,
#'     extraconstr = list(A = as.matrix(A), e = e)), # + NWPROP,
#'   data = d, E = EXP, family = "poisson",
#'   control.compute = list(config = TRUE), 
#'   control.predictor = list(compute = TRUE))
#'
#' summary(r)
#'
#' # Transformed parameters
#' r.hyperpar <- inla.MCAR.transform(r, k = 2, model = "IMCAR")
#' r.hyperpar$summary.hyperpar
#' 
#' #Get fitted data, i.e., relative risk
#' nc.sids$FITTED74 <- r$summary.fitted.values[1:100, "mean"]
#' nc.sids$FITTED79 <- r$summary.fitted.values[100 + 1:100, "mean"]
#'
#' #Display fitted relative risks
#' dev.new()
#' spplot(nc.sids, c("SMR74", "FITTED74", "SMR79", "FITTED79"))
#'
#' #Show marginals of tau_1, tau_2, rho
#'
#' marg.tau1 <- inla.tmarginal(
#'   function(x) exp(x),
#'   r$marginals.hyperpar[[1]])
#'
#' marg.tau2 <- inla.tmarginal(
#'   function(x) exp(x),
#'   r$marginals.hyperpar[[2]])
#'
#' marg.rho <- inla.tmarginal(
#'   function(x) (2*exp(x))/(1 + exp(x)) - 1,
#'   r$marginals.hyperpar[[3]])
#'
#'
#' dev.new()
#'
#' oldpar <- par(mfrow = c(2, 2))
#' plot(marg.tau1, main = "tau1", type = "l")
#' plot(marg.tau2, main = "tau2", type = "l")
#' plot(marg.rho, main = "rho", type = "l")
#'
#' par(oldpar)
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
#' spplot(nc.sids, c("SMR74", "FITTED74", "FITTED74.uni"))
#' spplot(nc.sids, c("FITTED74", "FITTED74.uni"),
#'        main=list(label="Relative risk estimation",cex=2))
#' dev.new()
#' plot(nc.sids$FITTED74.uni, nc.sids$FITTED74,
#'      main="Relative Risk estimations", xlab="Univariate RR estimations",
#'      ylab="Multivariate RR estimations", xlim=c(0.5, 2.5), ylim=c(0.5, 2.5))
#' abline(h=0, col="grey")
#' abline(v=0, col="grey")
#' abline(a=0, b=1, col="red")
#'
#' #Plot posterior mean of the spatial effects univ VS multi
#'
#' nc.sids$m.uni <- r.uni$summary.random$idx[, "mean"]
#' nc.sids$m.mult <- r$summary.random$idx[1:100, "mean"]
#' dev.new()
#' plot(nc.sids$m.uni, nc.sids$m.mult,
#'      main="Posterior mean of the spatial effect", xlab="Uni. post. means"
#'      , ylab="Mult. post. means", xlim=c(-1,1), ylim=c(-1,1))
#' abline(h=0, col="grey")
#' abline(v=0, col="grey")
#' abline(a=0, b=1, col="red")
#'
#' dev.new()
#' spplot(nc.sids, c("m.mult", "m.uni"),
#'        main=list(label="Post. mean spatial effect",cex=2))
#' }
#' }
#'
#' @export
#'
#' @usage inla.rgeneric.IMCAR.model(cmd, theta)

# Define previous variables as global to avoid warnings()
utils::globalVariables(c("k", "W", "alpha.min", "alpha.max"))

'inla.rgeneric.IMCAR.model' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
{
  ## IMCAR implementation MCAR(1, Lambda) = MCAR(Lambda)
  ##  = Improper CAR + dense Lambda
  ## k: number of diseases/blocks
  ## W: adjacency matrix
  ## alpha.min: minimum value for alpha
  ## alpha.max: maximum value for alpha


  #theta: (k + 1) * k / 2  for lower-tri matrix by col.
  interpret.theta <- function()
  {
    #Function for changing from internal scale to external scale
    # also, build the inverse of the matrix used to model in the external scale
    # the between-disease variability and then this matrix is inverted.

    # First k parameters are the marginal precisions,
    # the rest are the correlation parameters ordered by columns.
    mprec <- sapply(theta[as.integer(1:k)], function(x) { exp(x) })
    corre <- sapply(theta[as.integer(-(1:k))], function(x) {
      (2 * exp(x))/(1 + exp(x)) - 1 })

    param <- c(mprec, corre)

    #length non-diagonal elements
    n <- (k - 1) * k / 2

    # intial matrix with 1s at the diagonal
    M <- diag(1, k)

    #Adding correlation parameters (lower.tri) and (upper.tri)
    M[lower.tri(M)] <- param[k + 1:n]
    M[upper.tri(M)] <- t(M)[upper.tri(M)]

    #Preparing the st. dev matrix
    st.dev <- 1 / sqrt(param[1:k])

    # Matrix of st. dev.
    st.dev.mat <- matrix(st.dev, ncol = 1) %*% matrix(st.dev, nrow = 1)

    # Final inversed matrix
    M <- M * st.dev.mat

    # Inverting the matrix
    #PREC <- MASS::ginv(M)# Generalized inverse
    PREC <- solve(M)

    return (list(param = param, VACOV = M, PREC = PREC))
  }


  #Graph of precision function; i.e., a 0/1 representation of precision matrix
  graph <- function()
  {

    PREC <- matrix(1, ncol = k, nrow = k)
    G <- kronecker(PREC, Matrix::Diagonal(nrow(W), 1) + W)
    return (G)
  }

  #Precision matrix
  Q <- function()
  {
    #Parameters in model scale
    param <- interpret.theta()

    #Precision matrix
    Q <- kronecker(param$PREC, Matrix::Diagonal(nrow(W), apply(W, 1, sum)) -  W)

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
    #val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation))
    #+ 0.5 * log(1.0 - param$alpha^2)

    val <- numeric(0)
    return (val)
  }

  log.prior <- function() {
    ## return the log-prior for the hyperparameters.
    param <- interpret.theta()

    # Whishart prior for joint matrix of hyperparameters
    val <- 
      log(MCMCpack::dwish(W = param$PREC, v = 2 * k + 1, S = diag(rep(1, k)))) +
      sum(theta[as.integer(1:k)]) +  # This is for precisions
      sum(log(2) + theta[-as.integer(1:k)] - 2 * log(1 + exp(theta[-as.integer(1:k)]))) # This is for correlation terms

    return (val)
  }

  initial <- function() {
    ## return initial values

    #Initial values form a diagonal matrix
    return ( c(rep(log(1), k), rep(0, (k * (k - 1) / 2))))
  }

  quit <- function() {
    return (invisible())
  }

  # FIX for rgeneric to work on R >= 4
  # Provided by E. T. Krainski
  if (as.integer(R.version$major) > 3) {
    if (!length(theta))
      theta <- initial()
  } else {
    if (is.null(theta)) {
      theta <- initial()
    }
  }

  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

##' @rdname imcar
##' @param ...  Arguments to be passed to 'inla.rgeneric.define'.
##' @export
##' @usage inla.IMCAR.model(...)
inla.IMCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.IMCAR.model, ...)
}

