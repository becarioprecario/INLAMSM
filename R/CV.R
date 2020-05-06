#' @name CV
#' @rdname CV
#' @docType data
#' @keywords datasets
#' @title Multivariate mortality data from Comunidad Valenciana (Spain)
#' 
#' @format A \code{SpatialPolygonsDataFrame} with the boundaries of the 
#'  municipalities in Comunidad Valenciana with the following columns:
#' \describe{
#'   \item{CODMUNI}{Municipality code.}
#'   \item{NOMBRE}{Name of the municipality.}
#'   \item{Exp.Cirrhosis}{Expected number of cases of cirrhosis.}
#'   \item{Exp.Lung}{Expected number of cases of lung cancer.}
#'   \item{Exp.Oral}{Expected number of cases of oral cavity cancer.}
#'   \item{Obs.Cirrhosis}{Observed number of cases of cirrhosis.}
#'   \item{Obs.Lung}{Observed number of cases of lung cancer.}
#'   \item{Obs.Oral}{Observed number of cases of oral cavity cancer.}
#' }
#'
#' @usage data(CV)
#' @source The original data set is supplied as supplementary material of the
#' book: "Martinez-Beneito, M A & Botella Rocamora, P. Disease mapping: from
#' foundations to multidimensional modeling. CRC/Chapman & Hall, 2019". This
#' object has been built from several of the files available at the 
#' supplementary material repository of the book at:
#' \url{https://github.com/MigueBeneito/DisMapBook/tree/master/Data} 
#' 
#' @description Simulated multivariate mortality data from Comunidad Valenciana
#' (Spain). The data set contains (simulated) observed and expected deaths for
#' Cirrhosis, Lung cancer and Cirrhosis for the Valencian municipalities. The
#' supplied data have been simulated mimicking the original data set which has
#' privacy restrictions. Additional details on the generation of the supplied 
#' dataset can be found at the original book.
#'
#' @references Martinez-Beneito, M A & Botella Rocamora, P. Disease mapping: 
#' from foundations to multidimensional modeling. CRC/Chapman & Hall, 2019.
#'
#' @seealso CV.nb
#' 
#' @examples
#'
#' \donttest{
#' require(INLA)
#' require(spdep)
#' data(CV)
#' W <- as(nb2mat(CV.nb, style = "B"), "Matrix")
#'
#' #Data (two diseases only)
#' d <- list(OBS = c(CV$Obs.Cirrhosis, CV$Obs.Lung),
#'  EXP = c(CV$Exp.Cirrhosis, CV$Exp.Lung))
#'
#'  # Index for latent effect
#' d$idx <- 1:length(d$OBS)
#'
#' k <- 2  #Number of diseases
#'
#' # Three independent ICAR models
#' model <- inla.rgeneric.define(inla.rgeneric.indep.IMCAR.model,
#'   k = k, W = W)
#' r.simcar <- try(
#'   inla(OBS ~ 1 + f(idx, model = model),
#'     data = d, E = EXP, family = "poisson",
#'      # To run faster, REMOVE in real applications
#'      control.mode = list(theta = c(1.4, 2.1), restart = TRUE),
#'     control.predictor = list(compute = TRUE))
#' )
#' summary(r.simcar)
#'
#' # One IMCAR model
#' model <- inla.rgeneric.define(inla.rgeneric.IMCAR.model,
#'   k = k, W = W, alpha.min = 0, alpha.max = 1)
#' r.imcar <- try(
#'   inla(OBS ~ 1 + f(idx, model = model),
#'     data = d, E = EXP, family = "poisson",
#'      # To run faster, REMOVE in real applications
#'      control.mode = list(theta = c(1.77, 2.01, 0.93),
#'        restart = TRUE),
#'     control.compute = list(config = TRUE),
#'     control.predictor = list(compute = TRUE))
#' )
#' summary(r.imcar)
#'
#' # Transform parameters
#' summary.post <- inla.MCAR.transform(r.imcar, k = k)
#'
#' # Posterior of variance matrix
#' summary.post$VAR.p # Using point estimates
#' summary.post$VAR.m # Using posterior sampling
#'
#'}
NULL

#' @rdname CV_nb
#' @name CV.nb
#' @title Adjacency matrix of municipalities in Comunidad Valenciana (Spain)
#' @format An \code{nb} object with the adjacencies of the municipalities in
#' Comunidad Valenciana (Spain) to be used for the spatial models.
#'
#' @docType data
#' @keywords datasets
#' @usage data(CV)
#' @source The original data set is supplied as supplementary material of the
#' book: "Martinez-Beneito, M A & Botella Rocamora, P. Disease mapping: from
#' foundations to multidimensional modeling. CRC/Chapman & Hall, 2019". This
#' object has been built from several of the files available at the 
#' supplementary material repository of the book at:
#' \url{https://github.com/MigueBeneito/DisMapBook/tree/master/Data} 
#' 
#' @references Martinez-Beneito, M A & Botella Rocamora, P. Disease mapping: 
#' from foundations to multidimensional modeling. CRC/Chapman & Hall, 2019.
#' 
#' @seealso CV
#'
#' @examples
#' require(spdep)
#' data(CV)
#' plot(CV)
#' plot(CV.nb, coordinates(CV), pch = ".", col = "gray", add = TRUE)
#'
NULL

