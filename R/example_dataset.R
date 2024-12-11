#' Example Dataset for RST_GAM
#'
#' An example dataset containing all the necessary components for running the RST_GAM function.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{Y}{Matrix. Simulated response variable (500 rows x 1 columns).}
#'   \item{X}{Matrix. Simulated covariates (500 rows x 3 columns).}
#'   \item{X_t}{Matrix. Time-varying covariates (500 rows x 1 columns).}
#'   \item{S}{Matrix. Location points (500 rows x 2 columns).}
#'   \item{Tr}{Matrix. Triangulation matrix for the spatial domain (109 rows x 3 columns).}
#'   \item{V}{Matrix. Vertices of the triangulation (95 rows x 2 columns).}
#' }
#'
#' @usage data(example_dataset)
#' @examples
#' # Load the example dataset
#' data(example_dataset)
#'
#' # Inspect the structure of the dataset
#' str(example_dataset)
#'
#' # Access individual components
#' Y <- example_dataset$Y
#' X <- example_dataset$X
#' X_t <- example_dataset$X_t
#' S <- example_dataset$S
#' Tr <- example_dataset$Tr
#' V <- example_dataset$V
#' 
#' # The detailed univarite and bivariate components used in this simulated data
#' # are presented in the vignette.
#' 
"example_dataset"