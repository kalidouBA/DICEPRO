.onLoad <- function(libname, pkgname) {
  reticulate::py_require("numpy")
  reticulate::py_require("scikit-learn")
  reticulate::py_require("scipy")
  reticulate::py_require("pandas")
  reticulate::py_require("deflate")
  reticulate::py_require("rpy2")
  reticulate::py_require("matplotlib")
  reticulate::py_require("lightgbm")
  reticulate::py_require("seaborn")
  reticulate::py_require("hyperopt")
  reticulate::py_require("reservoirpy")

  # rpy2_import <- try(reticulate::import("rpy2"), silent = TRUE)
  # if(inherits(rpy2_import, "try-error")){
  #   stop("Python module `rpy2` is missing. ",
  #        "Try installing it by running the following R command in the console:\n",
  #         "reticulate::py_install('rpy2')")
  # }

}
