.onLoad <- function(libname, pkgname) {
  reticulate::py_require("numpy")
  reticulate::py_require("scikit-learn")
  reticulate::py_require("scipy")
  reticulate::py_require("pandas")
  reticulate::py_require("rpy2")
  reticulate::py_require("matplotlib")
  reticulate::py_require("lightgbm")
  reticulate::py_require("seaborn")
  reticulate::py_require("hyperopt")
  reticulate::py_require("reservoirpy")
}
