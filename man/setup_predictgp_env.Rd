% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{setup_predictgp_env}
\alias{setup_predictgp_env}
\title{setup_predictgp_env}
\usage{
setup_predictgp_env(
  env_name = "predictgp",
  python_version = "3.11",
  prefer_gpu = TRUE,
  cuda = c("auto", "cpu", "cu121", "cu118"),
  torch_version = NULL,
  index_url = NULL,
  install_scikit_sparse = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{env_name}{Parameter passed to the env_name workflow option.}

\item{python_version}{Parameter passed to the python_version workflow option.}

\item{prefer_gpu}{Logical; if \code{TRUE}, install/use a CUDA-capable PyTorch build when an NVIDIA GPU is available.}

\item{cuda}{Torch CUDA wheel selector. One of \code{"auto"}, \code{"cpu"}, \code{"cu121"}, or \code{"cu118"}.}

\item{torch_version}{Optional Torch version passed through to the Python installer.}

\item{index_url}{Optional custom Python package index URL for Torch wheels.}

\item{install_scikit_sparse}{Logical; if \code{TRUE}, attempt to install
\code{scikit-sparse} for the CHOLMOD MME factorization backend. Failure is
non-fatal because the GP backend can use its SciPy sparse-LU fallback.}

\item{quiet}{Parameter passed to the quiet workflow option.}

}
\description{
Documentation for \code{setup_predictgp_env}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
