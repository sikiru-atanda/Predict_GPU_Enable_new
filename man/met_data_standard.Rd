\name{met_data_standard}
\alias{met_data_standard}
\title{Describe the single-trait multi-environment input standard}
\usage{
met_data_standard() 
}
\description{
Returns the phenotype, genomic, cross-validation, supported-model, and
model-purpose contract for Gaussian, binary, multiclass, and ordinal
single-trait multi-environment analyses.
The returned scope distinguishes structured genetic/GxE models from
environment-aware predictive ML/DL models and lists models that are rejected
instead of being forced through an incompatible route.
}
\value{
A list describing the current single-trait MET input and model scope.
}
