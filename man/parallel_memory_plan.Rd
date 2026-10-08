\name{parallel_memory_plan}
\alias{parallel_memory_plan}
\title{Estimate PredictProR parallel worker memory}
\usage{
parallel_memory_plan(
  objects,
  n_tasks = 1L,
  workers = NULL,
  memory_budget_gb = NULL,
  payload_multiplier = NULL,
  worker_overhead_gb = NULL
)
}
\arguments{
\item{objects}{Named list of phenotype, marker, kernel, omic, cache, or
prepared-model objects used by tasks.}
\item{n_tasks}{Number of independent tasks.}
\item{workers}{Requested workers. Defaults to detected logical cores minus one,
bounded by the task count.}
\item{memory_budget_gb}{Optional worker-memory budget in GiB.}
\item{payload_multiplier}{Multiplier for transient copies and
model-specific transformations.}
\item{worker_overhead_gb}{Fixed GiB allowance for each R worker and its loaded
packages.}
}
\value{
A list containing the object inventory, serialized payload size, per-worker
estimate, memory budget, worker cap, recommended workers, and projected worker
memory.
}
\description{
Uses the same cross-platform memory assumptions as the package scheduler to
plan safe worker fan-out before a large genomic prediction run.
}
