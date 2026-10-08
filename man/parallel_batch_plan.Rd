\name{parallel_batch_plan}
\alias{parallel_batch_plan}
\title{Plan explicit capacity for concurrent PredictProR jobs}
\usage{
parallel_batch_plan(
  objects,
  jobs,
  workers_per_job,
  tasks_per_job = workers_per_job,
  data_mode = NULL,
  memory_budget_gb = NULL,
  total_worker_limit = NULL,
  payload_multiplier = NULL,
  worker_overhead_gb = NULL
)
}
\arguments{
\item{objects}{Named list of large phenotype, genotype, kernel, omic, cache,
shared-owner, or shared-reference objects used by every job.}
\item{jobs}{Number of top-level jobs requested concurrently.}
\item{workers_per_job}{Worker processes requested by each admitted job.}
\item{tasks_per_job}{Number of tasks in each job. Worker count is capped by
this value.}
\item{data_mode}{Required explicit choice: \code{"independent"},
\code{"per_job_shared"}, or \code{"supervisor_shared"}.}
\item{memory_budget_gb}{Total batch memory budget in GiB. The default uses the
current PredictProR parallel-memory budget.}
\item{total_worker_limit}{Maximum worker slots for the entire batch. The
default is detected logical cores minus one.}
\item{payload_multiplier}{Allowance for transient payload copies and model
transformations. The default is \code{GP_PAR_PAYLOAD_MULTIPLIER} or 4.}
\item{worker_overhead_gb}{Memory allowance for each supervisor, job, or worker
R process. The default is \code{GP_PAR_WORKER_OVERHEAD_GB} or 0.25 GiB.}
}
\value{
A \code{predictpror_batch_plan} list containing an object inventory, one-row
decision table, and one row of environment settings for every admitted job.
Jobs beyond \code{jobs_admitted} must remain queued by the caller.
}
\description{
Computes a conservative and inspectable worker, memory, and job-admission plan
before concurrent PredictProR work is launched. It does not start work or
change environment variables.
}
\details{
The three data modes model different ownership boundaries. In
\code{"independent"} mode, each worker may receive an ordinary serialized copy.
In \code{"per_job_shared"} mode, every top-level job owns one shared region. In
\code{"supervisor_shared"} mode, one living supervisor owns the region and final
workers map a lightweight reference directly.

The returned \code{job_environment} table exposes the per-job value for
\code{GP_PAR_MEMORY_BUDGET_GB}, the single-thread BLAS setting, the mori switch,
and the required data handoff. The caller remains responsible for launching
only admitted jobs and keeping queued jobs inactive.
}
\examples{
plan <- parallel_batch_plan(
  objects = list(genotypes = matrix(0, nrow = 20, ncol = 20)),
  jobs = 4,
  workers_per_job = 2,
  tasks_per_job = 5,
  data_mode = "per_job_shared",
  memory_budget_gb = 4,
  total_worker_limit = 4
)
plan$decision
plan$job_environment
}
