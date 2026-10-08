\name{parallel_shared_create}
\alias{parallel_shared_create}
\alias{parallel_shared_reference}
\alias{parallel_shared_map}
\alias{parallel_shared_status}
\title{Own and map immutable data shared across PredictProR jobs}
\usage{
parallel_shared_create(x)

parallel_shared_reference(owner)

parallel_shared_map(reference)

parallel_shared_status(x)
}
\arguments{
\item{x}{For \code{parallel_shared_create()}, an atomic vector, matrix, data
frame, or supported list that remains immutable while jobs use it. For
\code{parallel_shared_status()}, a shared owner or reference.}

\item{owner}{A shared owner returned by \code{parallel_shared_create()}, or an
existing lightweight reference.}

\item{reference}{A lightweight reference returned by
\code{parallel_shared_reference()}.}
}
\value{
\code{parallel_shared_create()} returns a \code{predictpror_shared_owner} that
contains the shared object and must remain reachable in the supervisor.
\code{parallel_shared_reference()} returns a small
\code{predictpror_shared_reference} without the data. \code{parallel_shared_map()}
returns the immutable shared object in the final worker.
\code{parallel_shared_status()} returns a one-row data frame.
}
\description{
Creates one supervisor-owned \code{mori} shared-memory region and exposes an
explicit reference-and-map handoff for jobs on the same machine.
}
\details{
Send only the lightweight reference between processes and call
\code{parallel_shared_map()} in the final process that reads the data. Keep the
owner alive until every job has completed. A reference does not extend the
lifetime of its region.

Do not send the owner or an already mapped object through a second Windows
PSOCK serialization boundary. Such forwarding can materialize the data and
defeat cross-job deduplication. Shared data are immutable and the region name
is a same-machine capability, not an authorization mechanism.
}
\examples{
owner <- parallel_shared_create(matrix(seq_len(20), nrow = 5))
reference <- parallel_shared_reference(owner)
mapped <- parallel_shared_map(reference)
stopifnot(identical(as.matrix(mapped), as.matrix(owner$data)))
parallel_shared_status(owner)
}
