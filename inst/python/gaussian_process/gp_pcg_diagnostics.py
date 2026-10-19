"""Collect convergence diagnostics for every RHS used by a GP prediction."""
import warnings
import torch


class PCGAudit:
    def __init__(self, solver):
        self.solver = solver
        self.results = []

    def solve(self, *args, **kwargs):
        result = self.solver(*args, **kwargs)
        if not bool(torch.isfinite(result.x).all()) or not torch.isfinite(torch.tensor(result.rel_resid)):
            raise RuntimeError('GP PCG returned a non-finite solution or residual')
        # Do not retain solution vectors: covariance requests can have thousands
        # of RHSs, each as long as the training panel.
        self.results.append((bool(result.converged), float(result.rel_resid)))
        return result

    def report(self):
        failed = sum(not converged for converged, _ in self.results)
        worst = max((residual for _, residual in self.results), default=0.)
        if failed:
            warnings.warn(f'GP PCG did not converge for {failed}/{len(self.results)} right-hand sides '
                          f'(worst relative residual {worst:.2e}); prediction or uncertainty may be unreliable.',
                          RuntimeWarning, stacklevel=2)
        return dict(pcg_all_converged=(failed == 0), pcg_failed_rhs=failed,
                    pcg_rhs_count=len(self.results), pcg_worst_residual_norm=worst)
