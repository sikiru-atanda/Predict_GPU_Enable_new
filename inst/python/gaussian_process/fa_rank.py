"""Necessary parameter-count bound for FA loadings plus specific variances.

Equality is saturated, not a guarantee of identification for particular data.
"""
import math


def max_fa_rank(groups):
    n = int(groups)
    if n < 1:
        return 0
    return max(0, math.floor((2*n+1-math.sqrt(8*n+1))/2+1e-9))


def validate_fa_rank(rank, groups, name='fa_rank'):
    try:
        k = int(rank)
        valid = k >= 1 and k == float(rank)
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        raise ValueError(f'{name} must be a positive integer')
    n = int(groups)
    maximum = max_fa_rank(n)
    if k > maximum:
        count = n*(k+1)-k*(k-1)//2
        entries = n*(n+1)//2
        alternative = (f'Use rank {maximum} or less, or unstructured covariance.' if maximum else
                       'Use unstructured covariance where supported, or a non-FA GP model.')
        raise ValueError(f'{name}={k} over {n} groups has {count} parameters for '
                         f'{entries} covariance entries and cannot be identified. {alternative}')
    return k
