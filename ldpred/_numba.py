"""Optional Numba acceleration.

The Gibbs sampler spends essentially all of its time in a per-SNP Python loop,
which is the single largest cost in a genome-wide run. Numba compiles that loop
to native code. It is an optional dependency: without it, ``_jit`` is the
identity and the same function runs as plain Python, so results and behaviour
are unchanged, only slower.

Set ``LDPRED_DISABLE_NUMBA=1`` to force the pure-Python path even when Numba is
installed. That exists so the two paths can be compared directly, which is how
the golden outputs were checked.
"""

import os

__all__ = ['HAVE_NUMBA', 'numba_status', '_jit']

_disabled = os.environ.get('LDPRED_DISABLE_NUMBA', '').strip() not in ('', '0')

if _disabled:
    HAVE_NUMBA = False
    _reason = 'disabled via LDPRED_DISABLE_NUMBA'
else:
    try:
        from numba import njit as _njit
        HAVE_NUMBA = True
        _reason = None
    except ImportError:
        HAVE_NUMBA = False
        _reason = 'numba is not installed'


if HAVE_NUMBA:
    def _jit(func):
        # cache=True writes the compiled function to __pycache__ so the ~1s
        # compile is paid once per machine rather than once per run.
        return _njit(cache=True)(func)
else:
    def _jit(func):
        return func


def numba_status():
    """Human-readable description of whether the fast path is active."""
    if HAVE_NUMBA:
        return 'Numba acceleration enabled'
    return 'Numba acceleration unavailable (%s); using the pure Python sampler' % _reason
