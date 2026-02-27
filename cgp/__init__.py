"""
Public API for the CGP project.

Top-level users should import from:
- cgp.CGPEvolver (and any other core classes you re-export)
- cgp.run (entrypoints)
- cgp.JobScheduler (optional)
"""

from __future__ import annotations

# ---- Core public API ----
from .cgp_evolver import CartesianGP
from .cgp_model import CGP


__all__ = [
    "CartesianGP",
    "CGP",
]

