"""lsdviztools3 — PyGMT-based visualisation tools for LSDTopoTools v3 output.

Version 1 provides a backend-agnostic IO/data core (``lsdviztools3.io``) and a
thin PyGMT render layer (``lsdviztools3.render``). Statistics (x-y) plotting is
deferred (see ``lsdviztools3.stats``).
"""

__version__ = "0.1.0"

# TODO: curate top-level re-exports once the render API stabilises, e.g.
#   from .render.hillshade import render_hillshade
#   from .render.style import MapStyle
