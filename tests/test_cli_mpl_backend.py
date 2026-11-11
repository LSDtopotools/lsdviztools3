"""The lsdviz3-* commands draw with matplotlib's Agg backend whatever MPLBACKEND says.

A notebook kernel exports MPLBACKEND for its own inline plots (Colab sets
module://matplotlib_inline.backend_inline); a command run from the notebook
with "!" inherits it, and in an environment without that module matplotlib
fails on import. lsdviztools3.cli sets Agg before matplotlib is imported.
"""

from __future__ import annotations

import os
import subprocess
import sys


def test_cli_package_overrides_an_unusable_mplbackend():
    env = dict(os.environ, MPLBACKEND="module://no_such_backend_module_for_lsdviztools3")
    code = (
        "import lsdviztools3.cli\n"
        "import matplotlib\n"
        "print(matplotlib.get_backend().lower())\n"
    )
    out = subprocess.run([sys.executable, "-c", code], env=env,
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "agg"
