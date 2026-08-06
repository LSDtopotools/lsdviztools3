"""GitHub-compatible heading slugify for MkDocs.

GitHub strips unicode sub/superscripts to their ASCII digits when building
heading anchors (e.g. ``log\u2081\u2080`` -> ``log10``); pymdownx's default
slugifier keeps the raw unicode.  This mirrors the GitHub behaviour so that
anchors written for GitHub (\`\`#42-channel-network-coloured-by-log10-k_sn\`\`)
resolve on the rendered docs site too.
"""

from __future__ import annotations

import re
import unicodedata

from pymdownx.slugs import slugify

_TRANSLATE = {
    "\u2080": "0", "\u2081": "1", "\u2082": "2", "\u2083": "3", "\u2084": "4",
    "\u2085": "5", "\u2086": "6", "\u2087": "7", "\u2088": "8", "\u2089": "9",
    "\u00b2": "2", "\u00b3": "3", "\u00b9": "1", "\u00b0": "",
}

_SUBSCRIPTS = re.compile("|".join(re.escape(c) for c in _TRANSLATE))


def github_slugify(text, sep):
    text = _SUBSCRIPTS.sub(lambda m: _TRANSLATE[m.group(0)], unicodedata.normalize("NFKC", text))
    return slugify(case="lower")(text, sep)
