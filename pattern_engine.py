"""Experimental pattern matching engine for POET.

Parses enriched pattern strings into element tuples and provides matching
functions that support:

- Concrete characters:     ``A``
- Single-char wildcards:   ``_``
- Character classes:       ``[AG]``  (match A or G)
- Variable-length gaps:    ``_{2,5}`` (match 2-5 of any character)

Element types (tuples)
----------------------
('c', ch)            – match exact character *ch*
('w',)               – match any single character
('cc', [ch, ...])    – match any character in the class
('vg', min, max)     – match *min* to *max* of any character
"""

import random as R

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_pattern(pattern_str):
    """Parse a pattern string into a list of element tuples."""
    elements = []
    i = 0
    n = len(pattern_str)
    while i < n:
        ch = pattern_str[i]
        if ch == "[":
            # Character class [AGK]
            try:
                end = pattern_str.index("]", i + 1)
            except ValueError:
                # Malformed: '[' without ']' – treat remaining as concrete
                for c in pattern_str[i:]:
                    elements.append(("c", c))
                break
            chars = sorted(set(pattern_str[i + 1 : end]))
            if len(chars) <= 1:
                # Degenerate class – convert to concrete or skip
                if chars:
                    elements.append(("c", chars[0]))
                # else empty brackets – skip
            else:
                elements.append(("cc", chars))
            i = end + 1
        elif ch == "_":
            # Check for variable gap _{min,max}
            if i + 1 < n and pattern_str[i + 1] == "{":
                try:
                    brace_end = pattern_str.index("}", i + 2)
                except ValueError:
                    # Malformed _{... without } – treat as plain wildcard
                    elements.append(("w",))
                    i += 1
                    continue
                parts = pattern_str[i + 2 : brace_end].split(",")
                min_g = int(parts[0])
                max_g = int(parts[1]) if len(parts) > 1 else min_g
                elements.append(("vg", min_g, max_g))
                i = brace_end + 1
            else:
                elements.append(("w",))
                i += 1
        else:
            elements.append(("c", ch))
            i += 1
    return elements


def render_elements(elements):
    """Reconstruct the pattern string from an element list."""
    parts = []
    for elem in elements:
        if elem[0] == "c":
            parts.append(elem[1])
        elif elem[0] == "w":
            parts.append("_")
        elif elem[0] == "cc":
            parts.append("[" + "".join(elem[1]) + "]")
        elif elem[0] == "vg":
            if elem[1] == elem[2]:
                parts.append("_{" + str(elem[1]) + "}")
            else:
                parts.append("_{" + str(elem[1]) + "," + str(elem[2]) + "}")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Element queries
# ---------------------------------------------------------------------------


def reverse_elements(elements):
    """Reverse the element list for reverse-direction matching."""
    return list(reversed(elements))


def elements_min_len(elements):
    """Minimum number of sequence characters this pattern can span."""
    total = 0
    for e in elements:
        if e[0] in ("c", "w", "cc"):
            total += 1
        elif e[0] == "vg":
            total += e[1]
    return total


def elements_max_len(elements):
    """Maximum number of sequence characters this pattern can span."""
    total = 0
    for e in elements:
        if e[0] in ("c", "w", "cc"):
            total += 1
        elif e[0] == "vg":
            total += e[2]
    return total


def has_variable_length(elements):
    """True if any element is a variable-length gap."""
    return any(e[0] == "vg" for e in elements)


def is_all_wildcard(elements):
    """True if every element is a wildcard or variable gap (no specificity)."""
    return all(e[0] in ("w", "vg") for e in elements)


def concrete_count(elements):
    """Number of concrete (non-wildcard) elements."""
    return sum(1 for e in elements if e[0] in ("c", "cc"))


# ---------------------------------------------------------------------------
# Matching – exact (boolean)
# ---------------------------------------------------------------------------


def match_at(elements, sequence, pos):
    """Try to match *elements* starting at *pos*.  Returns True on match."""
    return _match_from(elements, 0, sequence, pos, len(sequence))


def _match_from(elements, ei, sequence, pos, seq_len):
    """Recursive matcher with backtracking for variable gaps."""
    if ei == len(elements):
        return True

    elem = elements[ei]

    if elem[0] == "c":
        if pos < seq_len and sequence[pos] == elem[1]:
            return _match_from(elements, ei + 1, sequence, pos + 1, seq_len)
        return False

    if elem[0] == "w":
        if pos < seq_len:
            return _match_from(elements, ei + 1, sequence, pos + 1, seq_len)
        return False

    if elem[0] == "cc":
        if pos < seq_len and sequence[pos] in elem[1]:
            return _match_from(elements, ei + 1, sequence, pos + 1, seq_len)
        return False

    if elem[0] == "vg":
        min_g, max_g = elem[1], elem[2]
        for gap_len in range(min_g, max_g + 1):
            if pos + gap_len <= seq_len:
                if _match_from(elements, ei + 1, sequence, pos + gap_len, seq_len):
                    return True
        return False

    return False


# ---------------------------------------------------------------------------
# Matching – weighted / partial
# ---------------------------------------------------------------------------


def match_at_weighted(elements, sequence, pos, position_weights, threshold=0.7):
    """Weighted partial matching for **fixed-length** patterns only.

    Returns match quality (0.0–1.0) if >= *threshold*, else 0.0.
    Variable-gap patterns are not supported — falls back to exact matching.
    """
    n = len(elements)
    if n == 0:
        return 0.0
    if len(position_weights) != n:
        # Mismatched lengths – fall back to exact match
        return 1.0 if match_at(elements, sequence, pos) else 0.0
    if pos + n > len(sequence):
        return 0.0

    total_weight = sum(position_weights)
    if total_weight == 0:
        return 0.0

    matched_weight = 0.0
    for i, elem in enumerate(elements):
        ch = sequence[pos + i]
        matched = False
        if elem[0] == "c":
            matched = ch == elem[1]
        elif elem[0] == "w":
            matched = True
        elif elem[0] == "cc":
            matched = ch in elem[1]
        elif elem[0] == "vg":
            # Variable gaps break fixed-length assumption; fall back
            return 1.0 if match_at(elements, sequence, pos) else 0.0

        if matched:
            matched_weight += position_weights[i]

    quality = matched_weight / total_weight
    return quality if quality >= threshold else 0.0


# ---------------------------------------------------------------------------
# Random element generation helpers (used by init / mutations)
# ---------------------------------------------------------------------------


def random_element(codes, gap_chance=0.10, class_chance=0.15, max_class_size=None):
    """Generate one random element.

    *codes* is the list of amino-acid single-letter codes.
    *max_class_size* optionally caps the number of amino acids in a
    generated character class. ``None`` or values <= 0 mean unlimited
    (default upper bound of 4 still applies for initial classes).
    If ``max_class_size`` is 1, character classes are disabled and a
    concrete element is returned instead.
    """
    r = R.random()
    if r < gap_chance:
        return ("w",)
    elif r < gap_chance + class_chance:
        upper = min(4, len(codes))
        if max_class_size is not None and max_class_size > 0:
            upper = min(upper, max_class_size)
        if upper < 2:
            return ("c", R.choice(codes))
        size = R.randint(2, upper)
        chars = sorted(set(R.sample(codes, size)))
        return ("cc", chars)
    else:
        return ("c", R.choice(codes))


def random_variable_gap(max_gap=5):
    """Generate a random variable-length gap element."""
    lo = R.randint(1, max(1, max_gap - 1))
    hi = R.randint(lo, max_gap)
    return ("vg", lo, hi)


def is_gap_element(elem):
    """True if *elem* is a wildcard ('w') or variable-length gap ('vg')."""
    return elem[0] in ("w", "vg")


def strip_edge_gaps(elements):
    """Remove leading and trailing gap/wildcard elements.

    Returns a new list with all leading and trailing ``('w',)`` and
    ``('vg', lo, hi)`` elements stripped. Gaps are only allowed in the
    interior of a pattern.
    """
    start = 0
    end = len(elements)
    while start < end and is_gap_element(elements[start]):
        start += 1
    while end > start and is_gap_element(elements[end - 1]):
        end -= 1
    return elements[start:end]
