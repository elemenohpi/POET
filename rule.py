# Authors: Iliya "eLeMeNOhPi" Alavy - Department of Engineering - Michigan State University
# 		   Alexander Bricco - Department of Bioengineering -  Michigan State University


class Rule:
    def __init__(self, pattern, weight, status, tree_shape=None):
        self.pattern = pattern
        self.weight = weight
        self.status = status
        self.match_direction = ""
        # Regex mode fields
        self.tree_shape = (
            tree_shape  # list representing tree structure (regex mode only)
        )
        # ── Experimental fields ──
        self.position_weights = None  # list[float] | None — per-element weights
        self.group_id = 0  # 0 = independent (no composition group)
        self.group_op = "and"  # "and" | "or"  (for composition groups)
