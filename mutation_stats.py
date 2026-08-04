# Per-operator mutation instrumentation.
#
# Counts how often each mutation operator is *attempted* (its Bernoulli roll
# succeeded) versus how often it actually *changed* the genotype. The gap
# between the two is the "wasted roll" rate that motivates the
# `mutation_selection = feasible` ablation arm.

import os


class MutationStats:
    """Tallies attempted vs. changed counts per mutation operator."""

    def __init__(self):
        self.attempted: dict[str, int] = {}
        self.changed: dict[str, int] = {}

    def record(self, name, did_change):
        self.attempted[name] = self.attempted.get(name, 0) + 1
        if did_change:
            self.changed[name] = self.changed.get(name, 0) + 1

    def realized_fraction(self):
        """Fraction of all mutation rolls that actually changed the genotype.

        This is the factor to feed into `mutation_rate_scale` for a
        rate-matched control arm: multiplying every nominal rate by it makes
        the *realized* number of applied mutations match this run.
        """
        total = sum(self.attempted.values())
        if total == 0:
            return 1.0
        return sum(self.changed.values()) / total

    # ── Reporting ───────────────────────────────────────────────────────────

    def rows(self):
        """Return [(name, attempted, changed, wasted_fraction), ...]."""
        out = []
        for name in sorted(self.attempted, key=lambda n: -self.attempted[n]):
            att = self.attempted[name]
            chg = self.changed.get(name, 0)
            out.append((name, att, chg, (att - chg) / att if att else 0.0))
        return out

    def format_table(self, label=""):
        lines = []
        lines.append("=" * 66)
        lines.append("MUTATION STATS{}".format(" - " + label if label else ""))
        lines.append("=" * 66)
        lines.append(
            "{:<24s}{:>10s}{:>10s}{:>10s}{:>10s}".format(
                "operator", "attempted", "changed", "wasted", "wasted%"
            )
        )
        lines.append("-" * 66)
        for name, att, chg, frac in self.rows():
            lines.append(
                "{:<24s}{:>10d}{:>10d}{:>10d}{:>9.1f}%".format(
                    name, att, chg, att - chg, frac * 100
                )
            )
        total_att = sum(self.attempted.values())
        total_chg = sum(self.changed.values())
        lines.append("-" * 66)
        lines.append(
            "{:<24s}{:>10d}{:>10d}{:>10d}{:>9.1f}%".format(
                "TOTAL",
                total_att,
                total_chg,
                total_att - total_chg,
                ((total_att - total_chg) / total_att * 100) if total_att else 0.0,
            )
        )
        lines.append("=" * 66)
        lines.append(
            "Rate-matched control arm: set mutation_rate_scale = {:.4f}".format(
                self.realized_fraction()
            )
        )
        lines.append("=" * 66)
        return "\n".join(lines)

    def save_csv(self, path):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(path, "w") as f:
            f.write("operator,attempted,changed,wasted,wasted_fraction\n")
            for name, att, chg, frac in self.rows():
                f.write("{},{},{},{},{:.6f}\n".format(name, att, chg, att - chg, frac))
            total_att = sum(self.attempted.values())
            total_chg = sum(self.changed.values())
            f.write(
                "TOTAL,{},{},{},{:.6f}\n".format(
                    total_att,
                    total_chg,
                    total_att - total_chg,
                    ((total_att - total_chg) / total_att) if total_att else 0.0,
                )
            )
            f.write(
                "realized_fraction,,,,{:.6f}\n".format(self.realized_fraction())
            )
