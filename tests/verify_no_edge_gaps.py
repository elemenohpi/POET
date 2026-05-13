"""Test that no individual rule has leading/trailing gaps after evolution.

Runs a short evolutionary loop (10 generations) and inspects every rule of
every individual in the final population for leading/trailing wildcard ('_')
or variable-gap ('_{x,y}') elements. Reports violations.
"""
import os
import sys
import random as rand

# Make the project root importable when running from anywhere
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import eletility
import pop as population
import optimizer
import pattern_engine as PE


def has_edge_gap(pattern: str) -> tuple[bool, bool]:
    """Return (leading_is_gap, trailing_is_gap) for a pattern string."""
    if not pattern:
        return False, False
    elems = PE.parse_pattern(pattern)
    if not elems:
        return False, False
    return PE.is_gap_element(elems[0]), PE.is_gap_element(elems[-1])


def main(config_path: str, gens: int = 10, seed: int = 123):
    cp = eletility.ConfigParser()
    config = cp.read(config_path)
    config["runs"] = str(gens)
    config["seed"] = str(seed)

    rand.seed(seed)
    pop = population.Population(config)
    opt = optimizer.Optimizer(config, pop)
    opt.optimize()

    total_rules = 0
    violations = []
    for ind_idx, ind in enumerate(opt.P.pop):
        for rule_idx, rule in enumerate(ind.rules):
            total_rules += 1
            lead, trail = has_edge_gap(rule.pattern)
            if lead or trail:
                violations.append(
                    (ind_idx, rule_idx, rule.pattern, lead, trail)
                )

    print()
    print("=" * 60)
    print(f"Verification after {gens} generations (config: {config_path})")
    print(f"  Population size : {len(opt.P.pop)}")
    print(f"  Total rules     : {total_rules}")
    print(f"  Violations      : {len(violations)}")
    print("=" * 60)
    if violations:
        print("FAIL - rules with leading/trailing gaps:")
        for ind_idx, rule_idx, pat, lead, trail in violations[:20]:
            print(f"  ind={ind_idx} rule={rule_idx} lead={lead} trail={trail} pattern={pat!r}")
        if len(violations) > 20:
            print(f"  ... and {len(violations) - 20} more")
        sys.exit(1)
    else:
        print("PASS - no leading/trailing gaps in any rule.")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/exp_variable_gaps.ini"
    g = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    s = int(sys.argv[3]) if len(sys.argv) > 3 else 123
    main(cfg, g, s)
