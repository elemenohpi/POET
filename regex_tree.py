"""Tree-based regex pattern generation for POET regex matching mode.

Provides Grow/Full/Half initialization methods that build regex patterns
as binary trees and convert them to compilable regex strings.

Original author: N. Scalzitti (11/01/2022)
Adapted for POET2.0 integration.
"""

import random
import re
import csv

import numpy as np

# ── Module-level alphabet state (initialized lazily) ────────────────────────
OPCODE: dict[str, int] = {}
OPERATOR: list[str] = []
LAST: list[str] = []
ALPHABET: list[str] = []
MAX_IN_SQUARE: int = 0

_initialized = False


def init(alphabet_file: str) -> None:
    """Load the alphabet/operator definitions from a CSV file.

    File format: each line is ``symbol;arity`` where arity is 0 (leaf),
    1 (unary operator) or 2 (binary operator).
    """
    global OPCODE, OPERATOR, LAST, ALPHABET, MAX_IN_SQUARE, _initialized

    opcode: dict[str, int] = {}
    alphabet: list[str] = []
    last: list[str] = []
    operator: list[str] = []

    with open(alphabet_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(";")
            char = parts[0]
            arity = int(parts[1])
            opcode[char] = arity

    for k, v in opcode.items():
        if v == 1 or v == 2:
            operator.append(k)
        if v == 0:
            last.append(k)
        if v == 0 and k != ".":
            alphabet.append(k)

    OPCODE = opcode
    OPERATOR = operator
    LAST = last
    ALPHABET = alphabet
    MAX_IN_SQUARE = max(1, len(LAST) // 2)
    _initialized = True


def _ensure_init() -> None:
    if not _initialized:
        raise RuntimeError(
            "regex_tree module not initialized. Call regex_tree.init(alphabet_file) first."
        )


# ── Public API: tree construction ───────────────────────────────────────────


def indi_grow(depth: int, min_braces: int, max_braces: int):
    """Create a regex via the GROW method (heterogeneous branches).

    Returns (pattern_string, tree_list) or (None, tree_list) on failure.
    """
    _ensure_init()
    g = _Grow(depth, min_braces, max_braces)
    g.build_grow_tree()
    return tree2regex(g.tree), g.tree


def indi_full(depth: int, min_braces: int, max_braces: int):
    """Create a regex via the FULL method (homogeneous branches).

    Returns (pattern_string, tree_list) or (None, tree_list) on failure.
    """
    _ensure_init()
    f = _Full(depth, min_braces, max_braces)
    f.build_full_tree()
    return tree2regex(f.tree), f.tree


def indi_half(depth: int, min_braces: int, max_braces: int):
    """Create a regex via RAMPED HALF-AND-HALF (50/50 grow or full).

    Returns (pattern_string, tree_list) or (None, tree_list) on failure.
    """
    _ensure_init()
    if np.random.choice([0, 1], p=[0.5, 0.5]) == 1:
        return indi_grow(depth, min_braces, max_braces)
    else:
        return indi_full(depth, min_braces, max_braces)


# ── Tree → regex conversion ────────────────────────────────────────────────


def _build_regex(array: list) -> str:
    """Transform the explored tree array into a readable regex string."""
    for i, node in enumerate(array):
        if isinstance(node, list):
            array[i] = "[" + "".join(node) + "]"
    return "".join(array[1:-1]).replace("cat", "")


def _explore_tree(tree: list, arr: list, index_node: int = 0) -> None:
    """In-order traversal of a tree to produce a regex array."""
    if index_node >= len(tree) or tree[index_node] is None:
        return

    if tree[index_node] in OPERATOR:
        arr.append("(")

    left = index_node * 2 + 1
    if left < len(tree):
        _explore_tree(tree, arr, left)

    node = tree[index_node]
    if isinstance(node, list):
        arr.append(node)
    elif node not in ("[]", "[^]"):
        arr.append(node)

    right = index_node * 2 + 2
    if right < len(tree):
        _explore_tree(tree, arr, right)

    if tree[index_node] in OPERATOR:
        arr.append(")")


def tree2regex(tree: list) -> str | None:
    """Convert a tree (list) to a regex string. Returns None if invalid."""
    array: list = []
    _explore_tree(tree, array, 0)
    regex = _build_regex(array)
    try:
        re.compile(regex)
    except re.error:
        return None
    return regex


# ── Full initialization ────────────────────────────────────────────────────


class _Full:
    """Build a tree with the FULL method (all branches reach max depth)."""

    def __init__(self, depth: int, min_braces: int, max_braces: int):
        self.depth = depth
        self.max_nodes = (2**depth) - 1
        self.tree: list = [0] * self.max_nodes
        self.min_braces = min_braces
        self.max_braces = max_braces
        self.leaves: list[int] = []
        self.mid_layer: list[int] = []

    def _define_leaves(self) -> None:
        if self.depth <= 1:
            self.leaves.append(0)
        else:
            for i in range(2 ** (self.depth - 1)):
                self.leaves.append((self.max_nodes - i) - 1)
            for i in range(2 ** (self.depth - 2)):
                self.mid_layer.append(
                    ((self.max_nodes - 2 ** (self.depth - 1)) - i) - 1
                )
        self.mid_layer = list(reversed(self.mid_layer))
        self.leaves = list(reversed(self.leaves))

    def _add_operator_arity1(self, i: int) -> None:
        """Handle a node in the penultimate layer (creates leaf children)."""
        self.tree[i] = random.choice(OPERATOR)

        if self.tree[i] == "[":
            self.tree[i] = "[]"
            self.tree[i * 2 + 1] = random.sample(
                ALPHABET, random.randint(1, MAX_IN_SQUARE + 1)
            )
            self.tree[i * 2 + 2] = None
        elif self.tree[i] == "^":
            self.tree[i] = "[^]"
            self.tree[i * 2 + 1] = ["^"] + random.sample(
                ALPHABET, random.randint(MAX_IN_SQUARE, len(ALPHABET) - 1)
            )
            self.tree[i * 2 + 2] = None
        elif self.tree[i] == "+":
            self.tree[i * 2 + 1] = random.choice(LAST)
            self.tree[i * 2 + 2] = None
        elif self.tree[i] == "{":
            x = random.randint(self.min_braces, self.max_braces)
            self.tree[i] = "{" + str(x) + "}"
            self.tree[i * 2 + 1] = 0
            self.tree[i * 2 + 2] = None
        elif self.tree[i] == "cat":
            self.tree[i * 2 + 1] = random.choice(LAST)
            self.tree[i * 2 + 2] = random.choice(LAST)
        elif self.tree[i] == "|":
            self.tree[i * 2 + 1] = random.choice(LAST)
            self.tree[i * 2 + 2] = random.choice(LAST)

    def build_full_tree(self) -> list:
        self._define_leaves()
        for node in range(self.max_nodes):
            if node == 0:
                self.tree[node] = random.choice(["cat", "|"])
            elif self.tree[node] == 0:
                if node in self.leaves:
                    self.tree[node] = random.choice(LAST)
                elif node in self.mid_layer:
                    self._add_operator_arity1(node)
                else:
                    self.tree[node] = random.choice(["cat", "|"])
        return self.tree


# ── Grow initialization ────────────────────────────────────────────────────


class _Grow:
    """Build a tree with the GROW method (heterogeneous branch depths)."""

    def __init__(self, depth: int, min_braces: int, max_braces: int):
        self.depth = depth
        self.max_nodes = (2**depth) - 1
        self.tree: list = [0] * self.max_nodes
        self.min_braces = min_braces
        self.max_braces = max_braces

    def _new_node(self, i: int, spe_case: bool = False) -> str:
        if i == 0:
            return random.choice(["cat", "|"])
        if spe_case:
            tmp = list(LAST) + ["cat", "^", "["]
            return random.choice(tmp)
        return random.choice(list(OPCODE.keys()))

    def _add_leaf(self, i: int) -> str:
        """Populate children for a unary operator node."""
        op = self.tree[i]
        if op == "[":
            a = random.randint(1, MAX_IN_SQUARE)
            self.tree[i * 2 + 1] = random.sample(ALPHABET, a)
            return "[]"
        elif op == "^":
            self.tree[i * 2 + 1] = ["^"] + random.sample(
                ALPHABET, random.randint(MAX_IN_SQUARE, len(ALPHABET) - 1)
            )
            return "[^]"
        elif op == "{":
            x = random.randint(self.min_braces, self.max_braces)
            self.tree[i * 2 + 1] = 1
            return "{" + str(x) + "}"
        elif op == "+":
            self.tree[i * 2 + 1] = 1
            return "+"
        return op

    def _add_child(self, i: int, arity: int) -> None:
        if arity == 0:
            try:
                self.tree[i * 2 + 1] = None
                self.tree[i * 2 + 2] = None
            except IndexError:
                pass
        elif arity == 1:
            try:
                self.tree[i] = self._add_leaf(i)
                self.tree[i * 2 + 2] = None
            except IndexError:
                pass
        elif arity == 2:
            try:
                self.tree[i * 2 + 1] = 0
                self.tree[i * 2 + 2] = 0
            except IndexError:
                pass

    def _check_last_layer(self) -> None:
        for i in range(2 ** (self.depth - 1)):
            idx = -(i + 1)
            if isinstance(self.tree[idx], list):
                pass
            elif self.tree[idx] not in LAST and self.tree[idx] is not None:
                self.tree[idx] = random.choice(LAST)

    def build_grow_tree(self) -> list:
        for i in range(self.max_nodes):
            if self.tree[i] == 0:
                self.tree[i] = self._new_node(i)
                self._add_child(i, OPCODE[self.tree[i]])
            elif self.tree[i] == 1:
                self.tree[i] = self._new_node(i, spe_case=True)
                self._add_child(i, OPCODE[self.tree[i]])
            elif self.tree[i] is None:
                try:
                    self.tree[i * 2 + 1] = None
                    self.tree[i * 2 + 2] = None
                except IndexError:
                    pass
            else:
                self._add_child(i, 0)
        self._check_last_layer()
        return self.tree
