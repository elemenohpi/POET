import random as R
from functools import lru_cache

import pandas as pd


WILDCARD = "_"


def is_token_mode(config):
    return config.get("sequence_mode", "char").strip().lower() == "token"


def token_delimiter(config):
    raw = config.get("token_delimiter", "space").strip()
    if raw.lower() in ("space", "whitespace", ""):
        return None
    if raw.lower() == "tab":
        return "\t"
    return raw


def split_tokens(text, config):
    if not isinstance(text, str) or text == "":
        return []
    delimiter = token_delimiter(config)
    if delimiter is None:
        return text.split()
    return [part for part in text.split(delimiter) if part != ""]


def join_tokens(tokens, config):
    delimiter = token_delimiter(config)
    sep = " " if delimiter is None else delimiter
    return sep.join(tokens)


def pattern_size(pattern, config):
    if not pattern:
        return 0
    if is_token_mode(config):
        return len(split_tokens(pattern, config))
    return len(pattern)


@lru_cache(maxsize=32)
def _load_alphabet_cached(path, column):
    df = pd.read_csv(path)
    if column not in df.columns:
        raise ValueError(
            "Alphabet file '{}' does not contain column '{}'".format(path, column)
        )
    values = [str(v) for v in df[column].dropna().tolist()]
    if not values:
        raise ValueError("Alphabet file '{}' produced no symbols".format(path))
    return tuple(values)


def load_alphabet(config):
    path = config.get("alphabet_data", "data/translation/amino_to_amino.csv")
    column = config.get("alphabet_column", "code")
    return list(_load_alphabet_cached(path, column))


def _training_token_sequences(config):
    path = config.get("learn_data")
    if not path:
        return []
    df = pd.read_csv(path)
    seqs = []
    for value in df.iloc[:, 0].dropna().tolist():
        tokens = split_tokens(str(value), config)
        if tokens:
            seqs.append(tokens)
    return seqs


def random_observed_token_pattern(config, max_size):
    seqs = _training_token_sequences(config)
    if not seqs:
        alphabet = load_alphabet(config)
        return [R.choice(alphabet)]

    seq = R.choice(seqs)
    size = R.randint(1, min(max_size, len(seq)))
    start = R.randint(0, len(seq) - size)
    return seq[start : start + size]


def token_slot(token):
    if ":" in token:
        return token.split(":", 1)[0]
    if "=" in token:
        return token.split("=", 1)[0]
    return ""


def token_replacement_pool(token, alphabet):
    slot = token_slot(token)
    if not slot:
        return [item for item in alphabet if item != token]
    pool = [item for item in alphabet if item != token and token_slot(item) == slot]
    return pool or [item for item in alphabet if item != token]


def strip_edge_wildcards(tokens):
    start = 0
    end = len(tokens)
    while start < end and tokens[start] == WILDCARD:
        start += 1
    while end > start and tokens[end - 1] == WILDCARD:
        end -= 1
    return tokens[start:end]


def token_pattern_is_all_wildcard(tokens):
    return bool(tokens) and all(token == WILDCARD for token in tokens)


def token_pattern_matches(pattern_tokens, sequence_tokens, pos, gaps_enabled=True):
    if pos + len(pattern_tokens) > len(sequence_tokens):
        return False
    for offset, token in enumerate(pattern_tokens):
        if gaps_enabled and token == WILDCARD:
            continue
        if token != sequence_tokens[pos + offset]:
            return False
    return True
