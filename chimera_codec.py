from itertools import combinations, product, permutations


SPECIES = ("B2", "B3")
DEFAULT_REGION_COUNT = 12


def normalize_species(species):
    value = str(species).strip().upper()
    if value not in SPECIES:
        raise ValueError("Unknown species '{}'; expected B2 or B3".format(species))
    return value


def opposite_species(species):
    species = normalize_species(species)
    return "B3" if species == "B2" else "B2"


def make_token(target_slot, donor_species, donor_region):
    return "S{:02d}:{}_{:02d}".format(
        int(target_slot), normalize_species(donor_species), int(donor_region)
    )


def parse_token(token):
    slot_part, donor_part = token.split(":", 1)
    species, region = donor_part.split("_", 1)
    return {
        "target_slot": int(slot_part[1:]),
        "donor_species": normalize_species(species),
        "donor_region": int(region),
    }


def encode_construct(backbone, swaps=None, region_count=DEFAULT_REGION_COUNT):
    """Return fixed-order slot tokens for a chimera construct.

    swaps is an iterable of dictionaries or tuples. Each swap describes
    target_slot, donor_species, donor_region. Tuple order is:
    (target_slot, donor_species, donor_region).
    """
    backbone = normalize_species(backbone)
    slots = {
        slot: make_token(slot, backbone, slot)
        for slot in range(1, int(region_count) + 1)
    }

    for swap in swaps or []:
        if isinstance(swap, dict):
            target_slot = swap["target_slot"]
            donor_species = swap["donor_species"]
            donor_region = swap["donor_region"]
        else:
            target_slot, donor_species, donor_region = swap
        slots[int(target_slot)] = make_token(
            target_slot, donor_species, donor_region
        )

    return [slots[slot] for slot in range(1, int(region_count) + 1)]


def decode_construct(tokens):
    return [parse_token(token) for token in tokens]


def build_alphabet(region_count=DEFAULT_REGION_COUNT, species=SPECIES):
    tokens = []
    for target_slot in range(1, int(region_count) + 1):
        for donor_species in species:
            for donor_region in range(1, int(region_count) + 1):
                tokens.append(make_token(target_slot, donor_species, donor_region))
    return tokens


def manifest_rows(region_count=DEFAULT_REGION_COUNT, species=SPECIES):
    rows = []
    for token in build_alphabet(region_count, species):
        info = parse_token(token)
        rows.append(
            {
                "code": token,
                "target_slot": info["target_slot"],
                "donor_species": info["donor_species"],
                "donor_region": info["donor_region"],
                "is_native": info["target_slot"] == info["donor_region"],
            }
        )
    return rows


def generate_candidates(
    backbone,
    swap_count,
    region_count=DEFAULT_REGION_COUNT,
    donor_species=None,
    reuse_donors=True,
):
    """Yield fixed-order any-to-any chimera candidates.

    Each yielded item is (tokens, swaps), where swaps uses tuple order
    (target_slot, donor_species, donor_region).
    """
    backbone = normalize_species(backbone)
    donor_species = normalize_species(donor_species or opposite_species(backbone))
    slots = range(1, int(region_count) + 1)
    donor_regions = list(slots)

    for target_slots in combinations(slots, swap_count):
        if reuse_donors:
            region_iter = product(donor_regions, repeat=swap_count)
        else:
            region_iter = permutations(donor_regions, swap_count)
        for assigned_regions in region_iter:
            swaps = [
                (target_slot, donor_species, donor_region)
                for target_slot, donor_region in zip(target_slots, assigned_regions)
            ]
            yield encode_construct(backbone, swaps, region_count), swaps
