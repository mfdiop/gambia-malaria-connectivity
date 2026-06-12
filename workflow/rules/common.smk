# Shared wildcard constraints and config helpers.
configfile: "config/config.yaml"

wildcard_constraints:
    k       = r"\d+",
    m       = r"\d+",
    subset  = r"[a-zA-Z0-9_]+",

def canonical_k():
    return config["canonical_k"]

def canonical_m():
    return config["canonical_m"]

def canonical_subset():
    return config["canonical_subset"]

def all_subsets():
    return config["subsets"]
