"""Deterministic CD-HIT clustering and capacity-constrained high/low quotas."""
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from .calibration import fingerprint, subseed
from .io import write_fasta

def cluster_sequences(sequences, directory, identity=.5, coverage=.8):
    if len(set(sequences)) != len(sequences):
        raise ValueError("Clustering input must be unique")
    if not .4 <= identity <= 1 or not 0 < coverage <= 1:
        raise ValueError("Unsupported CD-HIT thresholds")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    ordered = sorted(sequences, key=lambda s: (-len(s), fingerprint(s)))
    source, target = directory / "input.fasta", directory / "clusters.fasta"
    write_fasta(source, ordered)
    word = 2 if identity < .5 else 3 if identity < .6 else 4 if identity < .7 else 5
    command = ["cd-hit", "-i", str(source), "-o", str(target), "-c", str(identity),
               "-n", str(word), "-aS", str(coverage), "-aL", str(coverage),
               "-G", "1", "-g", "1", "-l", "7", "-d", "0", "-T", "1", "-M", "0"]
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    groups, members = [], []
    for line in Path(str(target) + ".clstr").read_text().splitlines():
        if line.startswith(">"):
            if members:
                groups.append(members)
            members = []
        else:
            member = line.split(">seq_", 1)[1].split("...", 1)[0]
            members.append(ordered[int(member)])
    if members:
        groups.append(members)
    result = {}
    for group in groups:
        cluster = fingerprint(sorted(group))
        for sequence in group:
            result[sequence] = cluster
    if set(result) != set(sequences):
        raise ValueError("CD-HIT lost sequences, including possibly short peptides")
    return result, command

def rank_key(row, low=False):
    score = row["ranking"]
    return (score if low else -score, -row["weakest"], -row.get("motif_support", 0),
            fingerprint(row["sequence"]))

def assign_tiers(rows, seed):
    groups = defaultdict(list)
    for row in rows:
        if not row["passed"] or not math.isfinite(row["ranking"]) or not math.isfinite(row["weakest"]):
            raise ValueError("Only completely qualified finite-score candidates may enter selection")
        groups[row["cluster"]].append(dict(row))
    result = []
    for cluster in sorted(groups):
        group = sorted(groups[cluster], key=rank_key)
        if len(group) == 1:
            group[0]["tier"] = "high" if subseed(seed, "singleton", cluster) % 2 == 0 else "low"
        else:
            for index, row in enumerate(group):
                row["tier"] = "high" if index < (len(group) + 1) // 2 else "low"
        result.extend(group)
    return result

def allocate(capacities, total):
    """Iterative capped Hamilton apportionment with fixed sqrt(capacity) weights."""
    if total < 0 or any(c < 0 for c in capacities.values()) or total > sum(capacities.values()):
        raise ValueError("Insufficient tier capacity")
    result = {key: 0 for key in sorted(capacities)}
    active = {key for key, cap in capacities.items() if cap}
    remaining = total
    while remaining:
        weights = {key: math.sqrt(capacities[key]) for key in active}
        denominator = sum(weights.values())
        shares = {key: remaining * weights[key] / denominator for key in active}
        capped = [key for key in active if shares[key] >= capacities[key] - result[key]]
        if capped:
            for key in sorted(capped):
                amount = capacities[key] - result[key]
                result[key] += amount
                remaining -= amount
                active.remove(key)
            continue
        floors = {key: math.floor(shares[key]) for key in active}
        for key in active:
            result[key] += floors[key]
        remaining -= sum(floors.values())
        order = sorted(active, key=lambda key: (-(shares[key] - floors[key]), key))
        for key in order[:remaining]:
            result[key] += 1
        remaining = 0
    return result

def select(rows, count, low_fraction=.5):
    if not 0 <= low_fraction <= 1 or count < 1:
        raise ValueError("Invalid quota")
    if len({r["sequence"] for r in rows}) != len(rows):
        raise ValueError("Duplicate selection input")
    low = math.floor(count * low_fraction)
    selected = []
    for tier, needed in (("high", count - low), ("low", low)):
        groups = defaultdict(list)
        for row in rows:
            if row["tier"] == tier:
                groups[row["cluster"]].append(row)
        quotas = allocate({key: len(group) for key, group in groups.items()}, needed)
        for key in sorted(groups):
            selected.extend(sorted(groups[key], key=lambda r: rank_key(r, tier == "low"))[:quotas[key]])
    return sorted(selected, key=rank_key)
