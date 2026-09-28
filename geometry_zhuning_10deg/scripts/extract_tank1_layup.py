# -*- coding: utf-8 -*-
"""Extract Tank-1 hoop-section layer radii from a job INP (WCM continuum model).

For each given INP file, scan the *Part, name=Tank-1 block, keep nodes with
|y| < 5 mm (cylinder body) and cluster their radial distance r = sqrt(x^2+z^2)
into layer groups. Prints inner/outer radius, total thickness, group count and
per-group radii so the WCM layup can be documented/compared.

Usage:
    python extract_tank1_layup.py [inp1.inp inp2.inp ...]
Default (no args): the three successful jobs in the parent directory.
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def scan_tank1_radii(inp_path):
    radii = []
    in_part = False
    in_nodes = False
    with open(inp_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("*Part, name=Tank-1"):
                in_part = True
                continue
            if in_part and line.startswith("*End Part"):
                break
            if not in_part:
                continue
            if line.startswith("*Node"):
                in_nodes = True
                continue
            if in_nodes:
                if line.startswith("*"):
                    in_nodes = False
                    continue
                parts = line.split(",")
                if len(parts) < 4:
                    continue
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError:
                    continue
                if abs(y) < 5.0:
                    radii.append(math.hypot(x, z))
    return radii


def cluster(sorted_r, gap=0.35):
    groups = []
    cur = [sorted_r[0]]
    for r in sorted_r[1:]:
        if r - cur[-1] > gap:
            groups.append(cur)
            cur = []
        cur.append(r)
    groups.append(cur)
    return groups


def report(inp_path):
    radii = scan_tank1_radii(inp_path)
    if not radii:
        print(f"{inp_path.name}: no Tank-1 body nodes found")
        return
    uniq = sorted(set(round(r, 4) for r in radii))
    groups = cluster(uniq)
    means = [sum(g) / len(g) for g in groups]
    print(f"{inp_path.name}: body_nodes={len(radii)} distinct_r={len(uniq)} "
          f"r_min={uniq[0]:.3f} r_max={uniq[-1]:.3f} "
          f"total_thickness={uniq[-1] - uniq[0]:.3f} groups={len(groups)}")
    print("  group radii (mm): " + ", ".join(f"{m:.3f}" for m in means))
    th = [means[i + 1] - means[i] for i in range(len(means) - 1)]
    print("  group thickness (mm): " + ", ".join(f"{t:.3f}" for t in th))


def main():
    args = sys.argv[1:]
    files = [Path(a) if Path(a).is_absolute() else ROOT / a for a in args] if args else [
        ROOT / "Job-2.inp",
        ROOT / "Job-multi.inp",
        ROOT / "WCM_multi.inp",
    ]
    for f in files:
        report(f)


if __name__ == "__main__":
    main()
