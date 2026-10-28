#!/usr/bin/env python3
"""Plot TITAN storm outlines with Matplotlib."""

import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


TARGET = (2017, 8, 12, 16, 13, 32)
INPUT_FILE = "tracks_anna"
OUTPUT_FILE = "storm_outlines_20170812_161332_matplotlib.png"


def matching_records(path):
    """Read records matching TARGET from a TITAN ASCII track file."""
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        if tuple(map(int, fields[3:9])) != TARGET:
            continue

        yield {
            "simple": int(fields[2]),
            "xc": float(fields[13]),
            "yc": float(fields[14]),
            "rays": list(map(float, fields[-72:])),
        }


def outline_points(record):
    """Convert clockwise rays from north into native TITAN X/Y points."""
    points = []
    for index, radius in enumerate(record["rays"]):
        angle = math.radians(index * 5)
        x = record["xc"] + radius * math.sin(angle)
        y = record["yc"] + radius * math.cos(angle)
        points.append((x, y))
    return points


def main():
    base = Path(__file__).resolve().parent
    records = list(matching_records(base / INPUT_FILE))
    if not records:
        raise RuntimeError(f"No {INPUT_FILE} records found for {TARGET}")

    outlines = [outline_points(record) for record in records]
    all_x = [x for outline in outlines for x, _ in outline]
    all_y = [y for outline in outlines for _, y in outline]

    # Use equal X/Y ranges so the square axes also have equal km scaling.
    padding_km = 18
    xmin, xmax = min(all_x) - padding_km, max(all_x) + padding_km
    ymin, ymax = min(all_y) - padding_km, max(all_y) + padding_km
    span = max(xmax - xmin, ymax - ymin)
    xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2

    color = "#1769aa"
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200, constrained_layout=True)

    for record, outline in zip(records, outlines):
        x, y = zip(*(outline + [outline[0]]))
        ax.fill(x, y, color=color, alpha=0.16)
        ax.plot(x, y, color=color, linewidth=0.8)
        ax.plot(record["xc"], record["yc"], ".", color=color, markersize=2.5)
        ax.annotate(
            str(record["simple"]),
            (record["xc"], record["yc"]),
            ha="center",
            va="center",
            fontsize=5,
            color=color,
            path_effects=[path_effects.withStroke(linewidth=1.5, foreground="white")],
        )

    ax.set_xlim(xmid - span / 2, xmid + span / 2)
    ax.set_ylim(ymid - span / 2, ymid + span / 2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(1)
    ax.set_xlabel("TITAN X (km)")
    ax.set_ylabel("TITAN Y (km)")
    ax.set_title("TITAN storm outlines - 2017-08-12 16:13:32")
    ax.grid(True, color="0.86", linewidth=0.5)
    ax.plot([], [], color=color, linewidth=1.2, label=INPUT_FILE)
    ax.legend(loc="upper right", fontsize=7)
    ax.text(
        0.01,
        0.01,
        f"{len(records)} records",
        transform=ax.transAxes,
        fontsize=6,
        color="0.35",
    )

    output = base / OUTPUT_FILE
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(output)
    print(f"records={len(records)}")


if __name__ == "__main__":
    main()
