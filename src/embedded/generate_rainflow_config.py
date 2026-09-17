#!/usr/bin/env python3
"""
generate_rainflow_config.py

Writes rainflow_config.c with the damage lookup table RF_DAMAGE_LUT.
Values are computed as:

    RF_DAMAGE_LUT[i] = i ** k   for i = 0 .. RF_NUM_CLASSES-1

i.e. numpy.arange(RF_NUM_CLASSES) ** k
(RF_DAMAGE_LUT[0] is therefore always 0)

Usage:
    python generate_rainflow_config.py --num-classes 128 --k 5 \
        --output rainflow_config.c
"""

import argparse

import numpy as np

C_HEADER = '#include "rainflow.h"\n\n'

VALUES_PER_LINE = 8
# Insert a blank line after this many values (purely visual, as in the
# example with two blocks of 64). 0/None = no blank lines.
BLOCK_SIZE = 64


def compute_lut(num_classes: int, k: float) -> np.ndarray:
    """Compute LUT values: arange(N)**k, so RF_DAMAGE_LUT[0] == 0."""
    values = np.arange(num_classes, dtype=np.float64) ** k
    values = np.rint(values).astype(np.uint64)
    return values


def format_lut(values: np.ndarray, block_size: int = BLOCK_SIZE) -> str:
    lines = []
    row = []
    for idx, v in enumerate(values, start=1):
        row.append(f"{v}")
        if len(row) == VALUES_PER_LINE:
            lines.append("    " + ", ".join(row) + ",")
            row = []
        if block_size and idx % block_size == 0 and idx != len(values):
            lines.append("")
    if row:
        lines.append("    " + ", ".join(row) + ",")
    return "\n".join(lines)


def generate_c_file(num_classes: int, k: float, output_path: str) -> None:
    values = compute_lut(num_classes, k)
    body = format_lut(values)

    content = (
        C_HEADER
        + "const uint32_t RF_DAMAGE_LUT[RF_NUM_CLASSES] = {\n"
        + body
        + "\n};\n"
    )

    with open(output_path, "w") as f:
        f.write(content)

    print(f"Wrote: {output_path} ({num_classes} values, k={k})")


def main():
    parser = argparse.ArgumentParser(
        description="Write rainflow_config.c with RF_DAMAGE_LUT = arange(N)**k",
    )
    parser.add_argument(
        "--num-classes", "-n", type=int, default=128,
        help="Number of LUT classes (RF_NUM_CLASSES), default: 128",
    )
    parser.add_argument(
        "--k", "-k", type=float, default=5.0,
        help="Exponent k for the damage calculation, default: 5.0",
    )
    parser.add_argument(
        "--output", "-o", type=str, default="rainflow_config.c",
        help="Output file, default: rainflow_config.c",
    )
    args = parser.parse_args()

    generate_c_file(args.num_classes, args.k, args.output)


if __name__ == "__main__":
    main()
