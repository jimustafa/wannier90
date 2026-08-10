#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "matplotlib",
# ]
# ///
"""Plot wannier90 OPFM tutorial output files."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_wmat(path):
    with open(path) as f:
        header = f.readline().strip()
        num_proj, num_wann = (int(x) for x in f.readline().split())
        vals = np.loadtxt(f)

    w = vals[:, 0] + 1j * vals[:, 1]
    w = w.reshape((num_proj, num_wann), order="F")
    return w.T, header


def plot_w_matrix(args):
    """Plot a <seedname>_opfm_w.mat file (the M x N codiagonalization matrix W)."""
    w, header = read_wmat(args.wmat_file)

    num_wann, num_proj = w.shape

    fig, ax = plt.subplots(figsize=(6, 4.5), constrained_layout=True)

    im = ax.imshow(np.abs(w), cmap="viridis", vmin=0, vmax=1, aspect="equal")
    ax.set_title(r"$|[W^T]_{ij}|$")
    ax.set_xlabel("Projection Index")
    ax.set_ylabel("OPF Index")
    ax.set_xticks(range(num_proj), labels=range(1, num_proj + 1))
    ax.set_yticks(range(num_wann), labels=range(1, num_wann + 1))

    fig.colorbar(im, ax=ax, orientation="horizontal", location="top")

    fig.suptitle(f"{args.wmat_file}  ({header})")

    fig.savefig(args.output, bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="kind", required=True)

    w_matrix = subparsers.add_parser("w-matrix", help=plot_w_matrix.__doc__)
    w_matrix.add_argument("wmat_file", help="<seedname>_opfm_w.mat")
    w_matrix.add_argument(
        "-o", "--output", help="output image file (default: <wmat_file>.png)"
    )
    w_matrix.set_defaults(func=plot_w_matrix)

    args = parser.parse_args()
    if args.output is None:
        args.output = Path(args.wmat_file).with_suffix(".png")
    args.func(args)


if __name__ == "__main__":
    main()
