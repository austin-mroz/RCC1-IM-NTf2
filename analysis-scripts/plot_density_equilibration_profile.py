import argparse
import pathlib
import warnings
from dataclasses import dataclass

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv("10_ns_npt_equilibration_density_profile.csv")
    df = df[:-1]
    print(df.loc[:, "Density"].mean())
    print(df.loc[:, "Density"].std())

    fig, ax = plt.subplots(figsize=(4, 4))
    plt.plot(list(df["Step"]), list(df["Density"]))
    plt.ylim(0, 1.4)
    plt.xlim(0, 1e7)
    plt.xlabel("Timestep (fs)")
    plt.minorticks_on()
    ax.tick_params(axis="x", which="minor", direction="in")
    ax.tick_params(axis="y", which="minor", direction="in")
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(right=True, top=True)
    ax.tick_params(which="minor", right=True, top=True)
    plt.ylabel(r"Density (g/cm$^3$)")
    fig.tight_layout()
    plt.savefig(args.output_directory / "density_profile.png")
    plt.savefig(args.output_directory / "density_profile.pdf")


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "density_profile",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
