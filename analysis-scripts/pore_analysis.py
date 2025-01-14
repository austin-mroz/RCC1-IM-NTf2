import argparse
import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

warnings.filterwarnings("ignore")


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    print(args.output_directory)

    zeo_pp_df = _get_res_dataframe(pathlib.Path.cwd() / "zeo_pp_frame_analysis.txt")
    zeo_pp_df = zeo_pp_df.sort_values("frame_number")
    # zeo_pp_df = zeo_pp_df[zeo_pp_df["frame_number"] > 500]
    zeo_pp_df = zeo_pp_df[["largest_free_sphere"]]
    zeo_pp_df.columns = ["diameters"]
    zeo_pp_df["simulation"] = "full liquid simulation"

    pywindow_df = pd.read_csv("solo_cage_pore_diameters.csv")
    pywindow_df.columns = ["frame_number", "diameters"]
    pywindow_df.sort_values("frame_number")
    pywindow_df = pywindow_df[["diameters"]]
    pywindow_df["simulation"] = "solo cage"

    df = pd.concat([zeo_pp_df, pywindow_df])
    df_m = pd.melt(
        df, id_vars=["simulation"], var_name="diameters", value_name="values"
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.boxplot(x="simulation", y="values", data=df_m, hue="simulation")
    ax.tick_params(axis="y", which="minor", direction="in")
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(axis="y", direction="in")
    ax.set_ylabel(r"cage diameters ($\AA$)")
    plt.tight_layout()
    plt.savefig(args.output_directory / "diameters_boxplot.png")
    plt.savefig(args.output_directory / "diameters_boxplot.pdf")
    plt.clf()
    palette = {
        "full liquid simulation": "#ACDB7B",
        "solo cage": "#7BA8DB",
    }
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.violinplot(
        x="simulation",
        y="values",
        data=df_m,
        hue="simulation",
        palette=palette,
        split=True,
        inner="quart",
        linecolor="white",
        linewidth=0.5,
    )
    ax.tick_params(axis="y", which="minor", direction="in")
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(right=True, top=False)
    ax.tick_params(which="minor", right=True, top=True)

    ax.set_ylabel(r"cage cavity diameter ($\AA$)")
    plt.tight_layout()
    plt.savefig(args.output_directory / "diameters_violinplot.png")
    plt.savefig(args.output_directory / "diameters_violinplot.pdf")

    print(df)
    print(df_m)
    print(sns.__version__)

    fig, ax = plt.subplots(figsize=(4, 4))
    print(zeo_pp_df["diameters"].shape)
    print(zeo_pp_df["diameters"].isna().sum())
    l = zeo_pp_df["diameters"].values.flatten().tolist()
    sns.kdeplot(zeo_pp_df["diameters"])
    sns.kdeplot(pywindow_df["diameters"])
    # sns.kdeplot(list(np.array(pywindow_df["diameters"])))
    plt.ylim([0, 3.5])
    plt.xlim([0, 4])
    plt.tight_layout()
    plt.savefig(args.output_directory / "diameters_kdeplot.pdf")

    exit()
    plt.scatter(
        list(zeo_pp_df["frame_number"]),
        list(zeo_pp_df["largest_free_sphere"]),
        label="full liquid simulation",
    )

    plt.scatter(
        list(pywindow_df["frame_number"]),
        list(pywindow_df["diameters"]),
        label="solo cage simulation",
    )

    plt.legend()
    plt.savefig(args.output_directory / "pore_analysis.png")


def _get_res_dataframe(total_res_path: pathlib.Path) -> pd.DataFrame:
    df = pd.DataFrame()
    with open(total_res_path, "r") as res_file:
        lines = res_file.readlines()
    for line in lines:
        split_line = line.split(" ")
        if len(split_line[0]) > 0:
            nrow = pd.DataFrame(
                [
                    [
                        int(split_line[0].removesuffix(".res").removeprefix("frame_")),
                        float(split_line[4]),
                        float(split_line[5]),
                        float(split_line[7].removesuffix("\n")),
                    ]
                ]
            )
        else:
            nrow = pd.DataFrame(
                [
                    [
                        int(split_line[1].removesuffix(".res").removeprefix("frame_")),
                        float(split_line[5]),
                        float(split_line[6]),
                        float(split_line[8].removesuffix("\n")),
                    ]
                ]
            )
        tdf = pd.concat([df, nrow])
        df = pd.DataFrame(tdf).reset_index(drop=True)
    df.columns = [
        "frame_number",
        "largest_included_sphere",
        "largest_free_sphere",
        "largest_along_path",
    ]
    return df


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "pore_analysis",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
