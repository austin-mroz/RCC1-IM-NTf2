import argparse
import itertools
import math
import pathlib
import warnings

from tqdm import tqdm

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    perc_in_cage = []
    perc_in_window = []
    perc_near_cage = []
    perc_not_by = []

    for frame_xyz_files in tqdm(args.input_directory.glob("*.xyz")):
        frame_no = frame_xyz_files.name.removesuffix(".xyz").removeprefix("cluster_")
        frame_cluster_file = frame_xyz_files.parent / f"cluster_info_{frame_no}.txt"

        num_in_cage = 0
        num_in_window = 0
        num_near_cage = 0
        num_not_by = 0

        clusters = read_file(frame_cluster_file)
        clusters.columns = ["id", "size", "COM_X", "COM_Y", "COM_Z"]

        xyzs = read_file(frame_xyz_files)
        xyzs.columns = ["id", "cluster", "x", "y", "z"]

        cage_clusters = clusters[clusters["size"] > 100]
        solv_clusters = clusters[clusters["size"] < 100]

        total_solv = 0

        ppf_in_cage = []
        ppf_in_window = []
        ppf_near_cage = []
        ppf_not_by = []

        for cidx, cage_cluster in cage_clusters.iterrows():
            total_solv = 0
            num_not_by = 0
            num_in_window = 0
            num_near_cage = 0
            num_in_cage = 0
            for sidx, solv_cluster in solv_clusters.iterrows():
                total_solv += 1
                if (
                    math.dist(
                        np.asarray(cage_cluster[["COM_X", "COM_Y", "COM_Z"]]),
                        np.asarray(solv_cluster[["COM_X", "COM_Y", "COM_Z"]]),
                    )
                    < 20.0
                ):
                    cage_xyzs = xyzs[xyzs["cluster"] == int(cage_cluster["id"])]
                    solv_xyzs = xyzs[xyzs["cluster"] == int(solv_cluster["id"])]

                    dists = _get_dist(
                        np.asarray(cage_cluster[["COM_X", "COM_Y", "COM_Z"]]),
                        np.asarray(solv_xyzs[["x", "y", "z"]]),
                    )

                    if (dists < 3).sum() > 0:
                        num_in_cage += 1
                    elif (dists <= 5.5).sum() > 0:
                        num_in_window += 1
                    else:  # (dists < 7).sum() > 0:
                        num_near_cage += 1
                else:
                    num_not_by += 1
            ppf_in_cage.append(100 * (num_in_cage / 162))
            ppf_in_window.append(100 * (num_in_window / 162))
            ppf_near_cage.append(100 * (num_near_cage / 162))
            ppf_not_by.append(100 * (num_not_by / 162))

        perc_in_cage.append(np.mean(ppf_in_cage))
        perc_in_window.append(np.mean(ppf_in_window))
        perc_near_cage.append(np.mean(ppf_near_cage))
        perc_not_by.append(np.mean(ppf_not_by))

        plt.plot(perc_in_cage, "r", label="in cage")
        plt.plot(perc_near_cage, "b", label="near neighbor")
        plt.plot(perc_in_window, "k", label="window")
        plt.plot(perc_not_by, "g", label="not")
        plt.ylim([-5, 100])
        plt.legend()
        plt.savefig(args.output_directory / "percentage_solvent_locations_20A.pdf")
        plt.clf()
        plt.plot(perc_in_cage, "ro", label="in cage")
        plt.plot(perc_in_window, "ko", label="window")
        plt.ylim([-5, 5])
        plt.legend()
        plt.savefig(args.output_directory / "full_percentage_solvent_locash_20A.pdf")
        plt.clf()

    df = pd.DataFrame(
        {
            "perc_in_cage": perc_in_cage,
            "perc_near_cage": perc_near_cage,
            "perc_in_window": perc_in_window,
            "perc_not_by": perc_not_by,
        }
    )
    df.to_csv(args.output_directory / "DATAFRAME.csv")


def _get_dist(cage_com, solv_atoms):
    dists = []
    for solv_atom in solv_atoms:
        dists.append(math.dist(cage_com, solv_atom))
    return np.asarray(dists)


def read_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    df = pd.DataFrame()
    for line in lines[2:]:
        line = line.split(" ")
        nrow = pd.DataFrame(
            [
                [
                    int(line[0]),
                    int(line[1]),
                    float(line[2]),
                    float(line[3]),
                    float(line[4]),
                ]
            ]
        )
        tempdf = pd.concat([df, nrow])
        df = pd.DataFrame(tempdf)

    return df


def analyse_clusters(cluster_file_path, xyz_file_path, cutoff=3):
    def read_file(filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
        df = pd.DataFrame()
        for line in lines[2:]:
            line = line.split(" ")
            nrow = pd.DataFrame(
                [
                    [
                        int(line[0]),
                        int(line[1]),
                        float(line[2]),
                        float(line[3]),
                        float(line[4]),
                    ]
                ]
            )
            tempdf = pd.concat([df, nrow])
            df = pd.DataFrame(tempdf)

        return df

    def do_cages_fit_inside_eachother(other_cages, cage_COM, cutoff=3):
        # check that none of the COMs are too close
        other_cage_coords = np.asarray(other_cages[["x", "y", "z"]])
        cage_COM_coords = np.asarray(cage_COM[["COM_X", "COM_Y", "COM_Z"]])

        com_distance_matrix = np.sort(
            euclidean_distances(other_cage_coords, cage_COM_coords).flatten()
        )

        if min(com_distance_matrix) < cutoff:
            # print('COMS OVERLAPPING WITH CAGE R GROUPS')
            return True
        else:
            return False

    def do_solmols_fit_inside_cages(clusters, cutoff=3):
        cluster_coms = np.asarray(clusters[["COM_X", "COM_Y", "COM_Z"]])
        # check that none of the COMs are too close
        com_distance_matrix = np.sort(
            euclidean_distances(cluster_coms, cluster_coms).flatten()
        )
        com_distance_matrix = com_distance_matrix[com_distance_matrix != 0]
        if min(com_distance_matrix) < cutoff:
            # print('COMS OVERLAPPING')
            # solmols fit inside cages
            return False
        else:
            return True

    clusters = read_file(cluster_file_path)
    clusters.columns = ["id", "size", "COM_X", "COM_Y", "COM_Z"]

    xyzs = read_file(xyz_file_path)
    xyzs.columns = ["id", "cluster", "x", "y", "z"]

    cage_clusters = clusters[clusters["size"] > 100]

    cage_xyzs = xyzs[xyzs["cluster"].isin(list(cage_clusters["id"]))]

    cage_cage_porous = False
    # CHECK IF CAGES FIT INSIDE EACHOTHER
    for cage_num in range(1, max(set(list(cage_xyzs["cluster"]))) + 1):
        cage = cage_xyzs[cage_xyzs["cluster"] == cage_num]
        other_cages = cage_xyzs[cage_xyzs["cluster"] != cage_num]

        cage_COM = cage_clusters[cage_clusters["id"] == cage_num]

        if do_cages_fit_inside_eachother(other_cages, cage_COM, cutoff=3):
            # cages fit inside eachother
            # print('OHHHHH NOOOOOOOOOOOOO ... NOT POROUS WITH EACH OTHER')
            cage_cage_porous = False
        else:
            cage_cage_porous = True

    cage_sol_porous = do_solmols_fit_inside_cages(clusters, cutoff=3)

    if cage_sol_porous and cage_cage_porous:
        # porous liquid maintains porosity
        return {
            "total": True,
            "cage_sol": cage_sol_porous,
            "cage_cage": cage_cage_porous,
        }
    else:
        # porous liquid isn't porous
        return {
            "total": False,
            "cage_sol": cage_sol_porous,
            "cage_cage_porous": cage_cage_porous,
        }


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "thorogouh_percentage_occupation",
    )

    parser.add_argument(
        "-i",
        "--input_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "1_output" / "clusters",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
