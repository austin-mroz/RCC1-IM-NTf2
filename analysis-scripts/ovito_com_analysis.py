import argparse
import itertools
import pathlib
import warnings
from dataclasses import dataclass

from tqdm import tqdm

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import distances, rdf
from sklearn.metrics.pairwise import euclidean_distances


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    cluster_xyz_paths = args.input_directory.glob("*.xyz")
    print(args.input_directory)

    frames = []
    answers = {}

    for frame_xyz_files in tqdm(args.input_directory.glob("*.xyz")):
        frame_no = frame_xyz_files.name.removesuffix(".xyz").removeprefix("cluster_")
        frame_cluster_file = frame_xyz_files.parent / f"cluster_info_{frame_no}.txt"

        porous = analyse_clusters(frame_cluster_file, frame_xyz_files, cutoff=3)
        print(f"{frame_no}: {porous}")

        if not porous["total"]:
            frames.append(frame_no)
            answers[frame_no] = porous

    print(answers)
    print(f"check frame numbers: {frames}")


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
        default=_get_output_directory() / "2_output",
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
