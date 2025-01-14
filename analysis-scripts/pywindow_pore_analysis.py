import argparse
import itertools
import pathlib
import warnings
from dataclasses import dataclass

from tqdm import tqdm

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pywindow as pw


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    frame_nos = []
    diameters = []

    solo_xyz_paths = args.input_directory.glob("*.xyz")
    for solo_xyz in tqdm(solo_xyz_paths):
        print(solo_xyz)
        frame_no = solo_xyz.name[19:-4]
        pw_conformer = pw.MolecularSystem.load_file(solo_xyz)
        pw_cage = pw_conformer.system_to_molecule()
        pw_cage.calculate_pore_diameter()
        diameters.append(pw_cage.properties["pore_diameter"]["diameter"])
        frame_nos.append(frame_no)

    plt.scatter(frame_nos, diameters)
    plt.savefig(args.output_directory / "solo_cage_pore_diameters.png")
    print(np.mean(np.array(diameters)))
    print(np.std(np.array(diameters)))


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory that the outputs are written to",
        type=pathlib.Path,
        default=_get_output_directory() / "11_output",
    )

    parser.add_argument(
        "-i",
        "--input_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "10_output" / "solo_cage_maestro_cif_files",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
