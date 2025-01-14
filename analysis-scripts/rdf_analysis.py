import argparse
import pathlib
import warnings
from dataclasses import dataclass

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rdf


@dataclass(frozen=True)
class RDF:
    name: str
    bins: list
    rdf: list


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(
        "dump.npt_equilibration.lammpstrj",
        topology_format="LAMMPSDUMP",
    )

    nt_group = u.select_atoms("type 24")  # NBT
    cr_group = u.select_atoms("type 17")  # CR
    ci1_group = u.select_atoms("type 5")  # CI1
    ci2_group = u.select_atoms("type 6")  # CI2
    obt_group = u.select_atoms("type 22")  # OBT
    oy_group = u.select_atoms("type 11")  # OY

    nt_cr = _plot_rdf(nt_group, cr_group, args.output_directory, "CR-NBT.pdf")
    ci1_nt = _plot_rdf(nt_group, ci1_group, args.output_directory, "CI1-NT.pdf")
    ci2_nt = _plot_rdf(nt_group, ci2_group, args.output_directory, "CI2-NT.pdf")
    nt_nt = _plot_rdf(nt_group, nt_group, args.output_directory, "NBT-NBT.pdf")
    cr_cr = _plot_rdf(cr_group, cr_group, args.output_directory, "CR-CR.pdf")
    ob_cr = _plot_rdf(cr_group, obt_group, args.output_directory, "OBT-CR.pdf")
    oy_ob = _plot_rdf(oy_group, obt_group, args.output_directory, "OBT-OY.pdf")
    nt_oy = _plot_rdf(nt_group, oy_group, args.output_directory, "NBT-OY.pdf")

    fig, ax = plt.subplots(figsize=(4, 4))
    plt.xlim(0.5, 15)
    plt.ylim(0, 5)
    plt.xticks(np.arange(0, 15, 3.0))
    plt.minorticks_on()
    ax.tick_params(axis="x", which="minor", direction="in")
    ax.tick_params(axis="y", which="minor", direction="in")
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(right=True, top=True)
    ax.tick_params(which="minor", right=True, top=True)
    ax.set_ylabel("g(r)")
    ax.set_xlabel(r"r ($\AA$)")

    colors = ["#DB7B7B", "#ACDB7B", "#7BA8DB"]
    i = 0
    for rdf in [nt_cr, ob_cr, nt_oy]:
        plt.plot(
            rdf.bins, rdf.rdf, label=f"{rdf.name.removesuffix('.pdf')}", color=colors[i]
        )
        i += 1
    plt.legend()
    plt.savefig(args.output_directory / "total_rdf.pdf")

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    ax1.plot(nt_nt.bins, nt_nt.rdf, label="NBT-NBT")
    ax1.set_ylabel("g_Nan-Nan (r)")
    ax2.plot(nt_cr.bins, nt_cr.rdf, label="CR-NBT")
    ax2.set_ylabel("g_CR-Nan (r)")
    ax3.plot(cr_cr.bins, cr_cr.rdf, label="CR-CR")
    ax3.set_ylabel("g_CR-CR (r)")
    ax1.set_xlim(0.5, 16)
    ax2.set_xlim(0.5, 16)
    ax3.set_xlim(0.5, 16)
    ax1.set_ylim(0, 2)
    ax2.set_ylim(0, 4)
    ax3.set_ylim(0, 2)
    plt.savefig(args.output_directory / "rdf_lit_match.pdf")

    plt.clf()
    fig, ax = plt.subplots(figsize=(4, 4))
    plt.plot(nt_nt.bins, nt_nt.rdf, label="anion-anion (NBT-NBT)", color="#ACDB7B")
    # plt.set_ylabel("g_Nan-Nan (r)")
    plt.plot(nt_cr.bins, nt_cr.rdf, label="cation-anion (CR-NBT)", color="#7BA8DB")
    # plt.set_ylabel("g_CR-Nan (r)")
    plt.plot(cr_cr.bins, cr_cr.rdf, label="cation-cation (CR-CR)", color="#DB7B7B")
    plt.plot(
        ob_cr.bins,
        ob_cr.rdf,
        label=f"{ob_cr.name.removesuffix('.pdf')}",
        color="#ffd152",
    )
    plt.plot(
        nt_oy.bins,
        nt_oy.rdf,
        label=f"{nt_oy.name.removesuffix('.pdf')}",
        color="#d098fa",
    )
    plt.ylim(0, 5)
    plt.xticks(np.arange(0, 15, 3.0))
    plt.xlim(0.5, 15)
    plt.minorticks_on()
    ax.tick_params(axis="x", which="minor", direction="in")
    ax.tick_params(axis="y", which="minor", direction="in")
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(right=True, top=True)
    ax.tick_params(which="minor", right=True, top=True)
    ax.set_ylabel("g(r)")
    ax.set_xlabel(r"r ($\AA$)")
    plt.legend()
    plt.savefig(args.output_directory / "total_rdf_submission_fig.pdf")


def _plot_rdf(nt_group, c_group, path, name) -> RDF:
    rdfobj = rdf.InterRDF(nt_group, c_group)
    rdfobj.run()

    plt.plot(rdfobj.results.bins, rdfobj.results.rdf)
    plt.xlim(0.5, 15)
    plt.ylim(0, 4)
    plt.savefig(path / name)
    plt.clf()

    return RDF(name=name, bins=rdfobj.results.bins, rdf=rdfobj.results.rdf)


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "rdf_output",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
