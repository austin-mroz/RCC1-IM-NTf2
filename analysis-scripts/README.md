<style>
.circular_image {
  width: 200px;
  height: 200px;
  border-radius: 50%;
  overflow: hidden;
  background-color: blue;
  
  display:inline-block;
  vertical-align:middle;
}
.circular_image img{
  width:100%;
}
</style>

<h1 align="center">
    <br>
    <img class="circular_image" src="../imgs/starting_lattice2.png" alt="RCC1" width="200">
    <br>
    RCC1-IM-NTf2
    <br>
</h1>

<h4 align="center">Construction of an organic cage-based porous ionic liquid using an aminal tying strategy</h4>
<br>

The results presented in this work can be regenerated using the scripts in this repository.

Here, we describe them in the order they should be run:

| SCRIPT | DESCRIPTION |
| ------ | ----------- |
| `plot_density_equilibration_profile.py` | plots the density across the production trajectory |
| `rdf_analysis.py` | compures the radial distribution function of the production run |
| `pywindow_pore_analysis.py` | compute the solo cage pore distribution give the .xyz files from the solo cage simulation |
| `pore_analysis.py` | compares the distributions of the cage cavities from the full liquid simulation and a single cage. The full liquid simulation was analyzed using Zeo++, while the solo cage was analyzed using PyWindow. |
| `ovito_com_analysis.py` | extract all cage and solvent center-of-masses (COMs) over the course of the trajectory | 
| `ovito_percent_occupation_analysis_thorough_algo.py` | given the COMs from `ovito_com_analysis.py`, compute the percent occupation for the solvents |

Several results files are necessary:

| RESULTS | DESCRIPTION |
| ------- | ----------- |
| `solo_cage_pore_diameters.csv` | the pore diameters for each time step of the solo cage simulation | 
| `zeo_pp_frame_analysis.txt` | the Zeo++ analysis results for the full liquid simulation |
| `10_ns_npt_equilibration_density_profile.csv` | the 10 ns equiliration run to verify density converged |