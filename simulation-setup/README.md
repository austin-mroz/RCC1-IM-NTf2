# Simulation setup

## step_0_lattice_setup

The PIL simulation box was constructed using a combination of existing
softwares. The individual systems -- cage and anion -- were manually atom
typed. $NTf_2^-$ anions were atom-typed according to the [IL.ff]() forcefield
presented by the [Padua Group](). IL.ff is an OPLS-AA forcefield parameterized
for ILs. RCC1-IM was atom typed using a combination of IL.ff and OPLS-AA to
describe the charged and neutral species, respectively. Specific atom types
used for each of the PIL components are detailed in Figure 1, below.

The atom types for each component are detailed in the `*.ff` files of the
independent component directories. LAMMPS data files for a system containing
1 cage and 6 anions are generated using a combination of [fftool]() from the
Padua Group.

The generation procedure for the cage is:

```bash
    $ fftool 1 re_transformed.mol -b 32
    $ packmol < pack.inp
    $ fftool 1 re_transformed.mol -b 32 -l
```

The generation procedure for $NTf_2$ is:

```bash
    $ fftool 1 tf2n.xyz -b 8
    $ packmol < pack.inp
    $ fftool 1 tf2n.xyz -b 8 -l
```

These procedures generate a `data.lmp` file for the independent cage and anion.

We then use [moltemplate](https://www.moltemplate.org/) to generate the lattice
in the following way:

First, we use `ltemplify.py` to generate the `.lt` files for the independent
cage and anion. These are necessary to build the PIL lattice.

```bash
    $ python $MOLTEMPLATE_GITHUB/moltemplate/ltemplify.py -name rcc1 -molid "1" in.lmp data.lmp > cage.lt
    $ python $MOLTEMPLATE_GITHUB/moltemplate/ltemplify.py -name tfn2 -molid "2" in.lmp data.lmp > tfn2.lt
```

We then move the `*.lt` files to the `moltemplate_cage_10_ion_60` directory and
generate the `*.lt` files necessary to build the lattice using a procedure
that builds on the individual components. These files are built from the
"ground up"; each dependent on the ones previous. The hierarchy is as follows:

1. `tfn2.lt` -- individual anion
2. `tfn2_6.lt` -- lattice containing 6 anions (1)
3. `cage.lt` -- individual cage
4. `neutral_cage_1_tfn2_6.lt` -- 1 cage and 6-anions (2 and 3)
5. `system.lt` -- lattice of 48 neutral_cage_1_tfn2_6 lattices (4)

The lammps input files are then generated using `system.lt`

```bash

    $ moltemplate system.lt

```

The data files and input scripts for the density equilibration are described:

| INPUT     | TIMESTEP | RUN TYPE | TEMPERATURE | RUN STEP | TOTAL TIME (ns) | TOTAL TIME |
| --------- | -------- | -------- | ----------- | -------- | --------------- | ---------- |
| lattice_0 | 0.001    | npt      | 363         | 50000    | 5e-5            | 50 fs      |
| lattice_1 | 0.01     | npt      | 363         | 100000   | 0.001           | 1000 fs    |
|           | 0.1      | npt      | 363         | 50000    | 0.005           | 5000 fs    |
| lattice_2 | 1.0      | npt      | 363         | 900000   | 0.9             | 900000 fs  |


