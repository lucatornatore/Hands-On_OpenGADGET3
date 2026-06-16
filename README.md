# OpenGADGET3 Hands-On on Cooling, Star Formation and Chemical Enrichment

**gPLUTO / OpenGADGET3 Joint Symposium**

**Heidelberg, 15th-19th June 2026**.
**Instructor: Luca Tornatore, luca.tornatore@inaf.it**

This repository contains the materials used for the hands-on session about running cosmological simulations with Cooling, Star Formation and Chemical Enrichment using the OpenGADGET3 code.

> [!NOTE]
>
> The training event has been organized and funded in the frame of the [EU SPACE Center of Excellence](https://www.space-coe.eu/), *grant agreement No 101093441* and kindly hosted by  [*Haus der Astronomie*](https://www.haus-der-astronomie.de/), Heidelberg. 

## Welcome

In this session we will learn how to configure, compile and run the code, and how to visualize some basic physical outputs.
The user case is a small cosmological box, of size 5Mpc/h with $2\times 32^3$ particles that can be run on a laptop.

We will configure the code for different physical setup, and run it assessing the impact of the different physical processes and options.

A general introduction can be found in `Documentation/OpenGADGET3_cooling_SF_enrichment_lecture_note.md`.



## 1. Contents of this repo

The content of this repository is the following:

| Folder          | sub folder    | Content                                                      |
| --------------- | ------------- | ------------------------------------------------------------ |
| `Documentation` |               |                                                              |
| `Files`         | `Configs`     | Compile-time configuration files                             |
|                 | `ICs`         | Initial Conditions                                           |
|                 | `Input_files` | Run-time input files<br />`common_files.tar` == input files for all runs<br />`SEv_files.tar` == input files for Stellar Evolution |
| `Runs`          |               | The folders tree of all the runs - see below                 |
| `Utilities`     |               | Python routines for producing the plots used for the discussion |

Additional files can be found [HERE](https://drive.google.com/drive/folders/1VIV42mACVdQQzOft3z6mqOjD5kk3d3x9?usp=sharing), where you can download the following materials which are not on this git to limit its overall size:

- the metal-dependent Cooling Tables $\rightarrow$ **needed to run simulations with metal-dependent cooling**
- $z=0$ snapshots from higher resolution runs;
- the initial conditions for runs with higher resolutions, with $2\times 64^3$ and $2\times 128^3$ particlesù

## 2. The folders tree below the `Runs` folder

The contents of the `Runs` folder is as follows

| **Folder**              | sub folder         | Nr. of Part.    | Box (Mpc/h) |   Physics    |                   |              |              |                       |
| ----------------------- | ------------------ | --------------- | :---------: | :----------: | :---------------: | :----------: | :----------: | :-------------------: |
|                         | `box_`             |                 |             |  `COOLING`   | `METAL` `COOLING` |     `SF`     |    `SEv`     |        `WINDS`        |
| `COOL`                  | `32p_5Mpc`         | $2\times32^3$   |      5      | $\checkmark$ |                   |              |              |                       |
| `COOL_SF`               | `32p_5Mpc`         | $2\times32^3$   |      5      | $\checkmark$ |                   | $\checkmark$ |              |                       |
| `COOL_SF_WINDS`         | `32p_5Mpc`         | $2\times32^3$   |      5      | $\checkmark$ |                   | $\checkmark$ |              |     $\checkmark$      |
|                         | `64p_5Mpc`         | $2\times 64^3$  |      5      | $\checkmark$ |                   | $\checkmark$ |              |     $\checkmark$      |
| `ZCOOL_SEv`             | `32p_5Mpc`         | $2\times 32^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ |                       |
|                         | `64p_5Mpc`         | $2\times 64^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ |                       |
| `ZCOOL_SEv_WINDS`       | `32p_5Mpc`         | $2\times 32^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ |     $\checkmark$      |
|                         | `64p_5Mpc`         | $2\times 64^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ |     $\checkmark$      |
|                         | `128p_5Mpc`        | $2\times 128^3$ |      5      |              |   $\checkmark$    |              | $\checkmark$ |     $\checkmark$      |
|                         | `64p_10Mpc`        | $2\times 64^3$  |     10      |              |   $\checkmark$    |              | $\checkmark$ |     $\checkmark$      |
|                         | `128p_20Mpc`       | $2\times 128^3$ |     20      |              |   $\checkmark$    |              | $\checkmark$ |     $\checkmark$      |
|                         | `resolution_plots` |                 |             |              |                   |              |              |                       |
|                         | `volume_plots`     |                 |             |              |                   |              |              |                       |
| `ZCOOL_SEv_STRONGWINDS` | `32p_5Mpc`         | $2\times 32^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ | $\checkmark$ *strong* |
|                         | `64p_5Mpc`         | $2\times 64^3$  |      5      |              |   $\checkmark$    |              | $\checkmark$ | $\checkmark$ *strong* |



## 3. Content of every `box_` folder

In every `box_${NP}p_${L}Mpc` folder you find:

- the `sfr.txt` file, that logs the star formation rate history (not present in the `COOL/box_32p_5Mpc` folder)
- the `paramfile_box_${NP}p_${L}Mpc`: the run-time parameter files used for that run
- the `plots/` folder, in which you can already find the **relevant plots** and the final snapshot from the simulations and the log file of the star formation history:

| File                                        | Description                                                  | Notes           |
| ------------------------------------------- | ------------------------------------------------------------ | --------------- |
| `phasespace_with_marginals.mass.png`        | the $\rho-T$ phasespace color-coded by the mass fraction in every bin | all runs        |
| `sfr.png`                                   | the star formation rate history                              | SF and SEv runs |
| `phasespace_with_marginals.metals.png`      | the $\rho-T$ phasespace color-coded by the *metals* mass fraction in every bin | SEv runs        |
| `phasespace_with_marginals.metallicity.png` | the $\rho-T$ phasespace color-coded by the *mass-weighted metallicity* of every bin | SEv runs        |
| `gas_metals_rho.png`                        | the mass-weighted metallicity $\rho$ distribution of the gas | SEv runs        |
|                                             |                                                              |                 |
| `gas_metals_temp.png`                       | the mass-weighted metallicity $T$ distribution of the gas    | SEv runs        |
| `metallicity_both.cmp.png `                 | the mass distribution of gas and stars *per metallicity*     | SEv runs        |

