# Hands-On



## 0. First steps first

Let’s confirm that everybody:

- has downloaded the container from the previous sessions
- has got at least the relevant files from this repo
  - `ics_32p_5Mpc.?`
  - the configuration files from `Files/Configs/`
  - the run-time input files from `Files/Input_files/`
  - the cooling tables (needed when running simulations with metal-dependent cooling)
  - the run-time parameter file from every `box_${NP}p_${L}Mpc/` folder

### 0.1 Recap

To run OpenGADGET3 we need to:

1) setup the compile-time configuration file
2) compile the code
3) setup the run folder:
   - setup the run-time parameter file
   - prepare all the needed run-time files (depend on the active physics module)
4) successfully run the code



## 1. The first simulations: 

## primordial cooling, SH03 SF, Winds

#### Step 1 - COOL & SF

Our first attempt is running a simulation with:

**Gravity, SPH hydrodynamics, Primordial Cooling and SH03 Star Formation**

* Configuration file: `Config_COOL_SF.sh`
* Let’s comment together the parameter file 
* We need the following files in the run folder:
  * parameter file
  * `outputs.txt`  $\rightarrow$ the list of times at which a snapshot is written
  * `TREECOOL` $\rightarrow$ the UV background

#### Step 2 - Add the WInds

Now, let’s add the **winds**:

* Configuration file: `Config_COOL_SF_WINDS.sh`
* In the run folder we do not need anything more than before

#### Step 3 - Compare

Now, let’s compare:

* the star formation rate histories
* the phase-space plots, color-coded by mass

Let’s inspect together the differences with the `COOL` run: I have not suggested to run it because… ? guess.



## 2. Let’s leap forward with chemical enrichment

#### Step 1 - let’s recap the IMFs first

| ![imfs_0.1-100](/run/user/1000/doc/ffd2fb76/imfs_0.1-100.png) | ![imfs_relative_0.1-100](/run/user/1000/doc/a5869322/imfs_relative_0.1-100.png) |
| :----------------------------------------------------------: | :----------------------------------------------------------: |
|                    *Five different IMFs*                     |     *The relative weight of the IMFs, in Salpeter units*     |

#### Step 2 - Configure the run folder

What we need is something more than before:

- everything we had before
  - paramfile, `outputs.txt`, `TREECOOL`
- the yields files:
  Sn-II, Sn-Ia and LIMs yields
  $\rightarrow$ how to specify the yields files in the parameter files
- the `metals.dat` file: specify the elements to be tracked
- the `IMF` file: specify the Initial Mass Function
- the `SFs` file: specify the Star Formation Modes

Let’s inspect the parameter file for `ZCOOL + SEv + WINDS`.



#### Step 3 - Discussion

Now, let’s have a quick tour of the plots that are already there, and understand the impact of the cooling/metal cooling, the star formation and the SEv, and of the winds.