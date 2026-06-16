# OpenGADGET3 — Cooling, Star Formation and Chemical Enrichment
### A hands-on guide

---

## 0. About this note

**What this session is.** A ~1-hour practical in which you configure, compile and run
OpenGADGET3 with the cooling / star-formation / chemical-enrichment machinery switched
on, and inspect how the results change as you flip compile-time options and turn runtime
knobs. The gravity and hydro solvers are assumed known from the preceding sessions; black
hole formation, accretion and feedback come afterwards.

**How to read it.**

- **Part I** is the compact physics framework, one ingredient per section, each ending in
  the relevant switch.
- **Part II** is the consolidated *ingredient → Config flag → runtime parameter → effect* map.
- **Part III** is the hands-on workflow: configure, compile, run.
- **Part IV** is the run campaign and how to read the four plot families you will produce.
- **Part V** is the exercise menu.

**Conventions.** Compile-time switches (set in `Config.sh`, active via `#ifdef`) are written
`LIKE_THIS`. Runtime parameters (set in the parameter file, read at start-up) are written
`LikeThis`. Where a numeric value is given for a runtime parameter it is a *typical* value
from the literature and **must be checked against the parameter file shipped with this
practical** — runtime values are not baked into the binary.



[TOC]



---

# PART I — The framework

## 1. The physics we are playing with

The baseline test case is a small cosmological box:
side length: $5\textrm{Mpc/h}$
particles: $2\times 32^3$

When you see the number of particles stated as above it means that there are $32^3$ *gas* particles and $32^3$  $dark-matter$ particles.

In a 5 Mpc/h box sampled with 32³–128³ gas particles has a gas-particle mass of order 10⁶–10⁸ M⊙ and a gravitational softening of a few kpc. 
The interstellar medium — cold clouds, the star-forming cores, individual supernova remnants — lives orders of magnitude below that and as such we are definitely in a “sub-grid” limit. 

We will apply three families of ingredients, with three master switches:

| Ingredient | What it represents | Master switch |
|---|---|---|
| Radiative cooling | Energy losses of the gas | `COOLING` |
| Star formation + feedback | Conversion of cold gas to stars, and the energy returned | `SFR` |
| Stellar evolution + enrichment | Time-resolved metal production and metal-dependent cooling | `LT_STELLAREVOLUTION` (+ `LT_METAL_COOLING_WAL`) |

The campaign you will run is built around switching these on in sequence and watching the
star-formation history and the metal content respond.

## 2. Cooling

Cooling is what lets baryons dissipate energy, fall to the centres of halos, and become dense enough to form stars. 

**Primordial (H + He).** With `COOLING` on, OpenGADGET3 computes radiative losses for an optically-thin H/He plasma in ionization equilibrium following **Katz, Weinberg & Hernquist (1996)** — collisional excitation and ionization, recombination, free–free (bremsstrahlung) and Compton cooling off the CMB. The cooling function Λ(T) has the familiar shape: a cutoff below ~10⁴ K and two bumps near 10^4.3 and 10^5 K from H and He line cooling.

**UV background.** A spatially uniform, redshift-dependent photoionizing/photoheating background after **Haardt & Madau (1996)** is superimposed. Its main effect is to *photoheat* low-density gas and *suppress* its ability to cool — it sets a temperature floor
in the diffuse IGM and regulates when low-mass halos can accrete and cool gas.

**Metal-line cooling.** Once the gas is enriched, metals dominate the cooling above ~10⁵ K and add a forest of line-cooling channels below it. With `LT_METAL_COOLING_WAL` the code uses the **Wiersma, Schaye & Smith (2009a)** element-by-element tables: cooling rates interpolated in (T, n_H, redshift) for each tracked element, *consistently including the same UV background*. This is the physically important coupling — **the enrichment you compute in Part 5 feeds straight back into the cooling rate**, so cooling and stellar evolution are not independent modules.

> **Operational hook.** `COOLING` alone → primordial cooling + UVB (the `COOL` runs).
> `COOLING` + `LT_METAL_COOLING_WAL` → metal-dependent cooling (the `ZCOOL` runs). The
> comparison `sfr_vs_redshift.baseline__Winds.compare__Cooling.png` isolates exactly this
> difference at fixed winds and SF.

## 3. Star formation — the SH03 effective model

The model is **Springel & Hernquist (2003)** (hereafter SH03), enabled by `SFR`.

The idea: gas denser than a physical threshold `ρ_th` is declared *multiphase* — cold clouds (where stars form) embedded in a hot ambient medium, the two in pressure equilibrium.
Three processes are followed as a set of ODEs for the cold and hot mass fractions:

1. **Star formation** consumes cold clouds on a timescale $t_★$, with
   $dρ_★/dt = (1 − β) ρ_c / t_★$, where $β$ is the fraction promptly returned by short-lived
   massive stars and $t_★ ∝ (ρ/ρ_th)^{−1/2}$, giving a Schmidt-like $ρ_★ ∝ ρ^{3/2}$ law.
2. **Supernova feedback** heats the hot phase and **evaporates** cold clouds (efficiency
   `A₀`, density-dependent).
3. **Radiative cooling** of the hot phase grows the cold clouds back.

These balance into a **self-regulated** cycle with an analytic fixed point. The payoff is an **effective equation of state** `P_eff(ρ)`: star-forming gas is assigned a stiff effective pressure rather than its true (unresolved) thermal state. This is what stops dense gas from fragmenting numerically and keeps the ISM "puffed up" in a controlled way.
The density threshold `ρ_th` is fixed self-consistently by requiring the onset of multiphase behaviour to coincide with the thermal instability / EOS continuity.

Star formation itself is **stochastic**: a gas particle eligible by density spawns (or is converted into) a star particle with a probability matching its SFR over the timestep. `STELLARAGE` records the formation time of each star particle — needed by everything in
Part 5.

> **Operational hook.** Key runtime parameters (parameter file; *typical* values shown,
> confirm against the shipped file):
> `CritPhysDensity` / `CritOverDensity` (→ `ρ_th`), `MaxSfrTimescale` (`t★0`, ~1.5 Gyr),
> `TempSupernova` (~10⁸ K, the effective SN energy), `TempClouds` (~10³ K),
> `FactorSN` (`β`, ~0.1), `FactorEVP` (`A₀`, ~1000). Raising `ρ_th` delays and concentrates
> star formation; shortening `t★0` raises the normalisation of the SFR.

## 4. Galactic winds

SH03 also introduced **kinetic, energy-driven galactic winds**, enabled by `WINDS`. A fraction of star-forming gas is stochastically given a velocity kick `v_w` and (crucially) **temporarily decoupled from hydrodynamics** so it can leave the dense ISM without being stopped by it, recoupling once it reaches low density or a maximum time elapses.

The energy bookkeeping is

$\frac{1}{2}\dot{M}_w \cdot v_w^2 = x\cdot u_{SN}\cdot\dot{M_\star}$

with mass loading $\eta \equiv \dot{M_w}/\dot{M_\star}$.
Canonically
$\eta \sim 2\,,x \sim 1$

 giving

$v_w\sim 480$km/s.

In the stellar-evolution module the wind speed **can** set directly by `LT_WIND_VELOCITY` (350 km/s for the `W` runs, 700 km/s for the `SW` runs).

> **Note.** 
> Fixing `v_w` does *not* by itself fix the energy budget — it depends on whether the implementation holds the **mass loading η fixed** (in which case `SW` injects ~4× the kinetic energy per unit stellar mass of `W`, since energy ∝ v_w²) or holds the **energy fraction χ fixed** (in which case raising `v_w`
> *lowers* η as `η ∝ v_w^(−2)`, ejecting less mass but faster).

## 5. Stellar evolution and chemical enrichment — Tornatore et al. (2007)

This is the part you will spend the most time *playing* with, and it is the reason cooling and stellar evolution are bundled (`ZCOOL+SEv`).
Enabled by `LT_STELLAREVOLUTION`.
Each star particle is a **Simple Stellar Population** (an IMF-weighted ensemble born at one time). 
A star of mass *m* releases its ejecta only after its **lifetime** `τ(m)` — so enrichment is *time-resolved*, and different elements arrive on different clocks. This is what makes abundance *ratios* (not just total Z) closer to the real dynamics.

**The three channels.** 

| Channel | Progenitors | Timing | Principal products |
|---|---|---|---|
| **SNII** <br />(core-collapse) | massive stars, *m* ≳ 8 M⊙ | prompt (~Myr–tens of Myr) | α-elements: O, Mg, Si, Ne, S, + some Fe |
| **SNIa** | WD in a binary | delayed (~10⁸–10⁹⁺ yr) | Fe-peak: Fe, Ni |
| **AGB / LIMS** | low/intermediate mass, ~1–8 M⊙ | delayed (~10⁸–10¹⁰ yr) | C, N, s-process; bulk mass return |

The *delay* of SNIa relative to SNII is what produces the classic [α/Fe]–[Fe/H] "knee": early gas is α-enhanced (SNII only), and [α/Fe] drops once SNIa switch on. None of that exists under IRA.
The SNIa rate requires a **delay-time distribution** (or an explicit binary model, e.g. Greggio & Renzini); SNII and AGB timings follow directly from `τ(m)`.

**The three things you can change.** 
This is the experiment:

1. **The IMF** `φ(m)`. Salpeter, Kroupa, Chabrier, top-heavy, … Sets the *relative number*
   of massive vs low-mass stars → the SNII/SNIa balance, the total metal yield per unit mass
   formed, and the overall [α/Fe] level.
2. **The lifetime function** `τ(m)`. `LT_PM_LIFETIMES` selects **Padovani & Matteucci
   (1993)**; alternatives exist. Sets *when* each channel fires — i.e. the shape of the
   enrichment history at fixed IMF.
3. **The yield tables** `Y(m, Z)` — ejected mass per element as a function of stellar mass
   and initial metallicity, separately for SNII, SNIa and AGB. Sets *how much* of each
   element is made.

**The tracked metals.** `LT_NMET=N` means a (N-1)-element abundance vector is carried per particle (He, and, typically, the metals that matter for cooling and for abundance diagnostics: C, N, O, Mg, Si, S, Ca, Fe, …).
`H` is explicitly tracked from SN and LIMS, but not in the elements array (you can recover it just by difference from the article’s mass).

Which metals are tracked is specified in the file `metals.dat`, in which the elements are named by their chemical label.

Ejecta are spread over the neighbouring gas using the SPH kernel, and the resulting abundances drive the Metal-Dependent cooling of §2 — closing the loop.


> **Operational hook.** The IMF, lifetime and yield choices live in run-time *parameter files*, not in `Config.sh`. 
>
> The compile flags switch the *machinery* on (`LT_STELLAREVOLUTION`, `LT_NMET`, `LT_PM_LIFETIMES`); the *choices* you experiment with are the IMF and yield tables handed to you.
> Swapping a top-heavy for a Salpeter IMF, or one SNII yield set for another, changes the metallicity normalisation and the α/Fe pattern without recompiling.
>
> Main compilation flags that control the behaviour of spreading over gas particles:
>
> | Compilation Flag                  | Default | Meaning                                                      |
> | --------------------------------- | ------- | ------------------------------------------------------------ |
> | `USE_TRUE_NGB`                    | `ON`    | **ON:** the real count of neighbouring gas particles as defined in the run-time parameter file;<br />**OFF:**the number of neighbours is weightd by the SPH kernel. |
> | `LT_USE_KERNEL_IN_WEIGHT`         | `OFF`   | **OFF:** every particle receive a fraction of energy and metals proportional to its mass;<br />**ON:** the fraction of received energy and metals is proportional the the SPH-Kernel-weighted contribution to density |
> | `ALL_EXCLUDE_WINDS_FROM_FEEDBACK` | `OFF`   | **ON:** Wind particles are not receiving contributions from SNs |
> | `ALL_EXCLUDE_WINDS_FROM_COOLING`  | `ON`    | **ON:** Wind particles are not cooling down                  |



# PART II — The ingredient → flag → parameter map

Read this as the operational core. "Effect" is *what moves in the plots* when you change it.

| Physics ingredient | Compile flag(s) (`Config.sh`) | Runtime knob(s) (parameter file) | What changes when you vary it |
|---|---|---|---|
| H/He cooling + UVB | `COOLING` | UVB table path: `TREECOOL` | Enables gas to cool at all; sets IGM temperature floor |
| Metal-line cooling | `LT_METAL_COOLING_WAL` | (Wiersma table path) | Faster cooling of enriched gas → more/earlier SF; couples to enrichment |
| Star formation (SH03) | `SFR` | `CritPhysDensity`, `MaxSfrTimescale`, `TempSupernova`, `FactorSN`, `FactorEVP` | SFR normalisation & threshold; effective EOS stiffness |
| Star ages | `STELLARAGE` | — | Required for time-resolved enrichment; records formation time |
| Enrichment age | `LT_ZAGE`, `LT_ZAGE_LLV`, `LT_LOGZAGE` | — | for every particles you get a Zmass-weighted age of enrichment |
| Galactic winds | `WINDS` | `WindEfficiency` (η), `WindEnergyFraction` (χ), `WindFreeTravel…` | Suppresses & delays SF; ejects/redistributes gas & metals |
| Wind speed | `LT_WIND_VELOCITY=350/700` | (or via η, χ) | Strength of suppression; reach of metal redistribution |
| Stellar evolution | `LT_STELLAREVOLUTION` | IMF file, `τ(m)`, yield tables, SFs file | Turns on time-resolved 3-channel enrichment |
| # tracked species | `LT_NMET=16` | — | How many elements carried per particle |
| Lifetimes | `LT_PM_LIFETIMES` | —<br />(lifetime table if you want a custom one) | *When* each channel releases metals |
| Feedback temperature ceiling | `LT_MAX_TEMP_FEEDBACK=5e8` | — | Caps the temperature reachable by feedback heating |
| Neighbour search for spreading | `LT_USE_TRUENGB`, `LT_USE_KERNEL_IN_WEIGHT`, `LT_STARS_GUESSHSML` | (#neighbours) | Numerical spreading of ejecta; not physics per se |
| Grid for long-range gravity | `PMGRID=128 / 256` | — | Resolution of the PM force; scale it with particle number |

---

# PART III — Hands-on

## 6. Get the code and know the layout

You have been given the source tree and five compile configurations. The relevant pieces:

- `Config_XXXX.sh` — the active compile-time configuration (you copy one of the five over this).
- the parameter file (runtime; SF/wind/IMF/yield settings and I/O).
- the initial conditions for the 5 Mpc/h box at your resolution.
- the `Makefile` / build system.

## 7. Configure — anatomy of `Config_XXXX.sh`

A switch is *active* simply by appearing un-commented; prefix with `#` to disable. The five
configurations differ only in a handful of physics lines — everything else (SPH kernel,
TreePM settings, domain decomposition) is held fixed so the comparison is clean.

The physics deltas across the five:

| Config | Adds, relative to the previous | Net physics |
|---|---|---|
| `Config_COOL.sh` | `COOLING`, `STELLARAGE`; `PMGRID=128` | primordial cooling |
| `Config_COOL_SF.sh` | `COOLING`, `SFR`, `STELLARAGE`; `PMGRID=128` | primordial cooling + SH03 SF, **no winds** |
| `Config_COOL_SF_WINDS.sh` | `+ WINDS` | … + galactic **winds** |
| `Config_ZCOOL_SEv.sh` | `+ LT_STELLAREVOLUTION, LT_METAL_COOLING_WAL, LT_NMET=16, LT_PM_LIFETIMES, LT_MAX_TEMP_FEEDBACK, LT_USE_TRUENGB, LT_STARS_GUESSHSML`; `PMGRID=256` | metal cooling + 3-channel enrichment, **no winds** |
| `Config_ZCOOL_SEv_WINDS.sh` | `+ LT_WIND_VELOCITY=350.0` | … + standard **winds** |
| `Config_ZCOOL_SEv_STRONGWINDS.sh` | `LT_WIND_VELOCITY=700.0` | … + **strong** winds |

Note `PMGRID` jumps 128 → 256 between the `COOL_SF` and `ZCOOL_SEv` families: the long-range PM grid should scale with the particle load. If you push a `COOL_SF` run to 64³ you should raise `PMGRID` to match.

To select a configuration:

```bash
cp Config_ZCOOL_SEv_WINDS.sh Config.sh   # pick your physics
```

## 8. Compile

```bash
make EXEC=OG3.xxx clean
make -j Np CONFIG=Config_xxxx.sh EXEC=OG3.xxx 
```

Expect a couple of minutes. Practical points:

- **Dependencies**: MPI, plus FFTW (for `PMGRID`), GSL, HDF5 if you output HDF5. If `make`
  dies immediately it is almost always a missing library path in the `Makefile`'s
  system block, not a code error.
- **The binary is physics-specific.** Every switch in §7 is resolved at compile time via
  `#ifdef`. Changing physics ⇒ re-`make`. Keep one binary per configuration (rename them,
  e.g. `Gadget3_ZCOOL_SEv_W`) so you are never unsure which physics produced a snapshot.
- A clean way to work: build all the binaries you need *first*, then run them in turn.

## 9. The runtime parameter file

Read at start-up; **not** compiled in. Groups you will touch:

- **I/O & timing**: `OutputDir`, `TimeBegin/TimeMax`, `OutputListFilename`, snapshot cadence.
- **Cosmology & box**: `BoxSize` (5000 kpc/h for the 5 Mpc/h box), `Omega0`, `OmegaLambda`,
  `OmegaBaryon`, `HubbleParam`.
- **SF (SH03)**: `CritPhysDensity`, `CritOverDensity`, `MaxSfrTimescale`, `TempSupernova`,
  `TempClouds`, `FactorSN`, `FactorEVP`.
- **Winds**: `WindEfficiency` (η), `WindEnergyFraction` (χ), `WindFreeTravelLength`/`MaxTime`.
- **Stellar evolution**: paths to the **IMF**, **lifetime** and **yield** tables, and the SNIa
  DTD / binary-fraction settings.

> The IMF and yield files are the levers for the Part-V enrichment exercise. Note which
> directory they live in before you start, and keep a copy of the defaults so you can revert.

## 10. Run

```bash
mpirun -np <N> ./Gadget3_<config> <parameterfile>
```

While it runs, watch the on-the-fly logs — they are your first diagnostic:

- `sfr.txt` — the global SFR vs time. **This is the raw material for every SFR-vs-redshift
  plot.** If winds are working you should see the SFR depressed here in real time.
- `cpu.txt`, `timings.txt`, `energy.txt`, `info.txt` — load balance, step costs, energy
  bookkeeping. (`WRITE_INFO_BLOCK` is on in all five configs, so the INFO block is written.)
- snapshots at the requested output times — the input to the phase-space and metallicity
  plots.

A 32³ / 5 Mpc/h run is a laptop job (minutes to a small number of MPI ranks). 64³ is ~8×
the particles and the work scales worse than linearly; 128³ is a small-cluster job — you are
given the 64³ and 128³ *outputs* precomputed for the convergence study, you are not expected
to run them live.

---

# PART IV — The campaign and how to read the plots

## 11. The run matrix

All boxes are 5 Mpc/h unless a different size is noted. `noW / W / SW` = no winds / 350 / 700 km/s.

| Physics | noW | W (350) | SW (700) |
|---|---|---|---|
| `COOL+SF` | 32³ | 32³ | — |
| `ZCOOL+SEv` | 32³, 64³ | 32³, 64³, 128³ (+ 64³/10 Mpc, 128³/20 Mpc) | 32³, 64³ |

Two axes of variation are deliberately separated:

- **Physics axis** (at fixed 32³): no winds → winds → strong winds; primordial → metal cooling.
- **Numerics axis** (`ZCOOL+SEv+W` only): *resolution* at fixed volume (32³→64³→128³ in 5 Mpc)
  and *volume* at fixed mass resolution (5/10/20 Mpc).

## 12. Star-formation history (SFR vs redshift)

**What it is.** 
Total SFR in the box as a function of redshift (built from `sfr.txt`). The single most informative one-line summary of a run.

**How to plot it.** 
Use the `Utilities/plot_sfr.py` utility provided

**What each comparison isolates:**

- `…baseline__COOL_SF.compare__Winds.png` — winds vs no winds with *primordial* cooling.
  Winds should **suppress and delay** the SFR.
- `…baseline__Winds.compare__Cooling.png` — primordial vs metal cooling at fixed winds+SF.
  Metal cooling should **raise/advance** the SFR (enriched gas cools faster).
- `…baseline__ZCOOL_SEv.compare__Winds.png` (32³) and `…high_resolution…` (64³) —
  noW vs W vs SW with the full chemical model. The wind-strength sequence; expect monotonic
  suppression with wind speed *if* energy is the controlling variable (see §4's caveat).
- `resolution_plots/sfr.baseline__ZCOOL_SEv_W.compare__resolution.png` — same physics
  (`ZCOOL+SEv+W`), 32³/64³/128³. **This is a convergence test, not a physics result.**
- `volume_plots/sfr.baseline__ZCOOL_SEv_W.compare__volume.png` — fixed mass resolution,
  growing volume; tests whether the 5 Mpc box is a fair sample.

## 13. Metallicity distributions of gas and stars

**What they are.** Distributions of metallicity over gas particles and over star particles
(`gas_metallicity.*`, `stellar_metallicity.*`, and the combined `metallicity_both.cmp.png`).

**How to plot**

Use the `Utilities/plot_metallicity_distribution.py`.

**How to read them.** 
Position (median Z), width, and the difference between the gas and stellar distributions. Stars lock in the metallicity of the gas *at their formation time*; gas keeps evolving — so the two distributions encode different epochs of the same history.

**What the comparisons isolate.** The `compare__Winds` versions show how winds move metals: winds pull enriched gas out of star-forming regions, **lowering and broadening** the gas metallicity and redistributing metals to lower density. The `…resolution` and `…volume` versions test convergence of the metal content — historically the *most* resolution-sensitive quantity here, because enrichment depends on resolving the star-forming gas.

## 14. Phase-space diagrams ($ρ–T$), weighted three ways

**What they are.** The temperature–density plane, the workhorse diagram of the gas state, with each pixel coloured by a different weight:

- `phasespace.mass.*` — **mass fraction**: where the gas *is*. You should see the canonical
  branches — cold dense star-forming gas (on the effective EOS), the warm photoionized IGM,
  and shock-heated hot gas.
- `phasespace.metals.*` — **metal-mass fraction**: where the *metals* are.
- `phasespace.metallicity.*` — **mass-weighted mean metallicity** per pixel: how *enriched*
  the gas in each phase is.
- `phasespace_with_marginals.*` — same, with the 1-D ρ and T projections on the axes.
  *(Some 32³ files carry a `marginas` typo in the name — same plot.)*
- `gas_metals_rho.png`, `gas_metals_temp.png` — the metal content projected onto density and
  onto temperature alone.

**How to plot**

Use the `Utiolities/plot_phasespace_with_massdistribution.py`.

**How to read them, operationally.** Compare *mass* vs *metals* weighting: the metals will be
concentrated in a different part of the plane than the bulk mass — that offset *is* the
enrichment story. The effective EOS shows up as a tight diagonal locus in the star-forming
branch. When winds are on, look for metals carried up to lower density / higher temperature
than in the no-wind case.

---

# PART V — Exercise menu

A suggested progression; pick according to time.

1. **Build & first run.** Compile `Config_COOL_SF.sh`, run the 32³ box, watch `sfr.txt`,
   produce the SFR-vs-redshift curve. *Goal: the toolchain works.*

2. **Winds on.** Switch to `Config_COOL_SF_WINDS.sh`, re-`make`, re-run. Overlay the two SFR
   histories. *Goal: see suppression & delay.*

3. **Metal cooling + enrichment.** Switch to `Config_ZCOOL_SEv_WINDS.sh`, re-run. Compare its
   SFR to the primordial-cooling case at fixed winds. Produce the gas and stellar metallicity
   distributions and a `phasespace.metallicity` map. *Goal: see the cooling–enrichment loop.*

4. **Wind strength.** Run `…STRONGWINDS` (700 km/s) and compare noW / W / SW. **State which
   variable is held fixed (η or χ)** before interpreting the trend (§4). *Goal: separate a
   physics trend from a definitional artefact.*

5. **Play with the enrichment ingredients** (no recompile needed): swap the **IMF** (e.g.
   Salpeter ↔ top-heavy) and/or the **yield tables** in the parameter file and re-run the
   `ZCOOL+SEv+W` 32³ box. Predict *first*, then check: a top-heavier IMF → more SNII → higher
   overall Z and higher [α/Fe]; a different SNIa yield set → shifts the Fe normalisation.
   *Goal: connect a model choice to a measurable abundance pattern.*

6. **Convergence, if time.** Use the provided 64³/128³ and larger-volume outputs to ask: which
   of the conclusions from steps 2–5 are resolution- and volume-converged?

---

## References

- Katz N., Weinberg D. H., Hernquist L. 1996, ApJS, 105, 19 — H/He cooling, TreeSPH.
- Haardt F., Madau P. 1996, ApJ, 461, 20 — UV/X-ray background.
- Springel V., Hernquist L. 2003, MNRAS, 339, 289 — effective multiphase SF model and winds.
- Padovani P., Matteucci F. 1993, ApJ, 416, 26 — stellar lifetimes.
- Tornatore L., Borgani S., Dolag K., Matteucci F. 2007, MNRAS, 382, 1050 — time-resolved
  chemical enrichment.
- Wiersma R. P. C., Schaye J., Smith B. D. 2009a, MNRAS, 393, 99 — element-by-element cooling.

---

*Draft for review. Runtime parameter values are typical-literature placeholders pending the
shipped parameter file. Open question flagged in §4 (η vs χ in the wind implementation).*



$\partial_t \rho+\nabla\cdot(\rho\mathbf{v}) = 0$

$\partial_t u + \nabla\cdot[(u+P)\mathbf{v}] = \rho\mathbf{v}\cdot\mathbf{g}$

$\frac{1}{2}\dot{M}_w \cdot v_w^2 = x\cdot u_{SN}\cdot\dot{M_\star}$

$\eta \equiv \dot{M}_w/ \dot{M}_\star$



$\gamma > 4/3$
