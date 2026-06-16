# Gadget Format-2 I/O, Parallel Strategy, and Code-Unit Conversion

This document describes 
*(a)* the layout of the Gadget *Format-2* snapshot file as used by OpenGADGET3
*(b)* how code units are translated to physical / CGS units in a cosmological context.

## 1. Gadget Format-2 file layout

Format-2 is the *labelled* format: each data block is preceded by a small label block that
identifies it by a 4-character tag. The whole structure is wrapped in Fortran-style record markers (a 4-byte integer giving the size of the record both **before and after** the record itself).

### 1.1 Per-block on-disk layout

For every block (label + data):

```
   ┌──────────────────────────────────────────┐
   │  4 bytes : Fortran record marker  = 8    │   ← size of the label record
   │  4 bytes : ASCII tag, e.g. "POS "        │   ← block identifier, space-padded
   │  4 bytes : uint32  nextblock             │   ← bytes from here to start of next label
   │  4 bytes : Fortran record marker  = 8    │
   ├──────────────────────────────────────────┤
   │  4 bytes : Fortran record marker  = N    │   ← N = size of the data payload in bytes
   │  N bytes : the actual data               │
   │  4 bytes : Fortran record marker  = N    │
   └──────────────────────────────────────────┘
```

The very first record of the file is always the **HEAD** label, then the 256-byte header itself. After that, data blocks follow in an implementation-defined order; readers must therefore *seek by tag*, not by position.

Notes:
- `nextblock` in the label is `size_of_data_payload + 2*4` (the two
  Fortran markers), so it can be used to jump straight to the next
  label without parsing the data.
- The Fortran record markers are reflective; comparing the leading
  one against the value `8` (the expected size of the label record) is
  what OpenGADGET3 uses to **detect endianness**
  (`read_ic.cpp:1703-1722`): if it is not 8, the file is byte-swapped
  and `swap_file` is set so every subsequent read is swapped on the
  fly.

### 1.2 The 256-byte header

Canonical layout (matches `header` struct):

| offset | bytes | C type        | field                              | meaning                                                   |
|-------:|------:|---------------|------------------------------------|-----------------------------------------------------------|
|      0 |    24 | uint32 [6]    | `npart[6]`                         | particles per type **in this file**                       |
|     24 |    48 | double[6]     | `mass[6]`                          | mass per type; 0 if a `MASS` block is present             |
|     72 |     8 | double        | `time`                             | scale factor `a` (cosmological) or time (non-cosmological)|
|     80 |     8 | double        | `redshift`                         | redshift `z = 1/a − 1`                                    |
|     88 |     4 | int32         | `flag_sfr`                         |                                                           |
|     92 |     4 | int32         | `flag_feedback`                    |                                                           |
|     96 |    24 | uint32 [6]    | `npartTotal[6]`                    | total particles per type across all files (low 32 bits)   |
|    120 |     4 | int32         | `flag_cooling`                     |                                                           |
|    124 |     4 | int32         | `num_files`                        | number of files this snapshot is split across            |
|    128 |     8 | double        | `BoxSize`                          | comoving box size, code units                             |
|    136 |     8 | double        | `Omega0`                           | matter density today                                      |
|    144 |     8 | double        | `OmegaLambda`                      | dark-energy density today                                 |
|    152 |     8 | double        | `HubbleParam`                      | h ≡ H₀ / (100 km s⁻¹ Mpc⁻¹)                              |
|    160 |     4 | int32         | `flag_stellarage`                  |                                                           |
|    164 |     4 | int32         | `flag_metals`                      |                                                           |
|    168 |    24 | uint32 [6]    | `npartTotalHighWord[6]`            | high 32 bits, to support > 2³² particles                  |
|    192 |     4 | int32         | `flag_entropy_instead_u`           | 0 → `U` block holds u; 1 → it holds entropy A             |
|    196 |     4 | int32         | `flag_doubleprecision`             | 0 → float fields; 1 → double fields                       |
|    200 |     4 | int32         | `flag_ic_info`                     | IC generator flag                                         |
|    204 |     4 | float         | `lpt_scalingfactor`                | 2LPT scaling                                              |
|    208 |    48 | char[]        | padding / unused                   | brings struct to 256 bytes                                |

Total: **256 bytes**, wrapped in the usual two 4-byte Fortran markers.

Floats vs. doubles is controlled by `flag_doubleprecision`. The width of `ID` blocks can be inferred from `npart * sizeof(IDtype) == record_size` since both 32- and 64-bit IDs are in use.

### 1.3 The `INFO` block

If present, it contains the following structure per every block in the snapshot:

```c
char label[4];       --> the block label
char type[8];        --> the data type
int ndim;            --> how many elements per particle
int is_present[6];   --> which particle types are involved
```

The data type is as follows
`LONG`,`DOUBLE`, `FLOAT`, `LLONG`



### 1.4 Multi-file snapshots

When `header.num_files > 1`, the snapshot lives in files named `<base>.0`, `<base>.1`, ..., `<base>.N-1`. Each file carries its own header with its own `npart[]` (per-file count) and the same `npartTotal[]` / `npartTotalHighWord[]`. The total number of gas particles is

    Ngas_total = npartTotal[0] + (npartTotalHighWord[0] << 32) .



---

## 2. Code units → physical → CGS

Gadget's internal units are defined by three independent base units, chosen at configure time:

```
    UnitLength_in_cm       = 3.085678e21       // 1 kpc / h           (default)
    UnitMass_in_g          = 1.989e43          // 1e10 Msun / h       (default)
    UnitVelocity_in_cm_per_s = 1.0e5           // 1 km/s              (default)
```

From these, every other unit is derived (`allvars.hpp:1267-1269`, internal energies and densities computed at runtime):

```
    UnitTime_in_s       = UnitLength_in_cm / UnitVelocity_in_cm_per_s
    UnitDensity_in_cgs  = UnitMass_in_g    / UnitLength_in_cm^3
    UnitEnergy_in_cgs   = UnitMass_in_g    * UnitVelocity_in_cm_per_s^2
    UnitPressure_in_cgs = UnitMass_in_g    / UnitLength_in_cm / UnitTime_in_s^2
```

Cosmological snapshots use **comoving** lengths and **comoving** densities, with the Hubble parameter `h` baked into the units (the `/h` in the defaults above).
Let

```
    a = header.time                      // scale factor
    z = header.redshift = 1/a - 1
    h = header.HubbleParam
```

### 2.1 Density

Snapshot stores `ρ_code` in *comoving* units. To get the **physical, CGS** mass density used by the cooling tables:

```
    ρ_phys_cgs = ρ_code * UnitDensity_in_cgs * h^2 / a^3
```

The `h^2 / a^3` factor unwinds the `(kpc/h)^-3` comoving volume into a physical `cm^-3`. Equivalently, the hydrogen number density is

```
    n_H = X_H * ρ_phys_cgs / m_p
```

with `X_H ≈ 0.76` (or computed from the `Zs` array if it carries hydrogen explicitly).

### 2.2 Internal energy / entropy / temperature

The `U` block contains either specific internal energy `u` (when `flag_entropy_instead_u == 0`) or entropy
`A = P / ρ^γ = (γ-1) u / ρ^{γ-1}` (when the flag is 1). 
In the entropy case, recover `u` first:

```
    u_code = A_code * ρ_code^{γ-1} / (γ - 1) ,    γ = 5/3
```

The conversion `u_code → u_phys` is **independent of the scale factor** in Gadget's convention because velocities (and therefore `u`, which has units of velocity²) are stored as peculiar quantities in the comoving frame, with no extra `a` factor for gas internal energy:

```
    u_phys_cgs = u_code * UnitVelocity_in_cm_per_s^2     // erg / g
```

Temperature from `u`:

```
    T  = (γ - 1) * μ * m_p * u_phys / k_B
```

The mean molecular weight `μ` depends on the ionization state. Given the snapshot's electron abundance `n_e/n_H ≡ NE` and the hydrogen mass fraction `X_H`,

```
    1/μ  =  X_H * (1 + NE)  +  Y_He / 4 * (1 + Y_He_ionized_state)
```

Actually the block `TEMP` handles the Temperature in K, already.



## 3. Stellar fields

| Block name | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| `AGE`      | contains the value of the expansion factor `a` at which a star particle is born; effectively, this translate in a look-back time once you assume a Cosmological model, and hence you know the age of a star particle. |
| `MASS`     | the value of the star particle’s **mass** *at the moment of the snapshot*; that account for the evolution of the SSP and hence it monothonically decreases due to the mass loss |
| `iM`       | the value of the star particle’s **mass** *at birth*         |
| `Zs`       | the metal masses of the gas from which the star particle was born |

Note: to determine the metallicity of a star particle you must use `iM`:
$$
Z_\star^{mass} = \sum_{\substack{Z_i \ne He}} m_{Zi}\\
Z_\star = Z_\star^{mass} / {iM_\star}
$$

## 4. Metals

The block `Zs` contains the array of metals for every eligible particle type (`gas` and `stars`).
They are saved as *masses*, but for the `H` which can be inferred by
$$
M_H = M - \sum_{Zi} m_{Zi}
$$
