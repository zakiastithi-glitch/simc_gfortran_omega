
# Meson Electroproduction Kinematics — C++ Version

The original Fortran version was written by Garth Huber in 2000. Wenliang (Bill) Li converted the code to C++ in 2017.

* Original Fortran author: Garth Huber — [Garth.Huber@uregina.ca](mailto:Garth.Huber@uregina.ca).
* C++ conversion: Wenliang (Bill) Li — [wenliang.li@msstate.edu](mailto:wenliang.li@msstate.edu).

C++ program for calculating exclusive meson electroproduction kinematics for reactions of the form

[
p(e,e'p),\text{meson}
]

The program calculates electron-scattering, virtual-photon, proton, and meson kinematic quantities for user-specified meson mass, (Q^2), and invariant mass (W).

The code can be used for kinematic studies involving mesons such as:

* (\omega): 0.78265 GeV
* (\phi): 1.019461 GeV
* (J/\psi): 3.096916 GeV
* And other mesons

## Files

```text
kinematics.cpp
kinematics.h
```

### `kinematics.cpp`

Contains the main program and the kinematic calculations, including:

* electron scattering energy and momentum
* virtual-photon energy and momentum
* virtual-photon angle
* virtual-photon polarization (\epsilon)
* center-of-mass boost parameters
* proton momentum and angle
* meson momentum and angle
* Mandelstam variables (t) and (u)
* Bjorken (x)

The shared calculation state is private to this source file. Comments explain the units, two-body CM solution, boost to the laboratory frame, and table formatting.

### `kinematics.h`

Contains:

* particle masses
* beam-energy configuration
* angle conversion functions
* predefined proton-angle arrays
* function declarations

## Build

The program can be compiled directly with `g++`:

```bash
g++ kinematics.cpp -o main
```

Optionally, enable compiler warnings:

```bash
g++ -Wall -Wextra kinematics.cpp -o main
```

## Run

Run the executable with:

```bash
./main
```

The program will prompt for:

```text
Enter meson mass, Q^2, W in GeV
```

Enter three numbers separated by spaces:

```text
meson_mass  Q2  W
```

where:

* `meson_mass` is the meson mass in GeV
* `Q2` is (Q^2) in GeV(^2)
* `W` is the invariant mass of the photon-proton system in GeV

For example, a (J/\psi) mass can be entered as:

```text
3.096916  Q2  W
```

with the desired values of `Q2` and `W`.

## Kinematic Calculations

### Virtual photon

The virtual-photon energy is calculated from

[
\nu =
\frac{W^2 + Q^2 - M_p^2}{2M_p}.
]

Its three-momentum magnitude is

[
|\vec q| = \sqrt{\nu^2 + Q^2}.
]

The code then determines the boost between the laboratory and photon-proton center-of-mass frames.

### Electron scattering

For each configured beam energy, the program calculates:

* incoming electron energy and momentum
* scattered electron energy
* scattered electron momentum
* electron scattering angle
* virtual-photon polarization

The virtual-photon polarization is calculated as

[
\epsilon =
\frac{1}
{1 + 2\frac{|\vec q|^2}{Q^2}
\tan^2\left(\frac{\theta_e}{2}\right)}.
]

### Final-state proton and meson

The program calculates the proton-meson momentum in the center-of-mass frame and transforms the final-state particles into the laboratory frame.

Meson momentum is reported as a nonnegative magnitude. Its polar angle relative to the virtual photon is calculated from the boosted momentum components and lies between 0 and 180 degrees; backward mesons have angles greater than 90 degrees.

It determines quantities including:

* proton momentum (P_p)
* proton angle relative to the virtual photon
* meson momentum (P_m)
* meson laboratory angle
* (t)
* (u)
* Bjorken (x)

## Output

The output table contains columns similar to:

```text
Tinc   Te*   THe*   epsilon   -t   THq   Pp   THpq   Pm   Thmq   -u   x
```

The quantities correspond approximately to:

| Column    | Description                                 |
| --------- | ------------------------------------------- |
| `Tinc`    | Incident electron beam energy               |
| `Te*`     | Scattered electron kinetic energy           |
| `THe*`    | Electron scattering angle                   |
| `epsilon` | Virtual-photon polarization                 |
| `-t`      | Negative Mandelstam (t), in GeV²                     |
| `THq`     | Virtual-photon angle                        |
| `Pp`      | Proton momentum                             |
| `THpq`    | Proton angle relative to the virtual photon |
| `Pm`      | Meson momentum                              |
| `Thmq`    | Meson angle                                 |
| `-u`      | Negative Mandelstam (u), in GeV²                     |
| `x`       | Bjorken (x)                                 |

All table angle columns (`THe*`, `THq`, `THpq`, and `Thmq`) are reported in degrees (`deg`). Radian-valued angles are converted to degrees before printing. Internal trigonometric calculations use radians as required by the C++ math functions. Momenta and energies are primarily handled internally in MeV/MeV(c^{-1}), with selected quantities converted to GeV-based units for input and output.

## Beam Energy Configuration

Beam energies are configured in `kinematics.h`.

The current declaration is:

```cpp
const int beam_energy[] = {10730};
```

Thus the explicitly configured beam energy is

```text
10730 MeV = 10.730 GeV
```

To change or add beam energies, modify the `beam_energy` array in `kinematics.h`.

For example:

```cpp
const int beam_energy[] = {
    6600,
    7400,
    7900,
    8800,
    9900,
    10900
};
```

The program automatically determines the number of elements when beam energies are added or removed.

## Proton-Angle Scan

A generic array of proton angles relative to the virtual-photon direction is defined in `kinematics.h`:

```cpp
double theqp_generic[] = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20
};
```

These values are in degrees.

## Special (Q^2) Settings

The source retains commented-out historical proton-momentum settings for:

```text
Q^2 = 1.60 GeV^2
Q^2 = 2.45 GeV^2
```

These fixed momentum overrides are disabled. For every (Q^2) value, the proton laboratory momentum is calculated from the center-of-mass kinematics.

## Development Notes

Historical proton-angle arrays remain in the header as alternative scan configurations. Obsolete commented-out implementations and debugging output have been removed from the source.

The current code should be treated as research/development software. When changing beam energies, angular grids, or kinematic settings, verify that the selected (Q^2), (W), meson mass, and beam energy correspond to a physically allowed kinematic region.

## Repository Setup

To add this README to the Git repository:

```bash
git add README.md
git commit -m "Add README"
git push origin main
```
