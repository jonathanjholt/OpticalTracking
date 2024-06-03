# Optical Tracking using certus

## Overview
Optical tracking consists of two stages.

Define the coordinate systems:
1. Digitisation of landmarks (femur, tibia, patella,...)
2. Definition of bone coordinate system
3. Definition of tracker coordinate system
4. Calculation of bone-to-tracker transform

Process tracked data:
1. Import tracked data
2. Define dynamic tracker positions
3. Calculate dynamic bone positions

## Naming convention
Variables names convey the frame of reference, whether it's a transform or a position vector, the bone and the point in time.

A variable named `gTt0` is trying to convey the following notation:
$$
_gT_{t_0}
$$
which means that in the **global** ($_g$) frame of reference, we define a **transform** ($T$) that describes the **tibia**'s ($t$) frame of reference in the initial instant ($_0$).

Similarly, the variable `Pin1_r_tc` is:
$$
_{\text{Pin1}}r_{t_c}
$$
which means that in **Pin 1**'s ($_{\text{Pin1}}$) frame of reference, we define a position vector ($r$) for the **tibia**'s origin ($_t$).

## Visualising a frame of reference change
Taking an example operation from the code: `Pin1_T_tc = gT_Pin1_t0\gTt0`, the change of frame of reference can be visualised like so:
$$
\begin{align*}
_{\text{Pin1}}T_{t_c} &= {\left[ _gT_{\text{Pin1},{t0}} \right]}^{-1} \cdot \left[_gT_{t_0}\right] \\
    & = \left[ _{\text{Pin1}}T_{g,{t0}} \right] \cdot \left[_gT_{t_0}\right] \\
    & = \left[ _{\text{Pin1}}T_{\cancel{g},{t0}} \right] \cdot \left[\cancel{_g}T_{t_0}\right] \quad \text{\{ Change 0 to C\}} \\
    & = \left[ _{\text{Pin1}}T_{t_c}\right] \\
\end{align*}
$$
