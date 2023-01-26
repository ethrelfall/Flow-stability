# Repository for stability analysis of some fluid flows

script zucatti.py: this is an attempt at implementing perturbation analysis of convective flows as in *Assessment of reduced-order modeling strategies for convective heat transfer*, V. Zucatti et al.

In that paper it's shown that there are interesting harmonic perturbations for Rayleigh number $Ra = 3.4 \times 10^5$, shown here by modelling the same system using the *Nektar++* spectral / hp element code.  The figure shows (L-R) fluctuation in temperature, vertical velocity component, horizontal velocity component.

![Zucatti_flucts_nektarpp](png/Zucatti_flucts_nektarpp.png "Time-harmonic fluctuation in the temperature field for $Ra = 3.4 \times 10^5$ found by Zucatti et al and reproduced in *Nektar++*.")

The question is why the Firedrake perturbation analysis does not find an unstable or oscillatory eigenmode in the neighbourhood of $Ra = 3.4 \times 10^5$.  Admittedly the paper using a rather high element order compared to that in the script, but one would expect to see *something* at lower order (right?) ...

A very similar script (farrell.py) was used to find a bifurcation instability in the heated-from-below Rayleigh-Benard convection case, with results roughly agreeing with the analysis in Section IV.A of the paper *Bifurcation analysis of two-dimensional Rayleigh-Benard convection using deflation* by Boulle, Dallas, and Farrell.  Note this is not quite the same type of instability (bifurcation not turbulent transition).

Script vertical_convection_Pr0pt01.py: this is a parameter scan in $Ra$ for $Pr=0.01$ (the latter is a high-thermal conductivity / viscosity ratio, e.g. a liquid metal).  The scan seems to show a transition at $log_{10} Ra=5.95$ where stable eigenmodes become oscillatory (the eigenvalue becomes complex).  The graphs show agreement between runs using order-6 and order-5 elements (note the pressure uses one order fewer) and so it seems the system is converged (note there is no very thin boundary layer for these parameter choices).  Simulations using *Nektar++* seem to show a similar transition between a stationary laminar steady state and a time-dependent quasi-steady state.  Note this script uses the same techniques as in zucatti.py but seems to work as expected (why?).

![Pr0pt01_eigenvalues](png/Pr0pt01_eigenvalues.png "Real and imaginary parts of the (largest imaginary part) eigenvalue for linear perturbation over the laminar convecting state for $Pr=0.01$, apparently showing transition to a turbulent state.")
