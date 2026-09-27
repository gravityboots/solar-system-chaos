# Attribution

This project is a research application built on existing open-source software,
published data and a large body of scientific literature. The core numerical
work, including the N-body integrators, variational equations, MEGNO and
relativistic corrections, was done by the authors credited below. This
repository adds configuration, orchestration, diagnostics, visualisation and
experiment design on top of it.

If you use this repository, please also cite the REBOUND and REBOUNDx papers
listed here. Running `sim.cite()` on a configured simulation prints REBOUND's
recommended citations for that exact setup.

## REBOUND

An open-source N-body integration code by Hanno Rein and collaborators.

- Repository: <https://github.com/hannorein/rebound>
- Documentation: <https://rebound.readthedocs.io>
- License: GNU General Public License v3.0
- Version used for the included results: 4.6.0

Citations for the features this project uses:

- **REBOUND:** Rein, H. & Liu, S.-F. (2012). REBOUND: an open-source
  multi-purpose N-body code for collisional dynamics. *A&A* 537, A128.
  doi:10.1051/0004-6361/201118085
- **WHFast integrator:** Rein, H. & Tamayo, D. (2015). WHFast: a fast and
  unbiased implementation of a symplectic Wisdom-Holman integrator for
  long-term gravitational simulations. *MNRAS* 452, 376.
  doi:10.1093/mnras/stv1257
  and Wisdom, J. & Holman, M. (1991). Symplectic maps for the N-body problem.
  *AJ* 102, 1528. doi:10.1086/115978
- **High-order kernel and correctors (WHCKL, used in runs without MEGNO):**
  Rein, H., Tamayo, D. & Brown, G. (2019). High-order symplectic integrators
  for planetary dynamics and their implementation in REBOUND. *MNRAS* 489,
  4632. doi:10.1093/mnras/stz2503
- **Variational equations (used for MEGNO and the Lyapunov number):**
  Rein, H. & Tamayo, D. (2016). Second-order variational equations for N-body
  simulations. *MNRAS* 459, 2275. doi:10.1093/mnras/stw644
- **SimulationArchive:** Rein, H. & Tamayo, D. (2017). A new paradigm for
  reproducing and analyzing N-body simulations of planetary systems.
  *MNRAS* 467, 2377. doi:10.1093/mnras/stx232
- **MEGNO indicator** (implemented in REBOUND): Cincotta, P. M. & Simó, C.
  (2000). Simple tools to study global dynamics in non-axisymmetric galactic
  potentials – I. *A&AS* 147, 205. doi:10.1051/aas:2000108
  and Cincotta, P. M., Giordano, C. M. & Simó, C. (2003). Phase space
  structure of multi-dimensional systems by means of the mean exponential
  growth factor of nearby orbits. *Physica D* 182, 151.
  doi:10.1016/S0167-2789(03)00103-9

## REBOUNDx

A library of additional physical effects for REBOUND by Daniel Tamayo and
collaborators.

- Repository: <https://github.com/dtamayo/reboundx>
- Documentation: <https://reboundx.readthedocs.io>
- License: GNU General Public License v3.0
- Version used for the included results: 4.6.2

**Citation:** Tamayo, D., Rein, H., Shi, P. & Hernandez, D. M. (2020).
REBOUNDx: a library for adding conservative and dissipative forces to
otherwise symplectic N-body integrations. *MNRAS* 491, 2885.
doi:10.1093/mnras/stz2870

Effects used:

- **`gr_potential`** (on by default): Nobili, A. & Roxburgh, I. W. (1986).
  Simulation of general relativistic corrections in long term numerical
  integrations of planetary orbits. *IAU Symposium* 114, 105.
- **`gravitational_harmonics`** (optional solar J₂, off by default). The default
  value J₂ = 2.2 × 10⁻⁷ follows Mecheri, R. et al. (2004), *Solar Physics*
  222, 191.

## Planetary data

- **Initial orbital elements (default, built-in table):** J2000 heliocentric
  ecliptic elements from E. M. Standish and J. G. Williams, *Keplerian
  Elements for Approximate Positions of the Major Planets* (1800 AD–2050 AD
  table). Published by the JPL Solar System Dynamics group at
  <https://ssd.jpl.nasa.gov/planets/approx_pos.html>, and in Standish & Williams,
  "Orbital Ephemerides of the Sun, Moon, and Planets", *Explanatory Supplement
  to the Astronomical Almanac* (3rd ed., Urban & Seidelmann, eds., 2013).
  All included results were generated from this table.
- **Planetary masses:** reciprocal masses from the JPL DE405 ephemeris
  (Standish, E. M. 1998, JPL IOM 312.F-98-048), with Earth taken as the
  Earth–Moon barycentre.
- **Reference secular frequencies** (gᵢ, sᵢ, used to annotate spectra and for
  comparison): Laskar, J. et al. (2004). A long-term numerical solution for
  the insolation quantities of the Earth. *A&A* 428, 261.
  doi:10.1051/0004-6361:20041335
- **NASA/JPL Horizons (optional):** `ceam-sim horizons` can fetch initial
  conditions for any epoch from the JPL Horizons system through REBOUND's
  `sim.add(name, date=...)` interface. None of the included results use it.
  Service: <https://ssd.jpl.nasa.gov/horizons/>.
  Reference: Giorgini, J. D. et al. (1996). JPL's On-Line Solar System Data
  Service. *BAAS* 28, 1158.

## Python dependencies

| Package | Use | License |
|---|---|---|
| NumPy | arrays, FFTs, fitting | BSD-3-Clause |
| SciPy (optional) | spectral peak finding | BSD-3-Clause |
| Matplotlib | figures | Matplotlib License (PSF-based) |
| pandas | statistics tables | BSD-3-Clause |
| PyYAML | configuration files | MIT |
| tqdm | progress bars | MPL-2.0 AND MIT |

## Scientific literature

The scientific framing, target values and methodological choices follow the
literature on Solar System chaos. The following works most directly shaped
this project:

- Laskar, J. (1989). *Nature* 338, 237 — discovery of inner-Solar-System
  chaos and the ~5 Myr Lyapunov time this project reproduces.
- Laskar, J. (2008). Chaotic diffusion in the Solar System. *Icarus* 196, 1 —
  the g₁ − g₅ resonance and the role of GR.
- Batygin, K. & Laughlin, G. (2008). *ApJ* 683, 1207.
- Ito, T. & Tanikawa, K. (2002). *MNRAS* 336, 483.
- Laskar, J. & Robutel, P. (2001). High order symplectic integrators for
  perturbed Hamiltonian systems. *Celest. Mech. Dyn. Astron.* 80, 39.
- Mogavero, F. & Laskar, J. (2022). *A&A* 662, L3.
- Hoang, N. H., Mogavero, F. & Laskar, J. (2022). *MNRAS* 514, 1342.
- Mogavero, F., Hoang, N. H. & Laskar, J. (2023). *Phys. Rev. X* 13, 021018.
- Brown, G. & Rein, H. (2020). A repository of vanilla long term integrations
  of the Solar System. *RNAAS* 4, 221. arXiv:2012.05177 — the methodological
  template (WHFast + `gr_potential` in REBOUND).
- Brown, G. & Rein, H. (2023). General relativistic precession and the
  long-term stability of the Solar System. *MNRAS* 521, 4349.
  arXiv:2303.05567.
- Abbot, D. S. et al. (2023). Mercury's chaotic secular evolution as a
  subdiffusive process. arXiv:2306.11870.
- Laskar, J. (1997). Large scale chaos and the spacing of the inner planets.
  *A&A* 317, L75 — angular-momentum deficit.

Full references are listed in [`README.md`](README.md#references). These papers
are cited, not redistributed. Copies are not included in this repository.

## Development tools

Code and documentation were developed with assistance from Anthropic's
Claude (Claude Code).

## License compatibility

REBOUND and REBOUNDx are GPL-3.0 licensed, and this project imports both. This
project is therefore also released under the GNU General Public License,
version 3 or (at your option) any later version. See [`LICENSE`](LICENSE).
