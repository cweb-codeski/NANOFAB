# Equations Reference

Every equation used in `src/nanofab/`, with variables, units, assumptions and source. SI units
throughout. The code docstrings repeat these. If you change one, change both.

## Notation

| Symbol | Meaning | Unit |
| --- | --- | --- |
| d, l, h | pillar diameter, pitch, height | m |
| φ = πd²/(4l²) | areal solid fraction | – |
| γ | surface tension | N/m |
| θ_r | receding contact angle | rad |
| μ, ρ | liquid viscosity, density | Pa·s, kg/m³ |
| h_lv | latent heat of vaporization | J/kg |
| K | permeability | m² |
| L | feed length (feed to farthest evaporating point) | m |
| q'' | heat flux per footprint area | W/m² |
| s | accommodation coefficient | – |
| M, R | molar mass, universal gas constant | kg/mol, J/(mol·K) |

## Capillary pressure (`capillary.py`)

Young-Laplace:

$$\Delta P = \gamma\left(\frac{1}{R_1}+\frac{1}{R_2}\right)$$

Pillar array maximum capillary pressure (force balance at the receding angle; Xiao, Enright &
Wang, Langmuir 2010):

$$P_{c,max} = \frac{\gamma\cos\theta_r\,\pi d}{l^2-\pi d^2/4}$$

Assumes a static meniscus, uniform wettability, and pillars tall enough that the meniscus depins
before touching the floor.

## Permeability (`permeability.py`)

Darcy: \(U = -(K/\mu)\,dP/dx\), with U the superficial velocity over the full cross-section.

Infinite square array, transverse flow (Sangani & Acrivos 1982, in the Drummond & Tahir 1984
form), with \(a=d/2\):

$$K_s = \frac{a^2}{8\phi}\left(-\ln\phi-1.476+2\phi-1.774\phi^2+4.076\phi^3\right)$$

Finite height with a no-slip floor and a shear-free flat top (Brinkman 1949, μ_eff = μ; the
uniform-K Brinkman approach for pillar wicks is attributed to Adera et al. 2016 in V2019 Sec. 1;
the closed form below is our derivation, not checked against a paper):

$$K_{eff} = K_s\left[1-\frac{\tanh(h/\sqrt{K_s})}{h/\sqrt{K_s}}\right]$$

Limits: \(K_{eff}\to K_s\) for tall pillars; \(K_{eff}\to h^2/3\) (shear-free film) for very
short or sparse pillars. Assumptions: Stokes flow, flat top, φ below about 0.4.

Extraction from simulation: \(K=\mu (Q/A) L/\Delta P\).

## Capillary-limited dryout (`dryout.py`)

The wick is fed at x = 0 and closed at x = L, with uniform q'' fully converted to evaporation.

$$P(x)=P_{feed}-\frac{\mu q''}{\rho K h\,h_{lv}}\left(Lx-\frac{x^2}{2}\right),\qquad
\Delta P=\frac{\mu q'' L^2}{2\rho K h\,h_{lv}}$$

$$q''_{dry}=\frac{2\rho K h\,h_{lv}\,\Delta P_{avail}}{\mu L^2},\qquad \Delta P_{avail}=P_{c,max}\ (\text{minus gravity, feed losses})$$

Key scaling: \(q''_{dry}\propto K P_c/L^2\). The \(K P_c\) product is the figure of merit for a wick
geometry; \(L\) is set by the architecture (feed spacing).

## Kinetic evaporation (`evaporation.py`)

Schrage (Carey, *Liquid-Vapor Phase-Change Phenomena*, Ch. 16):

$$\dot m''=\frac{2s}{2-s}\sqrt{\frac{M}{2\pi R}}\left(\frac{p_{sat}(T_{lv})}{\sqrt{T_{lv}}}-\frac{p_v}{\sqrt{T_v}}\right),\qquad q''=\dot m'' h_{lv}$$

Linearized interfacial coefficient (saturated vapor, small superheat):

$$h_i=h_{lv}\frac{2s}{2-s}\sqrt{\frac{M}{2\pi R}}\frac{1}{\sqrt T}\left(\frac{dp_{sat}}{dT}-\frac{p_{sat}}{2T}\right)$$

Assumptions: no disjoining-pressure or curvature correction to p_sat, and near-equilibrium vapor
(no Knudsen-layer resolution).

Reduced linearized Boltzmann transport relation (form of the V2019 interface model, Eqs. 5-7;
Sone 2000), with the vapor pressure at the Knudsen-layer edge set to the far-field value. Here
\(R_s=R/M\). This is **not** a paper-verified implementation: the temperature-jump relation
(V2019 Eq. 6) is not used, \(\rho_k\approx p_{sat}(T_l)/(R_sT_l)\), and the nonunity-accommodation
correction below was written from memory; V2019 attributes its treatment to Meland & Ytrehus
(2003, DOI 10.1063/1.1564097), which is not held (open source gap):

$$\frac{p_k}{p_{sat}(T_l)}-1=-C_4\frac{u_k}{\sqrt{2R_sT_l}},\qquad q''_e=\rho_k u_k h_{lv}$$

$$\dot m''=\frac{\left(p_{sat}(T_l)-p_k\right)}{C_4(a)}\sqrt{\frac{2}{R_sT_l}},\qquad
C_4(a)=2.13204+2\sqrt{\pi}\,\frac{1-a}{a}$$

At a = 0.052 this is within ~3-5% of the Schrage flux (water, 100 °C): a cross-model check
between two approximate relations, not a verification of either.

## Conduction (`conduction.py`)

- Slab: \(R''=t/k\).
- Pillar wick bounds: parallel \(k=\phi k_s+(1-\phi)k_l\) (upper bound), series
  \(k=[\phi/k_s+(1-\phi)/k_l]^{-1}\) (lower bound).
- Spreading resistance: Lee, Song, Au & Moran (1995), average-temperature form. See the docstring.

## Dimensionless groups (`dimensionless.py`)

| Group | Definition | Small means |
| --- | --- | --- |
| Re | ρUL/μ | Stokes flow valid |
| Ca | μU/γ | meniscus shape set by statics |
| Pe | UL/α | liquid advection negligible in heat equation |
| Bo | Δρ g L²/γ | gravity negligible |
| Kn | λ/L, λ = k_B T/(√2 π d_m² p) | continuum gas valid |

## Unit-cell metrics (`analysis/heat_transfer_from_comsol.py`)

$$q''=\frac{Q_{in}}{A_{fp}},\quad h_{evap}=\frac{q''}{T_{base}-T_v},\quad
R''_{wick,app}=\frac{T_{base}-T_{lv}}{q''},\quad R''_{int,app}=\frac{T_{lv}-T_v}{q''}$$

The h_evap definition (which temperature, which area) must be matched to the reference paper
before comparing values. \(R''_{wick,app}\) and \(R''_{int,app}\) are apparent splits at the
area-averaged meniscus temperature \(T_{lv}\); with nonuniform interface temperature and flux they
depend on the averaging convention and are not unique resistances.

## Passive device model (`src/nanofab/device/`)

Half device \(0\le x\le L\), substrate \(-t_s\le z\le 0\), \(P_c=p_v-p_{liquid}\) (positive suction),
\(Q\) total flow across \(W\), transmissivity \(\tau(P_c)=hK_{bulk}(P_c)\):

$$\frac{Q}{W}=\frac{\tau(P_c)}{\mu}\frac{dP_c}{dx},\qquad \frac{dQ}{dx}=-\frac{W q''_e}{\rho h_{fg}},\qquad
P_c(0)=P_{c,feed},\ Q(L)=0$$

$$\nabla\cdot(k_s\nabla T)=0,\quad -k_s\partial_zT|_{-t_s}=q''_{in},\quad -k_s\partial_zT|_{0}=q''_e=h_e(P_c)\,(T_e-T_v),\quad
\partial_xT|_{0,L}=0$$

Kirchhoff transform \(\Psi(P_c)=\int\tau\,dP_c\) and the evaporation-moment identity (PROJECT_DERIVED):

$$\int_{P_c(0)}^{P_c(L)}\frac{\tau}{\mu}\,dP_c=\int_0^L\frac{x\,q''_e(x)}{\rho h_{fg}}\,dx$$

Full derivation, numerics, statuses and verification gates: `docs/passive_device_model.md`.

## Electroosmotic assistance (`src/nanofab/electroosmosis/`)

Thin-EDL Helmholtz-Smoluchowski mobility and superficial EO velocity (uniform zeta, insulating
solid, no surface conduction; F = formation factor, geometric):

$$\mu_{EO}=-\frac{\varepsilon_r\varepsilon_0\zeta}{\mu},\qquad U=\frac{K}{\mu}\left(-\frac{dp}{dx}\right)+\frac{\mu_{EO}E}{F},\qquad F_{Maxwell}=\frac{1+\phi}{1-\phi}$$

K is the same intrinsic permeability in the passive and assisted cases. Uniform-flux half-wick
capillary limit with EO, \(\Pi=\mu U_{EO}L/(KP_{c,max})=V\varepsilon|\zeta|/(FKP_{c,max})\):

$$\frac{q''_{dry,EO}}{q''_{dry,0}}=\frac{1+\Pi+\sqrt{1+2\Pi}}{2},\qquad \text{over-pumping (}P_c(L)<0\text{) when } U_0<2U_{EO}$$

Derivations, assumptions and validity diagnostics: `docs/electroosmosis_model.md`.
