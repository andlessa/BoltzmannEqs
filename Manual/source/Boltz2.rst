=========================================
Boltzmann Equations for the PQMSSM (v2.3)
=========================================

:Author: Andre Lessa

.. role:: math(raw)
   :format: html latex
..

General Formalism and Approximations
====================================

The general Boltzmann equation for the number distribution of a particle
species can be written as (assuming isotropy):

.. math:: {\frac{\partial F_{i}}{\partial t}} -H p {\frac{\partial F_{i}}{\partial p}} = C_{i}[F_{i},F_{j},p] \label{eq:d1}

where :math:`F_{i}(p)` is the number distribution of particle :math:`i`
as function of momentum :math:`p`, :math:`C` represents a source/sink
term and :math:`H` is the Hubble constant:

.. math:: H = \sqrt{\frac{8 \pi}{3} \frac{\rho_T}{M_P^2}} \label{H}

with :math:`M_{P} = 1.22\times 10^{19}` GeV and
:math:`\rho_T = \sum_{i} \rho_i`. The number, energy and pressure
densities are given in terms of :math:`F_{i}` as:

.. math::

   \begin{aligned}
   n_{i}(t) & = & \int \frac{dp}{2 \pi^2} p^2 F_i(p) \nonumber \\ 
   \rho_{i}(t) & = & \int \frac{dp}{2 \pi^2} p^2 E_i F_i(p) \label{beqs}\\
   P_{i}(t) & = & \frac{1}{3} \int \frac{dp}{2 \pi^2} \frac{p^4}{E_i} F_i(p) \nonumber\end{aligned}

where :math:`m_i` is the mass of particle :math:`i` and
:math:`E_i = \sqrt{p_i^2 + m_i^2}`. Using Eq.([eq:d1]), we obtain the
following equations for the number and energy densities:

.. math::

   \begin{aligned}
   {\frac{d n_i}{d t}} + 3H n_i & = & \int \frac{dp}{2 \pi^2} p^2 C_i \nonumber \\
   {\frac{d \rho_i}{d t}} + 3H (\rho_i + P_i) & = & \int \frac{dp}{2 \pi^2} p^2 E_i C_i \label{eq:meqs}\end{aligned}

The collision term, :math:`C_i`, for the process
:math:`i + f + \ldots \leftrightarrow a
+ b + c + \ldots` is given by:

.. math::

   \begin{aligned}
   C_i & = & \frac{1}{E_i} \int \prod_{j,a} \frac{d^3 p_j}{2 E_j (2 \pi)^3}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i + p_j + \ldots - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   &\times& \left[(1 \pm f_a) (1 \pm
   f_b)\ldots f_i f_j\ldots - f_a f_b \ldots (1 \pm f_i)(1 \pm f_j)\ldots \right]\end{aligned}

where the plus (minus) sign is for bosons (fermions). Below we always
assume :math:`f_{i,j,a,..} \ll 1`, so:

.. math::

   \begin{aligned}
   C_i & \simeq & \frac{1}{E_i} \int \prod_{j,a} \frac{d^3 p_j}{2 E_j (2 \pi)^3}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i + p_j + \ldots - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   &\times& \left[f_i f_j\ldots - f_a f_b \ldots \right]\end{aligned}

We will assume that :math:`C` is given by:

.. math:: C = C_{dec} + C_{prod} + C_{ann}

where :math:`C_{dec}` contains the contributions from decays and inverse
decays (:math:`i \leftrightarrow a + b + \ldots`), :math:`C_{prod}`
contains the contributions from decay injection and inverse decay
injection (:math:`a \leftrightarrow i + b + \ldots`) and :math:`C_{ann}`
from annihilations with the thermal plasma
(:math:`i + i \leftrightarrow a + b`). Below we compute each term
separately, under some assumptions.

Annihilation Term
-----------------

The annihilation term :math:`C_{ann}` for the
:math:`i + j \leftrightarrow a + b` process is given by:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{ann} = \int d\Pi_{i} d\Pi_{j} d\Pi_{a}
   d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \left[ f_a f_b -
   f_i f_j \right]

where :math:`d\Pi_{i} = d^{3} p_i/((2\pi)^3 2 E_i)`. Since we are
ultimately interested in Eqs.([eq:meqs]) for the number and energy
densities, we will consider the following integral:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{ann}  E_i^{\alpha} = \int d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 
   \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2
    \left[ f_a f_b - f_i f_j \right] E_i^{\alpha}

where :math:`\alpha = 0 (1)` for the number (energy) density. Here we
assume that the distributions can be approximated by [1]_:

.. math:: f_i \simeq \exp(-(E_i - \mu_i)/T)

so the annihilation term can then be written as:

.. math::

   \begin{aligned}
   & \int & \frac{dp}{2 \pi^2} p^2 C_{ann}  E_i^{\alpha} =  -\left( \exp((\mu_i + \mu_j)/T) -\exp((\mu_a + \mu_b)/T)\right) \nonumber \\
    & \times & \int  d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha} \nonumber\end{aligned}

where above we have used conservation of energy
(:math:`E_i + E_j = E_a + E_b`). Since for the cases of interest the
equilibrium distributions have zero chemical potential, we have:

.. math:: \frac{n_i}{\bar{n}_i} = \exp(\mu_i/T)

so:

.. math::

   \begin{aligned}
   & \int & \frac{dp}{2 \pi^2} p^2 C_{ann} E_i^{\alpha} = -\left( \frac{n_i n_j}{\bar{n}_i \bar{n}_j} - \frac{n_a n_b}{\bar{n}_a \bar{n}_b}\right) \nonumber \\
    & \times & \int  d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha} \nonumber\end{aligned}

In particular, for the process :math:`i + i \leftrightarrow a + b`,
where :math:`a` and :math:`b` are in thermal equilibrium
(:math:`\mu_a = \mu_b = 0`):

.. math::

   \begin{aligned}
   & \int & \frac{dp}{2 \pi^2} p^2 C_{ann} E_i^{\alpha} =  -\left( \frac{n_i^2}{\bar{n}_i^2} - 1 \right) \nonumber \\
   &  \times & \int d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha}  \nonumber \\
    & = & -\left( n_i^2 - \bar{n}_i^2 \right) \langle \sigma v E_i^{\alpha} \rangle\end{aligned}

For :math:`\alpha = 0`, the above equation is the well known
contribution from thermal scatterings to the annihilation term. To
estimate its value for :math:`\alpha = 1`, we assume:

.. math:: \langle \sigma v E \rangle \simeq \langle \sigma v \rangle \langle E_i \rangle = \langle \sigma v \rangle \frac{\rho_i}{n_i} \label{eq:app}

where :math:`\langle \;\; \rangle` represents thermal average. Thus:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{ann} E_i^{\alpha}  = \left( \bar{n}_i^2 - n_i^2 \right) \left\{ \begin{array}{rl}  
   \langle \sigma v \rangle & \mbox{, for $\alpha = 0$} \\
   \langle \sigma v \rangle \frac{\rho_i}{n_i} &\mbox{, for $\alpha = 1$}
   \end{array} \right. \label{eq:collfin}

Decay Term
----------

Now we derive a simplified expression for the decay (and inverse decay)
term, under approximations similar to the ones used in the last section.
The decay term includes the contributions from particle decay and
inverse decay:

.. math::

   C_{dec} \simeq \frac{1}{E_i} \int \prod_{a} \frac{d^3 p_a}{2 E_a (2 \pi)^3}
   (2 \pi)^4 \delta^{4}\left(p_i - p_a - p_b \ldots\right) |\mathcal{M}|^2 \left[f_i - f_a f_b \ldots \right]
   \label{eq:dec0}

As in the case of the annihilation term, we assume that the
distributions for :math:`a,b,\ldots` can be approximated by
:math:`f_x \simeq \exp(-(E_x -
\mu_x)/T)`, so we can write:

.. math::

   f_a f_b \ldots \simeq \exp\left(\frac{\mu_a +
   \mu_b + \ldots}{T}\right) \exp(-E_i/T) = \frac{n_a n_b \ldots}{\bar{n}_a
   \bar{n}_b \ldots} \exp(-E_i/T)  =  \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots}
   \bar{f}_{i}

where we used conservation of energy (:math:`E_a + E_b + \ldots = E_i`)
and :math:`\bar{f}_i` is the equilibrium distribution for the species
:math:`i`. Hence we can write Eq.([eq:dec0]) as:

.. math::

   \begin{aligned}
   C_{dec} & \simeq & \left[f_i - \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots}
   \bar{f}_{i} \right] \frac{1}{E_i} \int \prod_{a}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   & = & \mathcal{B}_{ab\ldots} \frac{\Gamma_i m_i}{E_i} \left[f_i -
   \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \bar{f}_{i} \right] \end{aligned}

where :math:`\Gamma_i` is the width for :math:`i` and
:math:`\mathcal{B}_{ab\ldots} \equiv BR(i \to a + b + \ldots)`

Once again we consider the integral:

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{dec}(p) E_i^{\alpha} = 
    & - & \Gamma_i \int \frac{dp}{2 \pi^2} p^2 \frac{m_i}{E_i} f_i E_i^{\alpha}
    \nonumber \\
    & + & \sum_{i \; decays} \mathcal{B}_{ab\ldots}
   \Gamma_i \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \int \frac{dp}{2 \pi^2}
   p^2 \frac{m_i}{E_i} \bar{f}_{i} E_i^\alpha \label{eq:dec2}\end{aligned}

where we have included the sum over all decay channels and
:math:`\alpha = 0 (1)` for the contribution to the number (energy)
density equation. Note that both integrals are identical, except for the
replacement :math:`f_i \to \bar{f_i}`. The first integral in
Eq.([eq:dec2]) gives:

.. math::

   -\Gamma_i \int \frac{dp}{2 \pi^2} p^2 \frac{m_i}{E_i} f_i(p) E_i^{\alpha} =
   \left\{ \begin{array}{rl} -\Gamma_i m_i n_i \langle \frac{1}{E_i} \rangle  & \mbox{, for $\alpha = 0$} \\
   -\Gamma_i m_i n_i &\mbox{, for $\alpha = 1$}
   \end{array} \right. \label{eq:dec1a}

where

.. math::

   \langle \frac{1}{E_i} \rangle \equiv \frac{1}{n_i} \int \frac{dp}{2 \pi^2} p^2
   \frac{1}{E_i} f_i(p)

Hence we can write Eq.([eq:dec2]) as:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{dec}(p) E_i^{\alpha} = -\Gamma_i m_i 
   \left\{ \begin{array}{ll} n_i \langle \frac{1}{E_i} \rangle - \bar{n}_i  \langle
   \frac{1}{E_i}
   \rangle_{eq} \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b\ldots}{\bar{n}_a \bar{n}_b\ldots}  & \mbox{, for $\alpha = 0$}  \\
    n_i - \bar{n}_i \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b\ldots}{\bar{n}_a \bar{n}_b\ldots}  & \mbox{, for $\alpha = 1$}
   \end{array} \right. \label{eq:decfin}

For the non-equilibrium average we assume:

.. math::

   \langle \frac{1}{E_i} \rangle \simeq \frac{1}{\langle E_i \rangle} =
   \frac{n_i}{\rho_i}

which is exact in the non-relativistic limit, but it is only an
approximation for the relativistic case. Although we can compute the
equilibrium average (:math:`\langle
\frac{1}{E_i}\rangle_{eq}`) explicitly, in order to have an exact
cancellation between the decay and inverse decay terms when :math:`i`,
:math:`a` and :math:`b` are all in equilibrium, we take:

.. math::

   \langle \frac{1}{E_i} \rangle_{eq} \simeq \langle \frac{1}{E_i} \rangle =
   \frac{n_i}{\rho_i}

With the above approximations we finally obtain:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{dec}(p) E_i^{\alpha} = 
    -\Gamma_i  m_i \left\{ \begin{array}{ll}\frac{n_i}{\rho_i}\left( n_i -
    \bar{n}_i \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \right)   &
    \mbox{, for $\alpha = 0$}
    \\
    n_i - \bar{n}_i \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots}  & \mbox{, for $\alpha = 1$}
   \end{array} \right. \label{eq:decfin}

where :math:`\mathcal{B}_{ab\ldots} \equiv BR(i\to a+b+\ldots)`.

Production Term
---------------

The decay and inverse decay of other particles
(:math:`a \to i + b + \ldots`) can also affect the species :math:`i`.
The contribution from these terms we label :math:`C_{prod}`, which is
given by:

.. math::

   C_{prod} \simeq \frac{1}{E_i} \int \frac{d^3 p_a}{2 E_a (2
   \pi)^3} \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_a
   - p_i - p_b \ldots\right) |\mathcal{M}|^2 \left[f_a - f_i f_b \ldots \right]

Using the same approximations of the previous section, we write:

.. math::

   f_i f_b\ldots \simeq  \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots}
   e^{-E_a/T} = \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots}
   \bar{f}_{a}

Hence:

.. math::

   C_{prod} = \frac{1}{E_i} \int \frac{d^3 p_a}{2 E_a (2 \pi)^3} \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} 
   (2 \pi)^4 \delta^{4}\left(p_a - p_i - p_b \ldots\right) |\mathcal{M}|^2
   \left(f_a - \bar{f}_a \frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right)

and

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) E_i^\alpha & = & 
   \int \frac{d^3 p_a}{E_a (2 \pi)^3} \left(f_a - \bar{f}_a \frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right) \nonumber \\
   & \times & \frac{d^3 p E_i^{\alpha}}{2 E_i (2 \pi)^3}
   \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_a - p_i - p_b \ldots\right) |\mathcal{M}|^2
   \label{eq:prod2}\end{aligned}

with :math:`\alpha = 0 (1)` for the contribution to the number (energy)
density equation. For :math:`\alpha = 0` we obtain:

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) & = & \Gamma_a  \mathcal{B}_{i} m_a 
   \int \frac{d^3 p_a}{E_a (2 \pi)^3} \left(f_a - \bar{f}_a \sum_b
   \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}}\frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & = & \Gamma_a \mathcal{B}_{i} m_a \frac{n_a}{\rho_a} \left( n_a - \bar{n}_a
     \sum_b \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i n_b
     \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)\end{aligned}

where :math:`\mathcal{B}_{ib\ldots} \equiv BR(a \to i + b + \ldots)`,
:math:`\mathcal{B}_i
= \sum_{b} \mathcal{B}_{ib\ldots}` and we have once again assumed
:math:`\langle 1/E_a
\rangle \simeq \langle 1/E_a \rangle_{eq} \simeq n_a/\rho_a`.

For :math:`\alpha = 1`, the integral in Eq.([eq:prod2]) does not take a
simple form. In order to compute it, we assume:

.. math:: E_i \simeq \frac{E_a}{2}

The above expression is only exact for 2-body decays and :math:`m_a \gg
m_i,m_b`. For the remaining cases, it is only an estimate.

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) E_i & \simeq & 
   \Gamma_a \mathcal{B}_{i}  \frac{m_a}{2} \int \frac{d^3 p_a}{(2
   \pi)^3} \left(f_a - \bar{f}_a \sum_b
   \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}}
    \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & = & \Gamma_a \mathcal{B}_{i}  \frac{m_a}{2} \left( n_a -
   \bar{n}_a \sum_b \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i n_b
   \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)\end{aligned}

Combining the results for :math:`\alpha = 0` and 1, we have:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) E_i^{\alpha} = 
   \Gamma_a \mathcal{B}_{i} m_a  \left( n_a - \bar{n}_a
   \sum_b \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i n_b
   \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right) \left\{ \begin{array}{ll}  \frac{n_a}{\rho_a}  & \mbox{, for $\alpha = 0$} 
   \\
    \frac{1}{2}  & \mbox{, for $\alpha = 1$}
   \end{array} \right. \label{eq:prodfin}

Number and Energy Density Equations
-----------------------------------

Using the results of Eqs.([eq:collfin]), ([eq:decfin]) and
([eq:prodfin]) in the Boltzmann equations for :math:`n_i` and
:math:`\rho_i` (Eq.([eq:meqs])), we obtain:

.. math::

   \begin{aligned}
   {\frac{d n_i}{d t}} + 3H n_i  & = &  \left( \bar{n}_i^2 - n_i^2 \right) \langle \sigma
   v \rangle - \Gamma_i m_i \frac{n_i}{\rho_i}\left(n_i - \bar{n}_i \sum_{i\to\ldots}
   \mathcal{B}_{ab\ldots} \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & + & \sum_a 
   \Gamma_a \mathcal{B}_i m_a \frac{n_a}{\rho_a} \left(n_a - \bar{n}_a \sum_{a \to
   i\ldots} \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)  + C_{i}(T) \label{eq:nieq} \\
   {\frac{d \rho_i}{d t}} + 3H (\rho_i + P_i) & = & \left( \bar{n}_i^2 - n_i^2 \right)
   \langle \sigma v \rangle \frac{\rho_i}{n_i} - \Gamma_i m_i \left( n_i -
   \bar{n}_i \sum_{i\to\ldots} \mathcal{B}_{ab\ldots} \frac{n_a n_b\ldots}{\bar{n}_a
   \bar{n}_b\ldots}\right) \nonumber \\
    & + & \sum_a \Gamma_a  \mathcal{B}_i \frac{m_a}{2} \left( n_a -
    \bar{n}_a \sum_{a \to i\ldots}  \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i
    n_b..}{\bar{n}_i \bar{n}_b..} \right) + \tilde{C}_{i}(T)
    \frac{\rho_i}{n_i}\end{aligned}

where :math:`\mathcal{B}_{ab\ldots} = BR(i \to a + b+ \ldots)`,
:math:`\mathcal{B}_{ib\ldots} =
BR(a \to i + b + \ldots)`,
:math:`\mathcal{B}_i = \sum_b \mathcal{B}_{ib\ldots}` and we have
included an extra term (:math:`C_i` and :math:`\tilde{C}_i`) to allow
for other possible sources for the number and energy densities. For
simplicity we assume :math:`C_i = \tilde{C}_{i}` from now on.

It is also convenient to use the above results to obtain a simpler
equation for :math:`\rho_i/n_i`:

.. math::

   {\frac{d \rho_i/n_i}{d t}} \equiv {\frac{d R_i}{d t}} = -3 H \frac{P_i}{n_i} + \sum_{a}
   \mathcal{B}_{i} \frac{\Gamma_a m_a}{n_i} \left( \frac{1}{2} - \frac{n_a}{\rho_a} \frac{\rho_i}{n_i} \right) \left(n_a -
   \bar{n}_a \sum_{a \to i\ldots} \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i
    n_b..}{\bar{n}_i \bar{n}_b..}\right) \label{eq:Rieq}

Besides the above equations, it is useful to consider the evolution
equation for entropy:

.. math:: dS \equiv \frac{dQ^{dec}}{T}

where :math:`dQ^{dec}` is the net energy injected from decays. With the
above definition we have:

.. math::

   \begin{aligned}
   \dot{S} & = & \frac{1}{T}\sum_i BR(i,X)
   \frac{d\left(R^3 \rho_i\right)^{dec}}{dt}  \nonumber \\
   \To \dot{S} & = & \frac{R^3}{T}\sum_i BR(i,X)
   \Gamma_i m_i\left(n_i - \bar{n}_i \sum_{i\to\ldots} \mathcal{B}_{ab\ldots} \frac{n_a n_b\ldots}{\bar{n}_a
   \bar{n}_b\ldots} \right) \label{Seq}\end{aligned}

where :math:`R` is the scale factor and :math:`BR(i,X)` is the fraction
of energy injected in the thermal bath from :math:`i` decays.

The above expressions can be written in a more compact form if we define
the following ”effective thermal densities” and ”effective BR”:

.. math::

   \begin{aligned}
   \mathcal{N}^{th}_{X} & \equiv &  \bar{n}_X \sum_{X \to \ldots} BR(X \to 1 + 2 +
   \ldots)
   \prod_{k}
   \frac{n_k}{\bar{n}_k} \nonumber \\
   \mathcal{N}^{th}_{XY} & \equiv & \frac{\bar{n}_X}{\mathcal{B}^{eff}_{XY}}
   \sum_{X \to Y + \ldots} g_Y BR(X \to g_Y Y + 1 + \ldots)
   \left(\frac{n_Y}{\bar{n}_Y}\right)^{g_Y} \prod_{k} \frac{n_k}{\bar{n}_k}
   \nonumber \\
   \mathcal{B}^{eff}_{XY} & \equiv & \sum_{X \to Y + \ldots} g_Y BR(X \to g_Y Y +
   1+\ldots) \nonumber\end{aligned}

where :math:`g_Y` is the :math:`Y` multiplicity in the final state of
:math:`X` decays. In addition, defining:

.. math:: x = \ln(R/R_0),\;\; N_i = \ln(n_i/s_0),\;\; {\rm and}\;\; N_S = \ln(S/S_0)

we can write Eqs.([Seq]), ([eq:nieq]) and ([eq:Rieq]) as:

.. math::

   \begin{aligned}
   N_S' & = & \frac{e^{(3 x - N_S)}}{HT} \sum_{i} BR(i,X) \Gamma_i m_i \left(n_i -
   \mathcal{N}_{i}^{th} \right) 
   \label{Seqb} \\
   N_i' & = & -3 + \frac{\sigv_i}{H} n_i [\left(\frac{\bar{n}_i}{n_i}\right)^2
   -1] -  \frac{\Gamma_i}{H} \frac{m_i}{R_i}\left(1 -
   \frac{\mathcal{N}_{i}^{th}}{n_i} \right) \nonumber  \\
    & + & \sum_{a} \mathcal{B}_{ai}^{eff} \frac{\Gamma_a}{H}
    \frac{m_a}{R_a}\left(\frac{n_a}{n_i} - \frac{\mathcal{N}_{ai}^{th}}{n_i}
     \right)
    \\
   R_i' & = &  -3 \frac{P_i}{n_i} + \sum_{a} \mathcal{B}_{ai}^{eff}
   \frac{\Gamma_a}{H} m_a \left( \frac{1}{2} - \frac{R_i}{R_a} \right) \left(\frac{n_a}{n_i} -
   \frac{\mathcal{N}_{ai}^{th}}{n_i} \right)
   \label{Nieq}\end{aligned}

where :math:`'=d/dx`.

The above equation for :math:`N_i` also applies for coherent oscillating
fields, if we define:

.. math:: N_i = \ln(n_i/s_0),\;\; {\rm and}\;\; n_i \equiv \rho_i/m_i

so

.. math::

   \begin{aligned}
   N_i' & = & -3 - \frac{\Gamma_i}{H}  \nonumber \\
   R_i'& = & 0 \label{Nico}\end{aligned}

where we assume that the coherent oscillating component does not couple
to any of the other fields.

Collecting Eqs.([Seqb])-([Nieq]) and ([Nico]) we have a closed set of
first order differential equations:

-  Entropy:

   .. math::

      N_S' = \frac{e^{(3 x - N_S)}}{HT} \sum_{i} BR(i,X) \Gamma_i m_i \left(n_i -
      \mathcal{N}_{i}^{th} \right) \label{eq:Sfin}

-  Thermal fields:

   .. math::

      \begin{aligned}
      N_i'& = & -3 + \frac{\sigv_i}{H} n_i [\left(\frac{\bar{n}_i}{n_i}\right)^2
      -1] -  \frac{\Gamma_i}{H} \frac{m_i}{R_i}\left(1 - \frac{\mathcal{N}_{i}^{th}}{n_i}
       \right)
        +  \sum_{a} \mathcal{B}_{ai}^{eff} \frac{\Gamma_a}{H}
       \frac{m_a}{R_a}\left(\frac{n_a}{n_i} - \frac{\mathcal{N}_{ai}^{th}}{n_i}
        \right) \nonumber
       \\
      R_i' & = &  -3 \frac{P_i}{n_i} + \sum_{a} \mathcal{B}_{ai}^{eff}
      \frac{\Gamma_a}{H} m_a \left( \frac{1}{2} - \frac{R_i}{R_a} \right) \left(\frac{n_a}{n_i} -
      \frac{\mathcal{N}_{ai}^{th}}{n_i} \right)\end{aligned}

-  Coherent Oscillating fields:

   .. math::

      \begin{aligned}
      N_i' & = & -3 - \frac{\Gamma_i}{H} \nonumber \\
      R_i' & = & 0 \label{eq:COeq}\end{aligned}

As seen above, the equation for :math:`R_i = \rho_i/n_i` depends on
:math:`P_i/n_i`. A proper evaluation of this quantity requires knowledge
of the distribution :math:`F_i(p,t)`. However, for relativistic (or
massless) particles we have :math:`P_i = \rho_i/3`, as seen from
Eq.([beqs]), while for particles at rest we have :math:`P_i = 0`. Hence
:math:`F_i(p,t)` is only required to evaluate the
relativistic/non-relativistic transition, which corresponds to a
relatively small part of the evolution history of particle :math:`i`.
Nonetheless, to model this transition we approximate :math:`F_i` by a
thermal distribution and take :math:`T_i, \mu_i \ll m_i`, where
:math:`T_i` is the temperature of the particle (which can be different
from the thermal bath’s). Under these approximations we have:

.. math::

   \begin{aligned}
   \frac{P_i}{n_i} & = & T_i \nonumber \\
   \frac{\rho_i}{n_i} & = & T_i \left[ \frac{K_1(m_i/T_i)}{K_2(m_i/T_i)} \frac{m_i}{T_i} + 3 \right] \label{eq:p1}\end{aligned}

where :math:`K_{1,2}` are the modified Bessel functions. In particular,
if :math:`m_i/T_i \gg 1`:

.. math:: \frac{\rho_i}{n_i} \simeq T_i \left[\frac{3}{2} + \frac{m_i}{T_i}  + 3 \right] \To \frac{P_i}{n_i} = T_i = \frac{2 m_i}{3}\left( \frac{R_i}{m_i} -1 \right)

As shown above, for a given value of :math:`R_i = \rho_i/n_i`,
Eq.([eq:p1]) can be inverted to compute :math:`T_i` (:math:`=P_i/n_i`):

.. math:: \frac{P_i}{n_i} = T_i(R_i)

Since we are interested in the non-relativistic/relativistic transition,
we can expand the above expression around :math:`R_i/m_i = 1`, so
:math:`P_i/n_i` can be written as:

.. math:: \frac{P_i}{n_i} = \frac{2 m_i}{3}\left( \frac{R_i}{m_i} -1 \right) + m_i \sum_{n >1} a_n \left(\frac{R_i}{m_i} -1 \right)^n

where the coefficients :math:`a_n` can be numerically computed from
Eq.([eq:p1]). The above approximation should be valid for
:math:`m_i/T_i \gtrsim 1` (or :math:`R_i \gtrsim m_i`). On the other
hand, for :math:`m_i/T_i \ll 1` (or :math:`R_i \gg m_i`), we have the
relativistic regime, with :math:`P_i/n_i = R_i/3`. Therefore we can
approximate the :math:`P_i/n_i` function for all values of :math:`R_i`
by:

.. math::

   \frac{P_i}{n_i} = \left\{ \begin{array}{rl}
   & \frac{2 m_i}{3}\left( \frac{R_i}{m_i} -1 \right) + m_i \sum_{n >1} a_n \left(\frac{R_i}{m_i} -1 \right)^n  \mbox{ , for $R_i < \tilde{R}$} \\
   & \frac{R_i}{3}  \mbox{ , for $R_i > \tilde{R}$} 
   \end{array} \right. \label{Pfin}

where the coefficients :math:`a_n` are given by the numerical fit of
Eq.([eq:p1]) and :math:`\tilde{R}` is given by the matching of the two
solutions.

Finally, to solve Eqs.([eq:Sfin])-([eq:COeq]) we need to compute
:math:`H` according to Eq.([H]), which requires knowledge of the energy
densities for all particles (:math:`\rho_i`) and for the thermal bath
(:math:`\rho_R`). The former are directly obtained from :math:`N_i` and
:math:`R_i`, while the latter can be computed from :math:`N_S`:

.. math:: T = \left(\frac{g_{*S}(T_R)}{g_{*S}(T)}\right)^{1/3} T_R \exp[N_S/3 -x] \To \rho_R = \frac{\pi^2}{30} g_{*}(T) T^4

Eqs.([eq:Sfin])-([eq:COeq]), with the auxiliary equations for :math:`H`
(Eq.([H])) and :math:`P_i/n_i` (Eq.([Pfin])) form a set of closed
equations, which can be solved once the initial conditions for the
number density (:math:`n_i`), energy density (:math:`\rho_i`) and
entropy (:math:`S`) are given. For thermal fluids we assume:

.. math::

   \begin{aligned}
   n_i(T_R) & = & \left\{ 
   \begin{array}{ll} 
   0 & , \mbox{ if $\sigv_i \bar{n}_i/H|_{T=T_R} < 10$} \\
   \bar{n}_i(T_R) & , \mbox{ if $\sigv_i \bar{n}_i/H|_{T=T_R} > 10$} 
   \end{array} \right. \label{ni0TP} \\
   \frac{\rho_i}{n_i}(T_R) & = & \frac{\bar{\rho}_i}{\bar{n}_i}(T_R)\end{aligned}

where :math:`\bar{\rho}_i` is the equilibrium energy density (with zero
chemical potential) for the particle :math:`i`. While for coherent
oscillating fluids the initial condition is set at the beginning of
oscillations:

.. math::

   \begin{aligned}
   n_i(T^{osc}_i) & = &\frac{\rho_i^{0}}{m_i(T^{osc}_i)} \\
   \frac{\rho_i}{n_i}(T^{osc}_i) & = & m_i\end{aligned}

where :math:`T^{osc}_i` is the oscillation temperature, given by
:math:`3H(T^{osc}_i) = m_i(T^{osc}_i)` and :math:`\rho_i^{0}` the
initial energy density for oscillations.

Finally, the initial condition for the entropy :math:`S` is trivially
obtained, once we assume a radiation dominated universe at
:math:`T=T_R`:

.. math:: S(T_R) = \frac{2 \pi^2}{45} g_*(T_R) T_R^3 R_0^3

Code
====

Here we describe how the above formalism is implemented in a numerical
code for solving the coupled Boltzmann equations. In Sec.[sec:In] we
describe how the input for specific models should be defined. Then, in
Sec.[sec:Main] we outline the procedure used to solve the Boltzmann
equations and to treat some of the discrete transitions required by the
formalism described above. Finally, in Sec.[sec:Out] we describe what is
the output of the code and how it can be controlled by the user.

Input
-----

In order to solve the Boltzmann equations for a particular model, the
user has to define in the main program the BSM components which must be
evolved. The components are defined using the

the SUBROUTINE INPUTBOLTZ(T), which, for a given (thermal bath)
temperature :math:`T`, fills the COMMON BLOCK:

| COMMON/INBOLTZ/BR(NP,NP),DEGF(NP),MASS(NP),GAM(NP),
| SIGV(NP),C(NP),COHOSC(NP),TRH,NCOMPS,LABEL(NP)

where NP = 20 and

-  NCOMPS (:math:`\leq 10`) = the number of particles (the first
   component must be radiation)

-  MASS(i) = mass for particle i (can be temperature dependent, as in
   the axion case)

-  DEGF(i) = +-number of degrees of freedom for particle i. A plus sign
   should be used for bosons, while a minus should be used for fermions,
   i.e. DEGF=-2 for neutralinos and DEGF=1 for axions. The value for the
   i=1 component (radiation) is never used, since the number of degrees
   of freedom in this case in given by the function GSTAR(T).

-  GAM(i) = decay width for particle i, in its rest frame.

-  BR(i,j) = branching ratio for the decay :math:`i \to j + X`,
   including the multiplicity factor, if the i particle decays into
   multiple j’s.

-  BR(i,1) = fraction of energy per i particle injected in the radiation
   fluid.

-  SIGV(i) = thermal averaged cross-section for the annihilation of i
   particles, as defined in the previous section.

-  C(i) = additional source term for particle i, as defined in the
   previous section.

-  COHOSC(i) = initial energy density for coherent oscillating
   particles. Must be zero for thermal (non-oscillating) components.

-  TRH = re-heat temperature.

-  LABEL(i) = label for particle i (optional)

Main Code
---------

Once the INPUTBOLTZ subroutine is provided, the user can compute the
solution for the Boltzmann equations from T=TRH to T=TF, calling:

| CALL INPUTBOLTZ(TRH) ! (iniatilization)
| CALL EQSBOLTZ(TF,IOUT) ! (compute solution)

where, if IOUT\ :math:`>0`, the scale factor (:math:`R`) and energy
densities as a function of :math:`T` are written to UNIT=IOUT. If TF =
0, the evolution proceeds until all unstable particles have decayed
and/or all coherent oscillating components have started to oscillate.
Before calling EQSBOLTZ, the user must define the parameters which
regulate the precision of the procedure, given by the BLOCK DEPARS:

COMMON/DEPARS/EPS,DX0,STEP,IERROR

where EPS is the relative precision for the solution :math:`N_i(TF)`,
DX0 is the :math:`x` interval for printing the solutions in IOUT and the
maximum :math:`\Delta x` step and STEP is the initial :math:`\Delta x`
step for the evolution. Failure to solve the equations (most likely due
to numerical instabilities) is indicated by IERROR\ :math:`<0`.

Specific components can be turned off using the COMMON BLOCK:

COMMON/SWITCHES/TURNOFF(NP)

If TURNOFF(I)=.TRUE., the :math:`i-`\ component will not be included in
the evolution of the Boltzmann equations.

The EQSBOLTZ is the main subroutine used to solve the equations, once
the appropriate input has been defined. Its main steps are:

#. Set initial conditions at :math:`T=TRH`: check which thermal
   particles are coupled/decoupled to the thermal bath and if coherent
   oscillating fluids are already oscillating at T=TRH. Then it sets the
   initial number densities and temperatures for each component, as
   defined in the previous section. Set X1=1 and X2=X1+DX0.

#. Solve the equations between X1 and X2.

#. Check if a particle has decayed. If the particle :math:`i` satisfies

   .. math:: \Gamma_i/H > 100\;\; {\rm and}\;\; \min_{j \neq i}(\rho_i/\rho_j)< 10^{-3}

   the particle is neglected from here on. The decay temperature
   (:math:`T_D`) is defined by the sudden decay approximation:
   :math:`\Gamma_i/\gamma_i = H(T_D)`, where :math:`\gamma_i` is the
   boost factor (:math:`\gamma_i \equiv \langle E_i \rangle/m_i`).
   Although this temperature is printed out in the output, it is never
   used in the code.

#. Check if a particle has decoupled from the thermal bath or started to
   oscillate in the interval (X1,X2). Decoupling is assumed if
   :math:`\sigv_i \bar{n}_i < H/10`, which also defines the freeze-out
   temperature. The oscillation temperature is given by
   :math:`3 H(T_{osc}) = m_i(T_{osc})` and defines the beginning of
   evolution for the oscillating components.

#. If a component has started to oscillate or if it has decoupled, loop
   over this interval with smaller steps until the decoupling or
   oscillation temperature converges (:math:`\Delta T_i/T_i < 0.1`).

#. Write temperature, scale factor and energy densities to IOUT, if
   IOUT\ :math:`>0`.

#. Set X1=X2 and X2=X1+DX0 and return to point 2. until :math:`T<TF` (or
   all unstable particles have decayed and all oscillating fluids have
   oscillated, if TF\ :math:`\leq 0`).

Output
------

The standard information printed after solving the Boltzmann equations
gives:

-  Freeze-out temperatures (:math:`T_{fr}`) for each thermal component.
   As mentioned in the last section, :math:`T_{fr}` is given by the
   decoupling condition: :math:`\sigv_i \bar{n}_i = H/10`. Since the
   decoupling is a continuous process, :math:`T_{fr}` is just an
   estimate for the decoupling temperature.

-  Decay temperatures (:math:`T_D`) for each thermal component. Once
   again the decay process is continuous and :math:`T_D` given in the
   print out is estimated by the sudden decay approximation
   (:math:`\Gamma_i/\gamma_i = H(T_D)`).

-  Oscillation temperature (:math:`T_{osc}`).

-  Entropy ratio (:math:`S/S_0`). In case of entropy injection from
   decays of unstable particles, :math:`S/S_0 > 1`.

-  Relic densities (:math:`\Omega_i h^2`) at :math:`T_0 = 2.725` K. In
   order to consistently compute the relic densities today, we evolve
   :math:`R_i = \rho_i/n_i` from :math:`TF` to :math:`T_0` assuming a
   trivial universe expansion:

   .. math:: R_i' = - 3 \frac{P_i}{n_i}

   Note that the result obtained above is insensitive to :math:`H`, so
   it does not matter if there is a transition from a radiation
   dominated to a matter dominated (or dark energy dominated) universe
   between :math:`TF` and :math:`T_0`. Once :math:`R_i(T_0)` is
   obtained, the relic density is given by:

   .. math:: \Omega_i h^2 = n_i(TF) \times \frac{g_{*S}(T_0) T_0^3}{g_{*S}(TF) TF^3} \times \frac{R_i(T_0)}{\rho_c/h^2}

-  Relic densities before decay. It may be relevant to compute the relic
   densities of an unstable particle as it would be given if it had not
   decayed. In particular, this value can be used to impose BBN bounds
   on the decays. This quantity is computed after the particle becomes
   non-relativistic and well before the decay starts
   (:math:`\Gamma_i/H(T) = 1/10`) and is given by:

   .. math:: \tilde{\Omega}_i h^2 = \frac{\rho_i(T)}{s(T)} \times \frac{s(T_0)}{\rho_c/h^2}

   Note that the above expression assumes a radiation dominated universe
   from :math:`T` to :math:`T_0` and should be used with caution.

-  Effective number of (new) neutrinos (:math:`\Delta N_{eff}`). Since
   neutrinos are still coupled for :math:`T > 1` MeV, this quantity is
   only compute below this temperature. :math:`\Delta N_{eff}` is given
   by:

   .. math:: \Delta N_{eff}(T) = \frac{\rho_{DR}(T)}{\rho_{\nu}}

   where :math:`\rho_{DR}` is the total energy density of relativistic
   particles (excluding radiation and neutrinos) and :math:`\rho_{\nu}`
   is the energy density of neutrinos after they freeze-out:

   .. math:: \rho_{DR} = \sum_{R_i/m_i > 2} \rho_i \mbox{ and } \rho_{\nu} = \frac{\pi^2}{15}\frac{7}{8}\left(\frac{4}{11}\right)^{4/3} T^4

   Note that :math:`\Delta N_{eff}` is in general a function of
   temperature, since :math:`\rho_{DR}` will decrease if massive
   particles become non-relativistic below 1 MeV.

Furthermore, if :math:`IOUT>0`, the scale factor (:math:`R`), the energy
densities and :math:`\Delta N_{eff}` are printed as a function of
:math:`T` in UNIT=IOUT. Also, the following quantities are stored in
COMMON BLOCKs:

-  Final relic densities and entropy ratio:

   COMMON/OUTPUT/OMEGA(NP),RS

-  Decoupling, oscillation and decay temperatures:

   COMMON/TEMPS/TDEC(NP),TOSC(NP),TDCAY(NP)

-  Relic density of unstable particles before decay
   (:math:`\tilde{\Omega} h^2`), temperature at the end of entropy
   injection (if any) and effective number of new neutrinos
   (:math:`\Delta N_{eff}`) *at the final temperature :math:`TF`*:

   COMMON/BBNINFO/UMEGA(NP),TSTAB,DNeff

PQMSSM
======

In order to apply the above formalism to the PQMSSM we need to define:

-  Masses: :math:`m_{\tz_1}`, :math:`m_{\ta}`, :math:`m_s`,
   :math:`m_{a}(T)`, :math:`m_{\tG}`

-  Decay Width: :math:`\Gamma_i`, with :math:`i = \tz_1,\ta,s,\tG`,

-  Branching Ratios: :math:`BR(i,j)` and :math:`BR(i,1)`, with
   :math:`i,j = \tz_1,\ta,s,\tG`,

-  Annihilation cross-sections: :math:`\sigv_i`, with
   :math:`i = a,\tz_1,\ta,s,\tG`,

-  Additional Source terms: :math:`C_i`, with
   :math:`i = a,\tz_1,\ta,s,\tG`,

-  Initial energy density for coherent oscillating fields:
   :math:`\rho_i^0`, with :math:`i=a,s`

-  and remaining SUSY spectrum (for computation of :math:`g_*` and
   gravitino/axino decays)

Below we describe how the quantities :math:`\sigv_i`, :math:`\Gamma_i`,
:math:`BR(i,j)`, :math:`C_i` and :math:`\rho_i^0` are computed.

:math:`\sigv`
-------------

The annihilation cross-section for the neutralino component is computed
in the SUBROUTINE ZSIG(T,MZ), where T is the temperature and MZ is the
neutralino mass. For efficiency purposes the calculation of
:math:`\sigv_{\tz_1}` can be controlled through the COMMON BLOCK:

COMMON/INSIGMA/INOMGZ,INMZ,INFLAG,INDATA

The main options are set by INFLAG:

-  If INFLAG=166, the subroutine returns a constant value for
   :math:`\sigv_{\tz_1}`, given by:

   .. math:: \sigv_{\tz_1} = 1.7\times 10^{-10}/\Omega_{\tz_1} h^2

   where the value for :math:`\Omega_{\tz_1} h^2` should be set in
   INOMGZ.

-  If INFLAG=266, then

   -  If INDATA :math:`<0`: generate file with :math:`\sigv_{\tz_1}`
      values for :math:`3\times 10^{-5} < T/m_{\tz_1} < 2` in
      UNIT=ABS(INDATA) for future extrapolation and set INDATA =
      ABS(INDATA). The number of points and specific values of :math:`T`
      are chosen as to properly describe the shape of
      :math:`\sigv_{\tz_1}`. The file must be open by the main program.

   -  If INDATA :math:`>0`: use values in UNIT=INDATA for a linear
      extrapolation in :math:`\log(T/m_{\tz_1})`

-  If INFLAG :math:`\neq` 166 and 266: compute :math:`\sigv_{\tz_1}(T)`
   (full integration). Due to convergence issues,
   :math:`\sigv_{\tz_1}(T > m_{\tz_1}/2) \equiv \sigv_{\tz_1}(T=m_{\tz_1}/2)`
   and
   :math:`\sigv_{\tz_1}(T < m_{\tz_1}/5\times 10^{-5}) \equiv \sigv_{\tz_1}(T=m_{\tz_1}/5\times 10^{-5})`.

The thermal production rate for axions, saxions, axinos and gravitinos
was computed in Refs.. However, the expressions derived in Refs. are
only valid for out of equilibrium production, where the (out of
equilibrium) thermal production rate (:math:`W_i`) is defined by:

.. math:: \frac{d n_i}{dt} + 3 H = W_i \label{eq:prodr}

In the out of equilibrium regime we have :math:`\bar{n}_i \gg n_i`, so
Eq.([eq:nieq]) becomes:

.. math:: \frac{d n_i}{dt} + 3 H = \bar{n}_i^2 \sigv_i

Comparing the above equation to Eq.([eq:prodr]), we identify:

.. math:: \sigv_i = \frac{W_i}{\bar{n}_i^2} \label{eq:sigw}

Although the above relation is only exact for the out of equilibrium
regime, we use it for all values of temperature. This is only a poor
approximation if the reheat temperature (:math:`T_R`) is very close to
the decoupling temperature (:math:`T_{dec}`). Since for
:math:`T_R \gtrsim T_{dec}`, Eq.([eq:nieq]) gives
:math:`n_i = \bar{n}_i` (independent of :math:`\sigv_i`), while for
:math:`T_R \lesssim T_{dec}` Eq.([eq:sigw]) is exact.

Using the expressions and Eq.([eq:sigw]), we obtain [2]_:

.. math::

   \begin{aligned}
   \sigv_{a} & = & \frac{9}{128 \pi^3 \xi(3)} \frac{g_s^6}{f_a^2}\ln(\frac{1.0126}{g_s})\\
   \sigv_{s} & = & \sigv_{a}\\
   \sigv_{\ta} & = & \frac{1}{576 \pi^3 \xi(3)^2} \frac{g_s^4}{f_a^2}F(g_s)\\
   \sigv_{\tG} & = & \frac{1.37}{M_P^2}\times \left[ 72 g_s^2 \ln(1.271/g_s)(1+\frac{M_3^2}{3 m_{\tG}^2})\right. \\
   &+& \left. 27 g^2 \ln(1.312/g)(1+\frac{M_2^2}{3 m_{\tG}^2}) + 11 g'^2 \ln(1.266/g')(1+\frac{M_1^2}{3 m_{\tG}^2}) \right]\end{aligned}

where :math:`F(g_s)` is numerically obatined from Ref.:

.. math:: F(x) = -0.365771 + 9.38897 x + 27.7315 x^2 - 20.1012 x^3 + 5.153 x^4

We note that while the expression for :math:`\sigv_{\ta}` includes
thermal decays such as :math:`g \to \tg + \ta` and is valid for all
values of :math:`g_s`, the expressions for :math:`\sigv_{a,s}`
correspond only to the hard themal loop calculation and assume
:math:`g_s \ll 1` (or :math:`T_R \gg 10^6` GeV). Nonetheless we use
these expressions for all values of :math:`g_s` (:math:`T_R`).

:math:`\Gamma`,\ :math:`BR`
---------------------------

The decay rates are computed through the SUBROUTINES:

| AXINOBR(MAXINO,FA,AXLT,AXWD,AXVIS)
| GRAVITINOBR(MGT,GLT,GWD)
| Z1BRS(MZ1,FA,Z1B,Z1WD,Z1LT)
| SAXIONBR(MSAXION,FA,XI,SAXWD,SAXLT)

where the axino, gravitino, neutralino and saxion decay rates are given
by AXWD, GWD, Z1WD and SAXWD, respectively. Neutralinos and axinos
decays into gravitinos (if kinematically allowed) are NOT included. So
the gravitino LSP case is NOT presently included. For the case of an
axino LSP, Z1BRS assumes MAXINO=0. Therefore, if
:math:`m_{\tz_1} < m_{\ta}`, the user should set Z1WD=0.

The branching ratios are given by:

.. math::

   \begin{aligned}
   BR(s \to X) & = &  1 - BR(s \to aa) - BR(s \to \ta\ta)\\
   BR(s \to \tilde{Z}_1) & = & 2\times BR(s \to \tg\tg) \\
   BR(s \to a) & = & 2\times BR(s \to aa)\end{aligned}

-  If :math:`m_{\tz_1} > m_{\ta}`:

   .. math::

      \begin{aligned}
      BR(\tz_1 \to \ta) & = &  1 \\
      BR(\tz_1 \to X) & = &  1\end{aligned}

-  If :math:`m_{\ta} > m_{\tz_1}`:

   .. math::

      \begin{aligned}
      BR(\ta \to \tz_1) & = &  1 \\
      BR(\ta \to X) & = &  1 \\\end{aligned}

-  If :math:`m_{\tG} > m_{\tz_1}`:

   .. math::

      \begin{aligned}
      BR(\tG \to \tz_1) & = &  1 \\
      BR(\tG \to X) & = &  1\end{aligned}

-  If :math:`m_{\tG} < m_{\tz_1}` (but :math:`m_{\tG} > m_{\ta}`, since
   :math:`\tG` LSP is not allowed!):

   .. math::

      \begin{aligned}
      BR(\tG \to \ta) & = &  1 \\
      BR(\tG \to a) & = &  1 \end{aligned}

and all other BRs are zero.

As discussed above, the branching ratio :math:`BR(i \to X)` is defined
as the fracion of energy injected in the thermal bath. For simplicity we
assume that most of the energy from decays goes into radiation, so
:math:`BR(i \to X) \sim 1`, except for saxion decays to axions and
axinos and gravitino decays to axion + axino.

:math:`C`
---------

The source term :math:`C_i` can be used to include other processes that
can not be described as decays or annihilations. However these are not
relevant for the PQMSSM. Therefore we set all these to zero.

:math:`\rho^0`
--------------

Finally, the initial energy densities for the oscillating saxion and
axion fields are given by:

.. math::

   \begin{aligned}
   \rho^0_a & = & 1.44  \frac{m_a(T)^2 f_a^2 \theta_i^2 }{2} f(\theta_i)^{7/6}\\
   \rho^0_s & = & \min\left[1.9 \times 10^{-8}\left(\frac{2\pi^2 g_*(T_R) T_R^3}{45}\right)\left(\frac{T_R}{10^5}\right)\left(\frac{s_i}{10^{12}}\right)^2,\frac{m_s^2 s_i^2}{2}\right]    \end{aligned}

where :math:`f(\theta_i) = \ln[e/(1-\theta_i^2/\pi^2)]` and
:math:`f_a \theta_i` and :math:`s_i` are the initial axion and saxion
amplitudes. The definition of :math:`\rho^0_s` accounts for the
possibility of beginning of saxion oscillations during inflation
(:math:`T_R < T_{osc}`).

99 E. Kolb and M. Turner, *The Early Universe*, Addison-Wesley Pub.
(1990). M. Kawasaki, G. Steigman and H. Kang, *Cosmological evolution of
generic early-decaying particles and their daughters*, Nucl.Phys. B402,
323-348 (1993). M. Kawasaki, G. Steigman and H. Kang, *Cosmological
evolution of an early-decaying particle*, Nucl.Phys. B403, 671-706
(1993). P. Graff and F. Steffen, *Thermal axion production in the
primordial quark-gluon plasma*, Phys.Rev. D 83, 070511 (2011). P. Graff
and F. Steffen, *Axions and saxions from the primordial supersymmetric
plasma and extra radiation signatures*, arXiv:1208.2951 (2012). A.
Strumia, *Thermal production of axino Dark Matter*, JHEP 06, 036 (2010).
J. Pradler and F. Steffen, *Constraints on the reheating temperature in
gravitino dark matter scenarios*, Phys.Lett. **B 648**, 224 (2007).

.. [1]
   This approximation is only valid for particles with a thermal
   distribution. However, since the annihilation term is responsible for
   keeping the particle :math:`i` in thermal equilibrium with the
   plasma, it is reasonable to assume a thermal distribution for
   :math:`i` while the annihilation term is relevant.

.. [2]
   According to Ref., the axion and saxion thermal rates are identical
   in supersymmetric models.
