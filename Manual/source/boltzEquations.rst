
.. _boltzEqs:

Boltzmann Equations
===================



General Formalism and Approximations
------------------------------------

The general Boltzmann equation for the number distribution of a particle
species can be written as (assuming isotropy):

.. math:: 
   :label: eq:d1
   :nowrap:
   
   {\frac{\partial F_{i}}{\partial t}} -H p {\frac{\partial F_{i}}{\partial p}} = C_{i}[F_{i},F_{j},p] 

where :math:`F_{i}(p)` is the number distribution of particle :math:`i`
as function of momentum :math:`p`, :math:`C` represents a source/sink
term and :math:`H` is the Hubble constant:

.. math:: H = \sqrt{\frac{8 \pi}{3} \frac{\rho_T}{M_P^2}} \label{H}

with :math:`M_{P} = 1.22\times 10^{19}` GeV and
:math:`\rho_T = \sum_{i} \rho_i`. The number, energy and pressure
densities are given in terms of :math:`F_{i}` as:

.. math::
         
   \begin{align}
   n_{i}(t) & = \int \frac{dp}{2 \pi^2} p^2 F_i(p) \nonumber \\ 
   \rho_{i}(t) & = \int \frac{dp}{2 \pi^2} p^2 E_i F_i(p) \label{beqs}\\
   P_{i}(t) & = \frac{1}{3} \int \frac{dp}{2 \pi^2} \frac{p^4}{E_i} F_i(p) \nonumber\end{align}

where :math:`m_i` is the mass of particle :math:`i` and
:math:`E_i = \sqrt{p_i^2 + m_i^2}`. Using Eq. :eq:`eq:d1`, we obtain the
following equations for the number and energy densities:

.. math::

   \begin{aligned}
   {\frac{d n_i}{d t}} + 3H n_i & = \int \frac{dp}{2 \pi^2} p^2 C_i \nonumber \\
   {\frac{d \rho_i}{d t}} + 3H (\rho_i + P_i) & = \int \frac{dp}{2 \pi^2} p^2 E_i C_i \label{eq:meqs}\end{aligned}

The collision term, :math:`C_i`, for the process
:math:`i + f + \ldots \leftrightarrow a
+ b + c + \ldots` is given by:

.. math::

   \begin{aligned}
   C_i & = \frac{1}{E_i} \int \prod_{j,a} \frac{d^3 p_j}{2 E_j (2 \pi)^3}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i + p_j + \ldots - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   &\times \left[(1 \pm f_a) (1 \pm
   f_b)\ldots f_i f_j\ldots - f_a f_b \ldots (1 \pm f_i)(1 \pm f_j)\ldots \right]\end{aligned}

where the plus (minus) sign is for bosons (fermions). Below we always
assume :math:`f_{i,j,a,..} \ll 1`, so:

.. math::

   \begin{aligned}
   C_i & \simeq \frac{1}{E_i} \int \prod_{j,a} \frac{d^3 p_j}{2 E_j (2 \pi)^3}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i + p_j + \ldots - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   &\times \left[f_i f_j\ldots - f_a f_b \ldots \right]\end{aligned}

We will assume that :math:`C` is given by:

.. math:: C = C_{dec} + C_{prod} + C_{ann}

where :math:`C_{dec}` contains the contributions from decays and inverse
decays (:math:`i \leftrightarrow a + b + \ldots`), :math:`C_{prod}`
contains the contributions from decay injection and inverse decay
injection (:math:`a \leftrightarrow i + b + \ldots`) and :math:`C_{ann}`
from annihilations with the thermal plasma
(:math:`i + i \leftrightarrow a + b`). Below we compute each term
separately, under some assumptions.

.. _annTerm:

Annihilation Term
-----------------

The annihilation term :math:`C_{ann}` for the
:math:`i + j \leftrightarrow a + b` process is given by:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{ann} = \int d\Pi_{i} d\Pi_{j} d\Pi_{a}
   d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \left[ f_a f_b - f_i f_j \right]

where :math:`d\Pi_{i} = d^{3} p_i/((2\pi)^3 2 E_i)`. Since we are
ultimately interested in Eqs.([eq:meqs]) for the number and energy
densities, we will consider the following integral:

.. math::

   \int \frac{dp}{2 \pi^2} p^2 C_{ann}  E_i^{\alpha} = \int d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 
   \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \\
    \times \left[ f_a f_b - f_i f_j \right] E_i^{\alpha}

where :math:`\alpha = 0 (1)` for the number (energy) density. Here we
assume that the distributions can be approximated by:

.. math:: f_i \simeq \exp(-(E_i - \mu_i)/T)

so the annihilation term can then be written as:

.. math::

   \begin{aligned}
   & \int \frac{dp}{2 \pi^2} p^2 C_{ann}  E_i^{\alpha} =  -\left( \exp((\mu_i + \mu_j)/T) -\exp((\mu_a + \mu_b)/T)\right) \nonumber \\
    & \times \int  d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha} \nonumber\end{aligned}

where above we have used conservation of energy
(:math:`E_i + E_j = E_a + E_b`). Since for the cases of interest the
equilibrium distributions have zero chemical potential, we have:

.. math:: \frac{n_i}{\bar{n}_i} = \exp(\mu_i/T)

so:

.. math::

   \begin{aligned}
   & \int \frac{dp}{2 \pi^2} p^2 C_{ann} E_i^{\alpha} = -\left( \frac{n_i n_j}{\bar{n}_i \bar{n}_j} - \frac{n_a n_b}{\bar{n}_a \bar{n}_b}\right) \nonumber \\
    & \times \int  d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha} \nonumber\end{aligned}

In particular, for the process :math:`i + i \leftrightarrow a + b`,
where :math:`a` and :math:`b` are in thermal equilibrium
(:math:`\mu_a = \mu_b = 0`):

.. math::

   \begin{aligned}
   & \int \frac{dp}{2 \pi^2} p^2 C_{ann} E_i^{\alpha} =  -\left( \frac{n_i^2}{\bar{n}_i^2} - 1 \right) \nonumber \\
   &  \times \int d\Pi_{i} d\Pi_{j} d\Pi_{a} d\Pi_{b} (2 \pi)^4 \delta^{(4)}(p_i + p_j - p_a - p_b) |M|^2 \exp(-(E_i + E_j)/T) \times E_i^{\alpha}  \nonumber \\
    & = -\left( n_i^2 - \bar{n}_i^2 \right) \langle \sigma v E_i^{\alpha} \rangle\end{aligned}

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

.. _decayTerm:

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
   \bar{n}_b \ldots} \exp(-E_i/T)  \\
   =  \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \bar{f}_{i}

where we used conservation of energy (:math:`E_a + E_b + \ldots = E_i`)
and :math:`\bar{f}_i` is the equilibrium distribution for the species
:math:`i`. Hence we can write Eq.([eq:dec0]) as:

.. math::

   \begin{aligned}
   C_{dec} & \simeq \left[f_i - \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots}
   \bar{f}_{i} \right] \frac{1}{E_i} \int \prod_{a}
   \frac{d^3 p_a}{2 E_a (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_i - p_a - p_b
   \ldots\right) |\mathcal{M}|^2 \nonumber \\
   & = \mathcal{B}_{ab\ldots} \frac{\Gamma_i m_i}{E_i} \left[f_i -
   \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \bar{f}_{i} \right] \end{aligned}

where :math:`\Gamma_i` is the width for :math:`i` and
:math:`\mathcal{B}_{ab\ldots} \equiv BR(i \to a + b + \ldots)`

Once again we consider the integral:

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{dec}(p) E_i^{\alpha} = 
    & - \Gamma_i \int \frac{dp}{2 \pi^2} p^2 \frac{m_i}{E_i} f_i E_i^{\alpha}
    \nonumber \\
    & + \sum_{i \; decays} \mathcal{B}_{ab\ldots}
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

   \int \frac{dp}{2 \pi^2} p^2 C_{dec}(p) E_i^{\alpha} = -\Gamma_i m_i 
   \left\{ \begin{array}{ll} \frac{n_i}{\rho_i}\left( n_i -
    \bar{n}_i \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \right)   &
    \mbox{, for $\alpha = 0$}  \\
    n_i - \bar{n}_i \sum \mathcal{B}_{ab\ldots}
    \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots}  & \mbox{, for $\alpha = 1$}
   \end{array} \right.

where :math:`\mathcal{B}_{ab\ldots} \equiv BR(i\to a+b+\ldots)`.

.. _prodTerm:

Production Term
---------------

The decay and inverse decay of other particles
(:math:`a \to i + b + \ldots`) can also affect the species :math:`i`.
The contribution from these terms we label :math:`C_{prod}`, which is
given by:

.. math::

   C_{prod} \simeq \frac{1}{E_i} \int \frac{d^3 p_a}{2 E_a (2
   \pi)^3} \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_a
   - p_i - p_b \ldots\right) |\mathcal{M}|^2\\
   \times \left[f_a - f_i f_b \ldots \right]

Using the same approximations of the previous section, we write:

.. math::

   f_i f_b\ldots \simeq  \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots}
   e^{-E_a/T} = \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots}
   \bar{f}_{a}

Hence:

.. math::

   C_{prod} = \frac{1}{E_i} \int \frac{d^3 p_a}{2 E_a (2 \pi)^3} \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} 
   (2 \pi)^4 \delta^{4}\left(p_a - p_i - p_b \ldots\right) |\mathcal{M}|^2 \\
   \times \left(f_a - \bar{f}_a \frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right)

and

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) E_i^\alpha & = 
   \int \frac{d^3 p_a}{E_a (2 \pi)^3} \left(f_a - \bar{f}_a \frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right) \nonumber \\
   & \times \frac{d^3 p E_i^{\alpha}}{2 E_i (2 \pi)^3}
   \prod_{b} \frac{d^3 p_b}{2 E_b (2 \pi)^3} (2 \pi)^4 \delta^{4}\left(p_a - p_i - p_b \ldots\right) |\mathcal{M}|^2
   \label{eq:prod2}\end{aligned}

with :math:`\alpha = 0 (1)` for the contribution to the number (energy)
density equation. For :math:`\alpha = 0` we obtain:

.. math::

   \begin{aligned}
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) & = \Gamma_a  \mathcal{B}_{i} m_a 
   \int \frac{d^3 p_a}{E_a (2 \pi)^3} \left(f_a - \bar{f}_a \sum_b
   \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}}\frac{n_i n_b \ldots}{\bar{n}_i
   \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & = \Gamma_a \mathcal{B}_{i} m_a \frac{n_a}{\rho_a} \left( n_a - \bar{n}_a
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
   \int \frac{dp}{2 \pi^2} p^2 C_{prod}(p) E_i & \simeq
   \Gamma_a \mathcal{B}_{i}  \frac{m_a}{2} \int \frac{d^3 p_a}{(2
   \pi)^3} \left(f_a - \bar{f}_a \sum_b
   \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}}
    \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & = \Gamma_a \mathcal{B}_{i}  \frac{m_a}{2} \left( n_a -
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

.. _bequations:

Number and Energy Density Equations
-----------------------------------

Using the results of Eqs.([eq:collfin]), ([eq:decfin]) and
([eq:prodfin]) in the Boltzmann equations for :math:`n_i` and
:math:`\rho_i` (Eq.([eq:meqs])), we obtain:

.. math::

   \begin{aligned}
   {\frac{d n_i}{d t}} + 3H n_i  & =  \left( \bar{n}_i^2 - n_i^2 \right) \langle \sigma
   v \rangle - \Gamma_i m_i \frac{n_i}{\rho_i}\left(n_i - \bar{n}_i \sum_{i\to\ldots}
   \mathcal{B}_{ab\ldots} \frac{n_a n_b \ldots}{\bar{n}_a \bar{n}_b \ldots} \right)
   \nonumber
   \\
   & + \sum_a 
   \Gamma_a \mathcal{B}_i m_a \frac{n_a}{\rho_a} \left(n_a - \bar{n}_a \sum_{a \to
   i\ldots} \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i n_b \ldots}{\bar{n}_i \bar{n}_b \ldots} \right)  + C_{i}(T) \label{eq:nieq} \\
   {\frac{d \rho_i}{d t}} + 3H (\rho_i + P_i) & = \left( \bar{n}_i^2 - n_i^2 \right)
   \langle \sigma v \rangle \frac{\rho_i}{n_i} - \Gamma_i m_i \left( n_i -
   \bar{n}_i \sum_{i\to\ldots} \mathcal{B}_{ab\ldots} \frac{n_a n_b\ldots}{\bar{n}_a
   \bar{n}_b\ldots}\right) \nonumber \\
    & + \sum_a \Gamma_a  \mathcal{B}_i \frac{m_a}{2} \left( n_a -
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

   \begin{aligned}
   {\frac{d \rho_i/n_i}{d t}} \equiv {\frac{d R_i}{d t}} & = -3 H \frac{P_i}{n_i} \\
   & + \sum_{a} \mathcal{B}_{i} \frac{\Gamma_a m_a}{n_i} \left( \frac{1}{2} - \frac{n_a}{\rho_a} \frac{\rho_i}{n_i} \right) \left(n_a -
   \bar{n}_a \sum_{a \to i\ldots} \frac{\mathcal{B}_{ib\ldots}}{\mathcal{B}_{i}} \frac{n_i
    n_b..}{\bar{n}_i \bar{n}_b..}\right) \label{eq:Rieq} \end{aligned}

Besides the above equations, it is useful to consider the evolution
equation for entropy:

.. math:: dS \equiv \frac{dQ^{dec}}{T}

where :math:`dQ^{dec}` is the net energy injected from decays. With the
above definition we have:

.. math::

   \begin{aligned}
   \dot{S} & = \frac{1}{T}\sum_i BR(i,X)
   \frac{d\left(R^3 \rho_i\right)^{dec}}{dt}  \nonumber \\
   \Rightarrow \dot{S} & = \frac{R^3}{T}\sum_i BR(i,X)
   \Gamma_i m_i\left(n_i - \bar{n}_i \sum_{i\to\ldots} \mathcal{B}_{ab\ldots} \frac{n_a n_b\ldots}{\bar{n}_a
   \bar{n}_b\ldots} \right) \label{Seq}\end{aligned}

where :math:`R` is the scale factor and :math:`BR(i,X)` is the fraction
of energy injected in the thermal bath from :math:`i` decays.

The above expressions can be written in a more compact form if we define
the following ”effective thermal densities” and ”effective BR”:

.. math::

   \begin{aligned}
   \mathcal{N}^{th}_{X} & \equiv  \bar{n}_X \sum_{X \to \ldots} BR(X \to 1 + 2 +
   \ldots)
   \prod_{k}
   \frac{n_k}{\bar{n}_k} \nonumber \\
   \mathcal{N}^{th}_{XY} & \equiv \frac{\bar{n}_X}{\mathcal{B}^{eff}_{XY}}
   \sum_{X \to Y + \ldots} g_Y BR(X \to g_Y Y + 1 + \ldots)
   \left(\frac{n_Y}{\bar{n}_Y}\right)^{g_Y} \prod_{k} \frac{n_k}{\bar{n}_k}
   \nonumber \\
   \mathcal{B}^{eff}_{XY} & \equiv \sum_{X \to Y + \ldots} g_Y BR(X \to g_Y Y +
   1+\ldots) \nonumber\end{aligned}

where :math:`g_Y` is the :math:`Y` multiplicity in the final state of
:math:`X` decays. In addition, defining:

.. math:: x = \ln(R/R_0),\;\; N_i = \ln(n_i/s_0),\;\; {\rm and}\;\; N_S = \ln(S/S_0)

we can write Eqs.([Seq]), ([eq:nieq]) and ([eq:Rieq]) as:

.. math::

   \begin{aligned}
   N_S' & = \frac{e^{(3 x - N_S)}}{HT} \sum_{i} BR(i,X) \Gamma_i m_i \left(n_i -
   \mathcal{N}_{i}^{th} \right) 
   \label{Seqb} \\
   N_i' & = -3 + \frac{{\langle \sigma v \rangle}_i}{H} n_i [\left(\frac{\bar{n}_i}{n_i}\right)^2
   -1] -  \frac{\Gamma_i}{H} \frac{m_i}{R_i}\left(1 -
   \frac{\mathcal{N}_{i}^{th}}{n_i} \right) \nonumber  \\
    & + \sum_{a} \mathcal{B}_{ai}^{eff} \frac{\Gamma_a}{H}
    \frac{m_a}{R_a}\left(\frac{n_a}{n_i} - \frac{\mathcal{N}_{ai}^{th}}{n_i}
     \right)
    \\
   R_i' & =  -3 \frac{P_i}{n_i} + \sum_{a} \mathcal{B}_{ai}^{eff}
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
   N_i' & = -3 - \frac{\Gamma_i}{H}  \nonumber \\
   R_i'& = 0 \label{Nico}\end{aligned}

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
      N_i'& = -3 + \frac{{\langle \sigma v \rangle}_i}{H} n_i [\left(\frac{\bar{n}_i}{n_i}\right)^2
      -1] -  \frac{\Gamma_i}{H} \frac{m_i}{R_i}\left(1 - \frac{\mathcal{N}_{i}^{th}}{n_i}
       \right) \\
        & +  \sum_{a} \mathcal{B}_{ai}^{eff} \frac{\Gamma_a}{H}
       \frac{m_a}{R_a}\left(\frac{n_a}{n_i} - \frac{\mathcal{N}_{ai}^{th}}{n_i}
        \right) \nonumber
       \\
      R_i' & =  -3 \frac{P_i}{n_i} + \sum_{a} \mathcal{B}_{ai}^{eff}
      \frac{\Gamma_a}{H} m_a \left( \frac{1}{2} - \frac{R_i}{R_a} \right) \left(\frac{n_a}{n_i} -
      \frac{\mathcal{N}_{ai}^{th}}{n_i} \right)\end{aligned}

-  Coherent Oscillating fields:

   .. math::

      \begin{aligned}
      N_i' & = -3 - \frac{\Gamma_i}{H} \nonumber \\
      R_i' & = 0 \label{eq:COeq}\end{aligned}

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
   \frac{P_i}{n_i} & = T_i \nonumber \\
   \frac{\rho_i}{n_i} & = T_i \left[ \frac{K_1(m_i/T_i)}{K_2(m_i/T_i)} \frac{m_i}{T_i} + 3 \right] \label{eq:p1}\end{aligned}

where :math:`K_{1,2}` are the modified Bessel functions. In particular,
if :math:`m_i/T_i \gg 1`:

.. math:: \frac{\rho_i}{n_i} \simeq T_i \left[\frac{3}{2} + \frac{m_i}{T_i}  + 3 \right] \Rightarrow \frac{P_i}{n_i} = T_i = \frac{2 m_i}{3}\left( \frac{R_i}{m_i} -1 \right)

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

.. math:: T = \left(\frac{g_{*S}(T_R)}{g_{*S}(T)}\right)^{1/3} T_R \exp[N_S/3 -x] \Rightarrow \rho_R = \frac{\pi^2}{30} g_{*}(T) T^4

Eqs.([eq:Sfin])-([eq:COeq]), with the auxiliary equations for :math:`H`
(Eq.([H])) and :math:`P_i/n_i` (Eq.([Pfin])) form a set of closed
equations, which can be solved once the initial conditions for the
number density (:math:`n_i`), energy density (:math:`\rho_i`) and
entropy (:math:`S`) are given. For thermal fluids we assume:

.. math::

   \begin{aligned}
   n_i(T_R) & = \left\{ 
   \begin{array}{ll} 
   0 & , \mbox{ if ${\langle \sigma v \rangle}_i \bar{n}_i/H|_{T=T_R} < 10$} \\
   \bar{n}_i(T_R) & , \mbox{ if ${\langle \sigma v \rangle}_i \bar{n}_i/H|_{T=T_R} > 10$} 
   \end{array} \right. \label{ni0TP} \\
   \frac{\rho_i}{n_i}(T_R) & = \frac{\bar{\rho}_i}{\bar{n}_i}(T_R)\end{aligned}

where :math:`\bar{\rho}_i` is the equilibrium energy density (with zero
chemical potential) for the particle :math:`i`. While for coherent
oscillating fluids the initial condition is set at the beginning of
oscillations:

.. math::

   \begin{aligned}
   n_i(T^{osc}_i) & =\frac{\rho_i^{0}}{m_i(T^{osc}_i)} \\
   \frac{\rho_i}{n_i}(T^{osc}_i) & = m_i\end{aligned}

where :math:`T^{osc}_i` is the oscillation temperature, given by
:math:`3H(T^{osc}_i) = m_i(T^{osc}_i)` and :math:`\rho_i^{0}` the
initial energy density for oscillations.

Finally, the initial condition for the entropy :math:`S` is trivially
obtained, once we assume a radiation dominated universe at
:math:`T=T_R`:

.. math:: S(T_R) = \frac{2 \pi^2}{45} g_*(T_R) T_R^3 R_0^3
