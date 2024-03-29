---
title: Hierarchical Methods for N-body Simulations
author: Jimmy Chen
date: 27 March 2024
geometry: margin=4cm
header-includes: |
    \usepackage{bm}
codeBlockCaptions: true
colorlinks: true
linkcolor: blue
subfigGrid: true

abstract: |

    Two classic algorithms for force calculations in N-body simulations with
    pairwise interactions, the Barnes-Hut (BH) algorithm and the Fast Multipole
    Method (FMM), along with a brute-force method have been implemented in
    `C++`. The time complexity of the algorithms and their scaling with
    different parameters have been studied. At large particle numbers, FMM was
    found to be the most efficient of the three, followed by BH, both of which
    showed close-to-linear scaling with particle number. The three methods were
    then used with a leapfrog integrator to simulate a number of physical
    scenarios, and their results were compared, showing good agreement in their
    qualitative behaviours.

---

# Background

N-body simulations have been a great tool to astrophysicists for the study of
galactic and celestial dynamics, and as the problems of interest grow
increasingly large in scale, the calculations are also becoming more
computationally demanding. In an age when Millennium Runs regularly simulate
over 10 billion particles [CITATION]{.mark}, execution efficiency of the program
is essential, and the naive approach of complexity $\mathcal{O}\left(n^2\right)$
proves completely untenable. Indeed, much work has been done to improve the
efficiency of force calculation, the typical bottleneck in the simulation. Many
such algorithms have been proposed, typically trading small amounts of error
to significantly boost performance.

One category of such algorithms is known as hierarchical methods, which involves
recursively subdividing space into boxes of octants, and applying special
treatment to particles of well-separated boxes. In this paper, I present
implementation of two well-known hierarchical approaches, the Barnes-Hut (BH)
method and Fast Multipole Method (FMM), which differ in their treatment of
distant boxes.

In the next section, I will discuss mathematical prerequisites, after which the
algorithms themselves are presented. [Section @sec:test_and_sim] provides
detailed complexity analyses and simulation results, and I discuss optimization
efforts in [Section @sec:optimizations]. Finally a summary is given in [Section
@sec:conclusions].

# Mathematical background {#sec:math_bg}

Before we can approach the algorithms, we must look at multipole expansions,
which lies at the heart of our implementation. For lack of space, the treatment
is necessarily very brief, and I refer interested readers to the appendix A of
[DEHNEN 2014]{.mark} for detailed mathematical derivations.

Our aim is to find for an inverse-square law of coupling constant $G$, an
expression for the potential $\psi$ at one cluster of charges $b$, with centre
of charge (COC) $r_b$, due to a far-away, well separated cluster $a$ with COC
$r_a$. For a particle $j$ in $b$ at $\bm{r}^\prime = \bm{r}_b + \bm{r}_j$ with
charge $q_j$, the potential of a particle $i$ in $a$ at $\bm{r} = \bm{r}_a +
\bm{r}_i$, charge $q_i$ is written
\ 
\begin{align}
\psi\left(\bm{r}^\prime\right) =
G \frac{q_i}{\left|\bm{r}_j + \bm{r}_b - \left(\bm{r}_i + \bm{r}_a\right)\right|}
\end{align}

Note that the denominator is in fact $\left(\bm{r}_b - \bm{r}_a\right) +
\left(\bm{r}_j - \bm{r}_i\right)$, where, since the clusters are well separated,
the first term in brackets dominates. This invites us to perform an Taylor
expansion in the small parameter (second term) and truncate at the $p$-th order,
the result of which is a _multipole expansion of order $p$_. So far we have
worked free of a basis: The resulting expression can be written in Cartesian
form, or in a more convenient form involving solid spherical harmonics, which I
adopt in my code.

It turns out that in spherical harmonics expansions, we only need two sets of
coefficients $M_n^m$ and $F_n^m$, $0 \le n \le p$, $-n \le m \le n$ localised at
the COC of clusters $a$ and $b$ to completely describe an expansion of order
$p$. Convenient expressions, known as _kernels_, exist for us to calculate
$M_n^m$ from $a$'s charge distribution, calculate $F_n^m$ from $M_n^m$ and force
on $b$ particles from $F_n^m$. We can additionally also translate $M_n^m$ and
$F_n^m$ to new centres. All these mathematical operations have been given
convenient names as shown in [Table @tbl:kernels]. Because most of them take a
form similar to matrix multiplications between one set of coefficients and the
set of solid harmonics of order $p$, the operations are generally of complexity
$\mathcal{O}\left(p^4\right)$.

| Name  | Purpose |
|------:|:--------|
|  P2M  | Find $M_n^m$ from local charge distribution |
|  M2M  | Shift $M_n^m$ expansion centre  |
|  M2L  | Find $F_n^m$ from distant $M_n^m$  |
|  M2P  | Find $\psi$, $\bm{g}$ on particle from distant $M_n^m$  |
|  L2L  | Shift $F_n^m$ expansion centre   |
|  L2P  | Find $\psi$, $\bm{g}$ on charge from local $F_n^m$  |

: Kernels in the multipole expansion. {#tbl:kernels}

# Algorithmic approach {#sec:algo}

## Barnes-Hut method

The Barnes-Hut approach starts by constructing an "octree" data structure, with
which it is easy to determine when a particle is well-separated from a region in
space. This allows us to make approximations by replacing the charges in that
space with a multipole expansion from which forces are derived efficiently.

### Octree construction

The octree is a tree structure representing 3D space. Each node represents a
cuboid in space and can have up to 8 children corresponding to the 8 cuboid
octants of the parent node. The root node represents the whole space of
interest.

We can store on each node the particles in the represented region, total mass,
COC and multipole coefficients $M_n^m$, $F_n^m$ at the COC up to order $p$
[^bhmultipole]. The storage of pre-processed information is key to our
speed-ups. A general application of the octree requires that each leaf node
holds no more than $n_\mathrm{max}$ particles; For Barnes-Hut, $n_\mathrm{max} =
1$.

[^bhmultipole]: The original Barnes-Hut paper stores only the total mass and so
    is in fact the zeroth order approximation. I extended this to work for
    higher order expansions.

Such a tree can be constructed with the recursive algorithm in [Listing
@lst:buildtree], starting from root. This can be done in time
$\mathcal{O}\left(n \log{n}\right)$, since the expected mean tree height is of
order $\log{n}$, and this is repeated for each particle. Nodes are created
whenever we visit a previously not-existing one.

```{.python}
def AddParticle(NODE, PARTICLE):
    if NODE is a leaf:
        if NODE has > n_max particles:
            split NODE into smaller octants
            AddParticle(correct child_node, PARTICLE and particles in NODE)
        else:
            add PARTICLE into NODE

    if NODE is not a leaf:
        AddParticle(correct child_node, PARTICLE)

    update COC, total mass on NODE
```

: Octree construction. {#lst:buildtree}


To initialise the $M_n^m$ coefficients, we then perform a depth-first search
(DFS), calculating $M_n^m$ at the leaf nodes first with the P2M kernel, then on
returning to the parent level, calculate the parent $M_n^m$ with the M2M kernel
passing information from the children.

### Barnes-Hut force calculation {#sec:bh_force}

With the octree ready, we are now able to calculate the forces. Introducing the
parameter $\theta$, we can define a _multipole acceptance criteria_(MAC):

> **Multipole Acceptance Criteria** We are allowed to multipole-expand a source
> box if $d / x < \theta$, where $d$ is the maximum side of the box, and $x$ is
> the distance from the particle to the box COC.

The Barnes-Hut calculation on one particle is then shown in [Listing
@lst:bh_interact]. We simply repeat this procedure to get forces on each
particle.

```{.python}
def Interact(NODE, PARTICLE):
    if NODE is a leaf:
        add force by single particle in NODE to PARTICLE
    if NODE meets MAC:
        calculate force by NODE on PARTICLE via M2P kernel
    else:
        Interact(child_nodes, PARTICLE)
```

: Barnes-Hut force calculation. {#lst:bh_interact}

The algorithm essentially simplifies calculation of force from distant boxes, by
replacing all those pairwise interactions with multipole expansion results. An
excellent argument in [BARNES 1986]{.mark} demonstrates that this is also of
complexity $\mathcal{O}(n \log{n})$, leading to an overall complexity of
$\mathcal{O}\left(n\log{n}\right)$. There is an additional dependence on
expansion order $p$, bounded above by $\mathcal{O}\left(p^4\right)$, since not
all interactions are calculated from the multipole expansion.


## Fast multipole method

We again build the octree at the start of FMM, this time setting $n_\mathrm{max}
\sim \mathcal{O}\left(1\right)$, which modifies the complexity of each step in
[Listing @lst:buildtree], and therefore the overall complexity, only by a
constant.

At the heart of FMM is the approximation of well-separated box-to-box
interactions with multipole expansions. This is made possible by using the M2L
kernel, an operation independent of particle numbers. We again need a multipole
acceptance criteria, this time between two boxes:

> **Multipole Acceptance Criteria** Two boxes are considered well-separated if
> the sum of the maximum sides of the boxes $l$ and the distance between box
> COCs $d$ satisfies $l / d < \theta$ where $\theta$ is a constant called the
> _opening angle_.

We can then make use of the M2L kernel with the algorithm in [Listing
@lst:fmm_interact].

```{.python}
def Interact(NODE1, NODE2):
    if NODE1 and NODE2 are leaf nodes or have few interaction pairs:
        Calculate pairwise forces directly
    if NODE1 and NODE2 satisfy MAC:
        Calculate F coefficients on both nodes via M2L kernels
    else:
        if NODE1 is NODE2:
            for all combinations of NODE1 child_nodes:
                # including child_node_1 == child_node_2
                Interact(child_node_1, child_node_2)
        else:
            for all child_nodes of larger NODE:
                Interact(child_node, smaller NODE)
```

: Dual-tree traversal {#lst:fmm_interact}

To understand the effect of [Listing @lst:fmm_interact], consider the chain of
nodes leading from the root to an arbitrary leaf node, $\{a_0, a_1 ..., a_n\}$.
As `Interact` descend to node $a_i$, large distant boxes $\{b_i\}$ will imprint
their influence in the $F_n^m$ of $a_i$. For descendants of $a_i$ then, the
children of $b_i$ are now irrelevant. Moving down to level $a_{i + 1}$, yet some
more nodes now satisfy the MAC and will be multipole-expanded. The number of new
nodes however should be _roughly constant_: only boxes sandwiched between nearby
boxes that don't satisfy MAC, and $\{b_i\}$, can be considered. We continue
descending, until it becomes necessary (both nodes are leaves) or more
economical (number of direct interactions below $n_\mathrm{direct}$) to perform
brute-force directly.

All the necessary interactions to $a_n$ are therefore imprinted in its
ancestors' $F_n^m$ coefficients. To calculate the total force on each particle
of $a_n$, we descend from $a_0$ to $a_n$, shift and add $F_n^m$ to the
child node COC at each step with the L2L kernel, and finally apply the L2P
kernel on particles in $a_n$. This can be done for all leaf nodes in one go
with a simple push-down-then-descend DFS after `Interact` has finished.

In building the tree, P2M kernel calculates $M_n^m$ coefficients on the root
node in time $\mathcal{O}\left(n\right)$. It is then passed upward in another
$\mathcal{O}\left(n\right)$ operations. The same is true the downward pass of
$F_n^m$ via L2L and the conversion of $F_n^m$ to forces. The intermediate
`Interact` function is also expected to be of order $\mathcal{O}\left(n\right)$
or even less, as demonstrated in [Dehnen 2002]{.mark}. `Interact` will however
still dominate the execution time, for it calls the time-consuming M2L kernel
much more than the other parts call their kernels.

# Tests and simulations {#sec:test_and_sim}

## Complexity and error tests

Various tests have been conducted where I change one parameter shared between
Barnes-Hut and FMM, while other parameters are fixed with values specified in
[Table @tbl:def_values]. For each parameter, $10$ random uniform distributions
of charges are generated from different (but known) seeds as repeats, then fed
to brute-force, Barnes-Hut and FMM. They are then timed by `std::chrono` calls
which sandwich the computational code. Finally the timing results and
information on all particles are saved to files, to be processed by `Python`
scripts. We define the error of a hierarchical method to be the relative error
against brute-force results, and use the error averaged over all three
components and over the $10$ runs in our figures below.

All tests are conducted on a MacBook Pro early 2015 with an Intel
i7-5557U@3.1GHz and 16GB RAM. The system is a fully upgraded Arch Linux as of
26th March with Linux kernel 6.8.1. All code has been compiled with `clang++`
version 17.0.6 and the following optimization flags:

> -std=c++17 -O3 -DNDEBUG -ftree-vectorize -flto -finline-small-functions
> -march=native

For best stability, Turbo Boost has been disabled, and CPU governor has been set
to `performance`, which disables dynamic scaling and fix the CPU frequency to
its highest (3.1GHz). The laptop is powered throughout, and all foreground
processes are closed. No throttling was observed during the runs.


| Variable  | Default |
|------:|:--------|
|  $n$       | $1000$ |
|  $\theta$  | $0.5$ |
|  $p$    |  $3$ |  
|  $n_\mathrm{direct}$ (FMM) | $3$ |
|  $n_\mathrm{max}$ (FMM) | $5$ |

: Default variable values. {#tbl:def_values}

### Effect of changing $n$

As expected, both Barnes-Hut and FMM show close-to-linear scaling against number
of masses, demonstrating their superiority at large $n$. Code profiling with
`valgrind` demonstrates that both code spends the most time in their `Interact`
routines ($97\%$ for Barnes-Hut, $82\%$ for FMM).

The Barnes-Hut method has a much larger constant of proportionality -
This is due to the need to compute multipole expansions for each sink particle,
which is very demanding. FMM does not face the same problem, due to its use of
box-to-box approximations.

At large $n$, the mean error is largely constant. the very small error of FMM at
small $n$ is likely because most interactions are allowed to proceed through
pairwise interactions.

<div id="fig:changing_n">
![Time against $n$](time_n.pdf){width=45%}
![Error against $n$](err_n.pdf){width=45%}

Response to variation in $n$.
</div>

### Effect of changing $\theta$

$\theta$ sets what we consider to be well-separated boxes in Barnes-Hut and FMM.
As we reduce $\theta$ from $1$ in [Figure @fig:time_theta], there is first a
rise in computation time, caused by FMM and BH having to do more multipole
expansions at lower levels of the octree, followed then by a drop, corresponding
to the reduction to brute-force (brute-force is still quick at $n = 1000$,
therefore the _speed-up_.) At the small $\theta$ limit, both methods have a
constant overhead above brute-force. Barnes-Hut is significantly slower, since
it needs to traverse the entire octree for each particle, whereas FMM benefit
from the one-shot dual-tree traversal and a shallower tree thanks to
$n_\mathrm{max} \ge 1$.

Barnes-Hut displays a slower drop in error: This is likely because its MAC,
which only concerns one cell, can allow multipole expansions at smaller
$\theta$. The low-$\theta$ limits in error is likely floating point errors in
the computation due to different calculation orders between BH, FMM and
brute-force.

<div id="fig:changing_theta">
![Time against $\theta$](time_theta.pdf){#fig:time_theta width=45%}
![Error against $\theta$](err_theta.pdf){#fig:err_theta width=45%}

Response to variation in $\theta$.
</div>

### Effect of changing $p$

As explained in [Section @sec:bh_force], we expect the complexity dependence $p$
to be bounded by $\mathcal{O}\left(p^4\right)$, but lower in practice. This is
indeed seen in the timing results of [Figure @fig:changing_p]. FMM has a
shallower power law, likely because its conditions for direct interaction is
more permissive.

<div id="fig:changing_p">
![Time against $p$](time_p.pdf){#fig:time_p width=45%}
![Error against $p$](err_p.pdf){#fig:err_p width=45%}

Response to variation in $p$.
</div>

As for the error, multipole expansions should converge quickly with reasonable
MACs. This is the case for both, but FMM is more sensitive to the value of $p$,
because it employs five multipole kernels in one calculation compared to three
in BH. Another interesting feature is that error distributions in BH for $p = 1$
and $p = 2$ are almost identical. This is not fully understood, but inspection
shows second order terms have been suppressed, possibly due to a mathematical
coincidence.

### Error distributions

For each combination of parameters, we have taken their error distributions. We
only show six distributions from BH with varying $\theta$, as most distributions
share the same features. On a logarithmic scale, tightening the constraint
(reducing $\theta$ here) causes the distribution to uniformly shift to the left,
but the same shape is maintained, until $\theta$ is so small the MAC simply
can't be met. This is when BH transitions to the brute-force regime ([Figures
@fig:bh_to_brute_1] to [@fig:bh_to_brute_3]).

Another feature is the long tail towards the right: in each case a small portion
of particles suffer from significant error. As noted in [DEHNEN 2014]{.mark},
this can be remediated by using a more sophisticated MAC.

<div id="fig:theta_dist">
![$\theta = 1.0$](dist/bh_hist_1.pdf){width=45%}
![$\theta = 0.316$](dist/bh_hist_0.316.pdf){width=45%}

![$\theta = 0.1$](dist/bh_hist_0.1.pdf){width=45%}
![$\theta = 0.0237$](dist/bh_hist_0.0237.pdf){width=45% #fig:bh_to_brute_1}

![$\theta = 0.0178$](dist/bh_hist_0.0178.pdf){width=45%}
![$\theta = 0.01$](dist/bh_hist_0.01.pdf){width=45% #fig:bh_to_brute_3}

Distribution of errors as $\theta$ is reduced.
</div>

## Simulation tests {#sec:sim_tests}

To verify the reliability of my algorithms, I have run simulations under three
initial conditions and compared BH, FMM results with brute force. After each
round of force calculation, a leap frog integrator evolves the system forward
with step size $\delta t = 0.01$, for a total time of $\Delta t = 10$. Tests are
run under the same conditions and default parameters as in the last section. The
three conditions are:

1. Cold start

   Generate a random uniform distribution of stationary charges in a sphere, and
   assign masses (and gravitational charge) to each with a normal distribution.

2. Single "galaxy"
    
   Generate a random uniform distribution of charges in a disk about set axis,
   excluding a central radius. Add a "supermassive black hole (SMBH)" - particle
   of very high mass at the centre, and set all particle velocities to that of
   circular motion around the SMBH.
\filbreak
3. Two "galaxies"

   Generate two galaxies as above, tilted at different angles, then boost them
   with velocities required for the SMBHs at the centre to orbit around each
   other.

Due to lack of space we present and discuss the most interesting third case only
here, while some snapshots of the other two can be found in the
[APPENDIX]{.mark}. Videos of full simulations are also available.

Snapshots of the third system, shown in [Figure @fig:two_galaxies], display very
good qualitative agreement across the three methods. As the galaxies rotate
around each other, many outer stars do not feel enough attraction to keep
following the galaxies' orbit, and are scattered outwards. This is not seen in
the simulation of a single galaxy, and is thus a unique feature of this model.
The stars at the centre on the other hand feel similar attractions from both
galaxies, and can be exchanged between the them. Interestingly, the galactic
centres seem further apart at half-period, likely because the attraction is
reduced when they are facing away from each other. At a full period the centres
are observed to return to almost the same initial locations.

<div id="fig:two_galaxies">
![](two-galaxies/brute_0.0.pdf){width=23%}
![](two-galaxies/brute_1.0.pdf){width=23%}
![](two-galaxies/brute_2.0.pdf){width=23%}
![](two-galaxies/brute_3.0.pdf){width=23%}

![](two-galaxies/bh_0.0.pdf){width=23%}
![](two-galaxies/bh_1.0.pdf){width=23%}
![](two-galaxies/bh_2.0.pdf){width=23%}
![](two-galaxies/bh_3.0.pdf){width=23%}

![](two-galaxies/fmm_0.0.pdf){width=23%}
![](two-galaxies/fmm_1.0.pdf){width=23%}
![](two-galaxies/fmm_2.0.pdf){width=23%}
![](two-galaxies/fmm_3.0.pdf){width=23%}

Simulation of two galaxies with total $n = 2000$. Top to bottom: brute-force,
BH, FMM. Left to right: $t = 0, 1, 2, 3$.
</div>

We can further inspect the change in energy and angular momentum against time
([Figure @fig:galaxies_E_L]). Several "steps" are seen as energy conservation is
severely violated, but these are likely caused by the leap frog integrator
evolving a particle with very large acceleration. In fact, in between every
jump, the energy has been reasonably well conserved for the two hierarchical
methods, demonstrating their sufficiency. In order to prevent the sharp "steps"
however, gravity softening will be needed, which demands modification of the
algorithms.

On the other hand, angular momentum is much better conserved. Note an
interesting "synchronised", "smooth" pattern in FMM's angular momenta
components, which is likely correlated with the relative positions of the
galactic centres, and which might even be understood with a perturbative
treatment.

<div id="fig:galaxies_E_L">
![](two-galaxies/E_t.pdf){width=45%}
![](two-galaxies/L_x_t.pdf){width=45%}

![](two-galaxies/L_y_t.pdf){width=45%}
![](two-galaxies/L_z_t.pdf){width=45%}

Change of energy and angular momentum as simulation progresses.
</div>

# Optimizations {#sec:optimizations}

Much effort has been made to speed up the running of the code, with the aid of
profiling tools such as `valgrind`. Here I briefly outline some techniques used,
and note some improvements that can be made later on.

**`-O3` and `-ftree-vectorize`**. The two compiler flags turn on the most
aggressive optimizations, and auto-vectorization of our program, which greatly
helps with loops encountered in vector operations. The combined effect is a near
two-fold increase compared to `-O2` optimizations.

**Memory allocation**. Because of stack limits, the `grid` carrying a static
list of all particles (containing many 3D vectors) has to be allocated on the
heap. However, all other particle instances are stored on the stack, reducing
the allocation and access times. This leads to a 10-fold speed-up compared to
when all vectors are dynamically allocated. At the end of each timestep, the
previous `grid` is written to disk and deallocated to save memory.

**`soul` labels**. To track the particles within a node, instead of storing
copies of the particle we merely store a label with which we can find the
particle in a separate list. This at least halves the memory requirement.

**Multipole kernels**. My spherical harmonics implementation of complexity
$\mathcal{O}\left(p^4\right)$ already beats the naive Cartesian expansion of
order $\mathcal{O}\left(p^6\right)$. However, as covered by [DEHNEN
2014]{.mark}, the complexity is further reduced to $\mathcal{O}\left(p^3\right)$
if we rotate the vector connecting the COC onto the z-axis. This might not be
better for computations of my size however, because of the extra overhead of
rotation operations.

**Recurrence relations**. We can pre-process all solid harmonics of order $p$ by
harnessing their recurrence relations in time $\mathcal{O}\left(p^2\right)$ with
simple math operations, instead of directly evaluating their computationally
costly expressions based on Legendre polynomials in the 4th order loop of the
kernels. This change brought a 10x speed-up to BH and FMM.

# Conclusions {#sec:conclusions}

In this paper we present implementations of two hierarchical algorithms,
Barnes-Hut and Fast Multipole Method for the accelerated calculation of N-body
inverse-square law interactions. Their time complexities, in particular, the
close-to-linear dependence on number of particles, were shown to be consistent
with theoretical descriptions. Much effort is put into optimization, allowing
the program to calculate forces on $10^5$ particles in $10 \,\mathrm{s}$ with
FMM, where it would have taken brute-force $120 \,\mathrm{s}$.

The algorithms have additionally been used to perform simulations on simplified
models such as collision of two galaxies, and good agreements were seen between
the hierarchical approaches and a brute-force calculation.

As a next step, energy and angular momentum conservation can be improved by
using higher order integrators and incorporating softened gravity, and we can
extend the range of applicability by building in e.g. the Stokes and Helmholtz
kernels.

\appendix

# Appendix

## Simulation results

Below we list partial results to the cold start and single galaxy initial
conditions. [TODO]{.mark} Side-by-side comparisons of results from the three
algorithms, for all three initial conditions, are included in the submission.

### Cold start

The initial condition for the cold start has been described in [Section
@sec:sim_tests]. In this run, $n = 500$. Snapshots are shown in [Figure
@fig:cold_start]; energy and angular momentum change, in [Figure @fig:cold_E_L].

<div id="fig:cold_start">
![](cold/brute_0.0.pdf){width=23%}
![](cold/brute_1.0.pdf){width=23%}
![](cold/brute_2.0.pdf){width=23%}
![](cold/brute_3.0.pdf){width=23%}

![](cold/bh_0.0.pdf){width=23%}
![](cold/bh_1.0.pdf){width=23%}
![](cold/bh_2.0.pdf){width=23%}
![](cold/bh_3.0.pdf){width=23%}

![](cold/fmm_0.0.pdf){width=23%}
![](cold/fmm_1.0.pdf){width=23%}
![](cold/fmm_2.0.pdf){width=23%}
![](cold/fmm_3.0.pdf){width=23%}

Simulation of cold start with total $n = 500$. Top to bottom: brute-force, BH,
FMM. Left to right: $t = 0, 1, 2, 3$.
</div>

<div id="fig:cold_E_L">
![](cold/E_t.pdf){width=45%}
![](cold/L_x_t.pdf){width=45%}

![](cold/L_y_t.pdf){width=45%}
![](cold/L_z_t.pdf){width=45%}

Change of energy and angular momentum in cold start.
</div>

### Single galaxy

For this we use $n = 1000$. Snapshots are shown in [Figure @fig:galaxy]. All the
masses away from centre are in fact circulating the centre, and despite the
initial condition of a uniform angular velocity, after some time the angular
velocity of outer masses drop, agreeing with galactic rotation curves. These
details are not at all obvious from snapshots alone, and the reader is therefore
encouraged to check out the actual animations.

For energy and angular momentum change, see [Figure @fig:galaxy_E_L].

<div id="fig:galaxy">
![](galaxy/brute_0.0.pdf){width=23%}
![](galaxy/brute_1.0.pdf){width=23%}
![](galaxy/brute_2.0.pdf){width=23%}
![](galaxy/brute_3.0.pdf){width=23%}

![](galaxy/bh_0.0.pdf){width=23%}
![](galaxy/bh_1.0.pdf){width=23%}
![](galaxy/bh_2.0.pdf){width=23%}
![](galaxy/bh_3.0.pdf){width=23%}

![](galaxy/fmm_0.0.pdf){width=23%}
![](galaxy/fmm_1.0.pdf){width=23%}
![](galaxy/fmm_2.0.pdf){width=23%}
![](galaxy/fmm_3.0.pdf){width=23%}

Simulation of a single galaxy with total $n = 1000$. Top to bottom: brute-force,
BH, FMM. Left to right: $t = 0, 1, 2, 3$.
</div>

<div id="fig:galaxy_E_L">
![](galaxy/E_t.pdf){width=45%}
![](galaxy/L_x_t.pdf){width=45%}

![](galaxy/L_y_t.pdf){width=45%}
![](galaxy/L_z_t.pdf){width=45%}

Change of energy and angular momentum in a single galaxy.
</div>

