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
    pairwise interactions, the Barnes-Hut(BH) algorithm and the Fast Multipole
    Method(FMM), along with a brute-force method have been implemented in `C++`.
    The time complexity of the algorithms and their scaling with different
    parameters have been studied. At large particle numbers, FMM was found to be
    the most efficient of the three, followed by BH, both of which showed
    close-to-linear scaling with particle number. The three methods were then
    used with a leapfrog integrator to simulate a number of physical scenarios,
    and their results were compared, showing good agreement in their qualitative
    behaviours.

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
implementation of two well-known hierarchical approaches, the Barnes-Hut(BH)
method and Fast Multipole Method(FMM), which differ in their treatment of
distant boxes.

In the next section, I will discuss mathematical prerequisites, after which the
algorithms themselves are presented. [Section @sec:res] provides detailed
complexity analyses and simulation results, and I discuss the code in [Section
@sec:discussions]. Finally a summary is given in [Section @sec:summary].

# Mathematical background {#sec:math_bg}

Before we can approach the algorithms, we must look at multipole expansions,
which lies at the heart of our implementation. For lack of space, the treatment
is necessarily very brief, and I refer interested readers to the appendix A of
[DEHNEN 2014]{.mark} for detailed mathematical derivations.

[CHANGE COM TO COC]{.mark}

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

--- 
# \begin{align}
# \label{eq:multipole}
# \psi\left(\bm{x_b} - \bm{x_a}\right)
# = \sum^p_{\left|n\right|=0}
# \sum^{p-{\left|n\right|}}_{\left|m\right| = 0}
# \frac{\left(-1\right)^{\left|m\right|}}{\bm{n}!\bm{m}!}
# \bm{r}_b^{\bm{n}} \bm{r}_a^{\bm{m}} \bm{\nabla}^{\bm{n} + \bm{m}}
# \psi\left(\bm{r}\right) + R_p
# \end{align}
---

It turns out that in spherical harmonics expansions, we only need two sets of
coefficients $M_n^m$ and $F_n^m$, $0 \le n \le p$, $-n \le m \le n$. $M_n^m$,
localised at the COC of clusters $a$ and $b$ to completely describe an expansion
of order $p$. Convenient expressions, known as _kernels_, exist for us to
calculate $M_n^m$ from $a$'s charge distribution, calculate $F_n^m$ from $M_n^m$
and force on $b$ particles from $F_n^m$. We can additionally also translate
$M_n^m$ and $F_n^m$ to new centres. All these mathematical operations have been
given convenient names as shown in [Table @tbl:kernels]. Because most of them
take a form similar to matrix multiplications between one set of coefficients
and the set of solid harmonics of order $p$, the operations are generally of
complexity $\mathcal{O}\left(p^4\right)$.

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

First proposed in [BARNES 1986]{.mark}, the Barnes-Hut method starts by
constructing an efficient representation of the particle distribution in space,
known as an octree, which stores subdivisions of space ("box") of increasingly
small volume. Through this data structure, by means of a experimenter-defined
criteria, it is very easy to determine when a test particle is "distant" from a
box. The Barnes-Hut method then leverages this distance by replacing all
particles in the distant box with an _effective image_, therefore bypassing the
need to calculate the interaction pair by pair. Below I will discuss the
algorithm in the aforementioned order.

### Octree construction

The octree is a tree structure representing 3D space. Each node represents a
cuboid in space and can have up to 8 children corresponding to the 8 cuboid
octants of the parent node. The root node represents the whole space of
interest. As we descend down the octree we get a _hierarchy of scales_,
therefore the name of this category of algorithms.

Stored on each node are the particles in the represented region, total mass,
centre of mass (COM) and multipole coefficients $M_n^m$, $F_n^m$ at the COM up
to order $p$ [^bhmultipole]. The storage of pre-processed information is key to
our speed-ups. A general application of the octree requires that it has
$n_\mathrm{max}$ particles in each leaf node only; For Barnes-Hut,
$n_\mathrm{max} = 1$.

[^bhmultipole]: The original Barnes-Hut paper stores only the total mass and so
    is in fact the zeroth order approximation. I extended this to work for
    higher order expansions.

Such a tree can be constructed with the recursive algorithm in [Listing
@lst:buildtree], starting from root. This can be done in time
$\mathcal{O}\left(n \log{n}\right)$, since similar to a binary tree, we expect a
mean tree height of order $\log{n}$, and we are inserting a total of $n$
particles. Nodes are created whenever we visit a previously not-existing one.

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

    update COM, total mass on NODE
```

: Octree construction. {#lst:buildtree}


So far the $M_n^m$ coefficients remain uninitialised; For this we perform a
depth-first search (DFS), calculating $M_n^m$ at the leaf nodes first with the
P2M kernel, then on returning to the parent level, calculate the parent $M_n^m$
with the M2M kernel passing information from the children. Note that this cannot
be done simultaneously with the `AddParticle` routine, since the COM of each
node, at which the $M_n^m$ is calculated is not known then.

### Barnes-Hut force calculation

With the octree ready, we are now able to calculate the forces. Introducing the
parameter $\theta$, we can define a _multipole acceptance criteria_(MAC):

> **Multipole Acceptance Criteria** We are allowed to multipole-expand a source
> box if $d / x < \theta$, where $d$ is the maximum of the three sides of a box,
> and $x$ is the distance from the particle to the box COM.

The Barnes-Hut calculation on one particle is then shown in [Listing
@bh_interact]

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

The first step of FMM, octree construction is nearly identical to that in
Barnes-Hut, with the exception that we relax the leaf node particle limit to a
small constant $n_\mathrm{max}\sim\mathcal{O}\left(1\right)$. This doesn't
affect the overall complexity, for the complexity of each step in [Listing
@lst:buildtree] is merely modified by a constant.

The essential idea of an FMM calculation is the approximation of distant
box-to-box interactions with multipole expansions. This is made possible by
using the M2L kernels, an operation independent of particle numbers. For us to
decide when boxes are well-separated, we again need a multipole acceptance
criteria (MAC):

> **Multipole Acceptance Criteria** Two boxes are considered well-separated if
> the sum of the maximum sides of the boxes $l$ and the distance between box
> COMs $d$ satisfies $l / d < \theta$ where $\theta$ is a constant called the
> _opening angle_.

The dual-tree traversal style implementation in my code is then summarised as
follows, starting at the root node:

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

: Core of FMM force calculation. {#lst:fmm_interact}

To understand the effect and complexity of [Listing @lst:fmm_interact], consider
an arbitrary leaf node $\alpha$ and the chain of nodes leading up to the root.
As `Interact` descend from root, large distant boxes will imprint their effect
on an ancestor of $\alpha$. For $\alpha$ then, the descendants of those distant
boxes are now irrelevant. moving one level down, yet some more nodes now satisfy
the MAC so interact with the new $\alpha$ ancestor via M2L. Importantly, the
number of nodes that interact at each level is **roughly a constant**, since the
vast number of more distant small boxes have all been taken care of at the
previous levels with bigger boxes. This continues until we hit nodes with which
it becomes necessary or economical (number of pairs below a limit
$n_\mathrm{direct}$)to perform pairwise interactions. 

All the necessary interactions to $\alpha$ are therefore imprinted in its
ancestor's $F_n^m$ coefficients. To calculate the total force on each particle
of $\alpha$, we descend from root down to $\alpha$, shift and add $F_n^m$ to the
child node COM at each step with the L2L kernel, and finally apply the L2P
kernel on particles in $\alpha$. This can be done for all leaf nodes in one go,
with a simple push-down-then-descend DFS.

In building the tree, P2M kernel calculates $M_n^m$ coefficients on the root
node in time $\mathcal{O}\left(n\right)$. It is then passed upward in another
$\mathcal{O}\left(n\right)$ operations. The same is true the downward pass of
$F_n^m$ via L2L and the conversion of $F_n^m$ to forces. The intermediate
`Interact` function is also expected to be of order $\mathcal{O}\left(n\right)$
or even less, as demonstrated in [Dehnen 2002]{.mark}. `Interact` will however
still dominate the execution time, for it calls the time-consuming M2L kernel
much more than the other parts call their kernels.

# Tests and simulations

## Complexity and error tests

Various tests have been conducted where I change one parameter shared between
Barnes-Hut and FMM, while other parameters are fixed with values specified in
[Table @tbl:def_values]. For each parameter, $10$ random uniform distributions
of charges are generated from different (but known) seeds as repeats, then fed
to brute force, Barnes-Hut and FMM. They are then timed by `std::chrono` calls
which sandwich the computational code. Finally the timing results and
information on all particles are saved to files, to be processed by `Python`
scripts. We define the error of a hierarchical method to be the relative error
against brute-force results.

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
$\theta$. The low-$\theta$ limits in error is likely floating point error due to
different orders of calculation between BH, FMM and brute-force.

<div id="fig:changing_theta">
![Time against $\theta$](time_theta.pdf){#fig:time_theta width=45%}
![Error against $\theta$](err_theta.pdf){#fig:err_theta width=45%}

Response to variation in $\theta$.
</div>

### Effect of changing $p$

<div id="fig:changing_p">
![Time against $p$](time_p.pdf){#fig:time_p width=45%}
![Error against $p$](err_p.pdf){#fig:err_p width=45%}

Response to variation in $p$.
</div>

### Error distributions
