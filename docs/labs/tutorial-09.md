<script type="text/x-thebe-config">
  {
      requestKernel: true,
      mountActivateWidget: true,
      mountStatusWidget: true,
      binderOptions: {
      repo: "mmosmond/executable-cells",
      ref: "main",
      },
  }
</script>
<script src="https://unpkg.com/thebe@latest/lib/index.js"></script>
<link rel="stylesheet" href="https://unpkg.com/thebe@latest/lib/thebe.css">

# Tutorial 9: Evolutionary diversification of competitors I

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

When does competition for resources favor a single dominant *evolutionarily stable strategy* that wins over all others? And when does it favor evolutionary diversification, allowing coexistence? Here we'll use evolutionary invasion analysis to address this question.

## Model

As a biological example, consider a population of zooplankton, such as daphnia, that feeds on algal species of different sizes. For simplicity we'll assume that each daphnid has a feeding strategy characterized by the alga that it is best able to consume. We'll also assume that being able to feed well on large algae comes at the expense of being able to feed on small algae (i.e., there is a *trade-off* between being able to feed on alga types).

Consider a resident, with population size $n$ and strategy $s$, that competes with a mutant, with population size $n_m$ and strategy $s_m$, via Lotka-Volterra competition in discrete time,

$\displaystyle{n(t+1) = n(t) \left(1 + r \left(1 - \frac{n(t) + \alpha(s,s_m) n_m(t)}{K(s)}\right)\right)}$

$\displaystyle{n_m(t+1) = n_m(t) \left(1 + r \left(1 - \frac{n_m(t) + \alpha(s_m,s) n(t)}{K(s_m)}\right)\right)}$.

where carrying capacities $K$ are determined by the strategy (e.g., some algae sizes are more abundant) and $\alpha(s,s_m)$ represents how competition between residents with algae preference $s$ and mutants with algae preference $s_m$ reduces the growth of residents, and vice versa for $\alpha(s_m,s)$.

## Problem

Solve for the non-zero equilibrium population size of a resident with trait value $s$ in the absence of the mutant.

Construct the Jacobian and evaluate it at this resident equilibrium. 

Find the two eigenvalues:

- What do we need to be true about the resident for stability (in the complete absence of the mutant)? Assume this is true.

- What is invasion fitness, $\lambda(s_m,s)$? 
