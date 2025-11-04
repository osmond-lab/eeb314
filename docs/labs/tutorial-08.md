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

# Tutorial 8: Rosenzweig-MacArthur predator-prey model

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

## Model

The Rosenzweig-MacArthur predator-prey model assumes that the density of prey, $N$, grows logistically in the absence of predators, with per capita growth rate $r$ and carrying capacity $k$. Predators, with density $P$, consume prey at a per prey per predator rate of $a/(1+b N)$. The form of this per capita consumption rate implies that the rate of prey eaten per predator asymptotes to $a/b$ as $N$ gets large; essentially, we assume it takes some time for a predator to eat a prey, limiting the rate prey are consumed even when there are many of them. Finally, we assume each consumed prey is converted into a fraction $c$ new predators (i.e., consuming $1/c$ prey is required to make a new predator), and that predators die at per capita rate $d$. Putting all this together in a continuous-time model we have

$\frac{\mathrm{d}N}{\mathrm{d}t} = r N (1 - N / k) - \frac{a}{1 + b N} N P$

$\frac{\mathrm{d}P}{\mathrm{d}t} = c\frac{a}{1 + b N} N P - d P$.

Based on the model description we can assume all parameters are positive.

## Problem

Solve for the three equilibria. When are they biologically valid?

Determine the conditions for local stability for each equilibria, assuming they are biologically valid. (Hint: use the Routh-Hurwtiz criteria.)

Show that the equilibrium with both species present is unstable when the prey's carrying capacity is large enough, $k>\frac{ac+bd}{b(ac-bd)}$. (Hint: use the trace.) This is called the paradox of enrichment, because enriching the resources available for prey (increasing $k$) might be (naively) expected to lead to more prey and hence more stable predator-prey coexistence.
