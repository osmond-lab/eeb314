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

# Tutorial 10: Superinfection and viral co-existence

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

Recall the epidemiological model from Lecture 18, which we used to examine the evolution of virulence:

$\frac{\mathrm{d}S}{\mathrm{d}t} = \theta - dS - \beta(v)SI - \beta(v_m)SI_m$

$\frac{\mathrm{d}I}{\mathrm{d}t} = \beta(v)SI - (d+v)I$

$\frac{\mathrm{d}I_m}{\mathrm{d}t} = \beta(v_m)SI_m - (d+v_m)I_m.$

As a reminder, $S$ is the number of uninfected (susceptible) hosts and $I$ and $I_m$ are the numbers of hosts infected by resident and mutant strains of the virus. Susceptible hosts immigrate at rate $\theta$, die at per capita rate $d$, and are infected upon meeting infected hosts at rate $\beta$, which is a function of the strain's virulence, $v$, which elevates the death rate of infected hosts.

Recall also that, in the absence of mutants ($I_m = 0$), the resident equilbrium is $\hat{S} = \frac{d+v}{\beta(v)}$, $\hat{I} = \frac{\theta}{d+v}-\frac{d}{\beta(v)}$.

From here, we determined the invasion fitness of the mutant, $r(v_m,v) = -d - v_m + \hat{S}\beta(v_m)$, which implied that the mutant would invade and supplant the resident whenever 

$r(v,v_m) = -d - v_m +S^*\beta(v_m) > 0 \Rightarrow  \frac{\beta(v)}{d+v}<\frac{\beta(v_m)}{d+v_m}.$

I.e., evolution leads to the strain that maximizes $\frac{\beta(v)}{d+v}$.

In this model, the mutant and resident strains can coexist. However, in nature, we often see multiple co-existing strains within a host population, even strains with different $R_0$ values.

## Model

Consider a modification of the model above that allows for 'super-infection', by which the mutant pathogen can supersede resident infections:

$\frac{\mathrm{d}S}{\mathrm{d}t} = \theta - dS - \beta(v)SI-\beta(v_m)SI_m$

$\frac{\mathrm{d}I}{\mathrm{d}t} = \beta(v)SI - (d+v)I -\beta(v_m)II_m$

$\frac{\mathrm{d}I_m}{\mathrm{d}t} = \beta(v_m)SI_m + \beta(v_m)II_m - (d+v_m)I_m$.


## Problem

Find the resident equilbrium, Jacobian, and invasion fitness of the new system:

- When can the mutant invade?
- How does this criteria differ from the model without superinfection, $\frac{\beta(v)}{d+v}<\frac{\beta(v_m)}{d+v_m}$?

Now, see if the resident can invade a population of mutants. That is, set $I=0$, find the new mutant-only equilibrium $S$ and $I_m$, and evaluate the Jacobian at the new equilbrium to find the invasion fitness of $I$. 

- Are there parameter regimes where both 1) the mutant can invade the resident and 2) the resident can invade the mutant? You might show this by checking various combinations of parameters.

Bonus: Solve for the coexistence equilibrium (on paper or with SymPy). 


<pre data-executable="true" data-language="python">

</pre>
