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

# Lecture 16: Evolutionary invasion analysis I

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

## Lecture overview

1. [Invasion fitness](#section1)
2. [Evolutionarily singular strategies](#section2)
3. [Summary](#section3)

In the models we've discussed so far we've taken the parameters to be fixed. In reality, many of these parameters can evolve. For example, in our model of exponential growth in discrete time, $n(t+1)=n(t) R$, we took $R$ to be the same for all individuals for all time. But any mutation causing a larger $R$ would increase in frequency, causing the value of $R$ to increase over time. In this lecture we'll explore how to determine the direction of evolution and the stability of evolutionary endpoints for more complex models using a technique called <b>evolutionary invasion analysis</b> (also known as adaptive dynamics).

The idea:

  - determine which parameter(s) of our model an evolving trait affects
  - take the population to be fixed for some "resident" trait value
  - determine the equilibria and stability of the system with only the resident trait
  - derive the growth rate of a new individual with a "mutant" trait value
  - ask when the mutant trait value will invade
  - look for potential evolutionary endpoints
  - determine the stability of those endpoints (next lecture)

!!! note "The evolution of dispersal"

    To motivate evolutionary invasion analysis, let's consider a model for the evolution of dispersal. 
    
    Imagine there are $S$ sites, with a most one individual reproducing at each. We census the population at the time of reproduction. A reproducing individual has a large number of offspring, $B$, and then dies. A fraction $d$ of those offspring disperse and a fraction $1-c$ of those survive. The survivors then equally divided among all sites. One individual in each site is then chosen at random to reproduce, which begins the life-cycle anew.
    
    The question is, how should dispersal, $d$, evolve? There is a cost, $c$, which selects against dispersal but dispersal also allows offspring to avoid competiting with their kin, which could select for more dispersal. We'll use evolutionary invasion analysis to sort this out.

<span id='section1'></span>
## 1. Invasion fitness
<hr>

Let's think about this analysis very generally (in discrete time).

Let the number of individuals with the resident trait value be $n$ and the number of individuals with the mutant trait value $n_m$. We'll assume asexual reproduction for simplicity, so that residents only produce residents and mutant only produce mutants. Let the potentially nonlinear dynamics of these two groups of individuals depend on their respective trait values, $z$ and $z_m$,

$$
\begin{aligned}
n(t+1) &= n(t) R(n(t), n_m(t), z, z_m)\\
n_m(t+1) &= n_m(t) R_m(n(t), n_m(t), z, z_m).
\end{aligned}
$$

The Jacobian of this system is

$$
\begin{aligned}
\mathbf{J} &= 
\begin{pmatrix}
\frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n(t+1)}{\mathrm{d}n_m(t)} \\
\frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)} 
\end{pmatrix}.
\end{aligned}
$$

Now consider some non-zero resident equilibrium, $\hat{n}>0$, without the mutant, $\hat{n}_m=0$.
Given the resident does not continually produce mutants, $\frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n(t)}\big|_{n_m=0}=0$,
the Jacobian evaluated at this equilibrium simplifies to

$$
\begin{aligned}
\mathbf{J}\big|_{n_m=0,n=\hat{n}} &= 
\begin{pmatrix}
\frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n(t+1)}{\mathrm{d}n_m(t)} \\
0 & \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)} 
\end{pmatrix}_{n_m=0,n=\hat{n}}.
\end{aligned}
$$

We can immediately see that the two eigenvalues of this upper triangular matrix are

$$
\begin{aligned}
\lambda_1 &= \frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)}\bigg|_{n_m=0,n=\hat{n}}\\
\lambda_2 &= \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)}\bigg|_{n_m=0,n=\hat{n}}.
\end{aligned}
$$

The first, $\lambda_1$, determines whether the resident equilibrium, $\hat{n}>0$, is stable in the absence of mutants. We'll take $0 < \lambda_1 < 1$ as a given (i.e., we're only interested in stable resident equilibria, we can disregard the others).

The second, $\lambda_2$, determines whether the resident equilibrium is stable in the presence of a small number of mutants. We call $\lambda_2$ the **invasion fitness**, $\lambda(z_m,z)$. The mutant will invade whenever $\lambda(z_m,z) > 1$ (in continuous-time we'd need $\lambda(z_m,z) > 0$).

In some simple cases we might be able to use the invasion criterium, $\lambda(z_m,z)>1$, to determine what values of $z_m$ (relative to $z$) can invade. In most cases, however, $\lambda(z_m,z)$ will be complex enough that this will not be possible, so we rely on a simple approximation.

!!! note "The evolution of dispersal"

    Let a fraction $d$ of the resident offspring disperse and a fraction $d_m$ of the mutant offspring disperse. And let there be $n(t)$ residents, $n_m(t)$ mutants, and $S-n(t)-n_m(t)\geq0$ empty sites. Then the probability a resident offspring replaces a resident is the number of resident offspring in a resident patch divided by the total number of offspring in that patch,

    $$
    \begin{aligned}
    p_{rr} &= \frac{B(1-d) + (n(t)-1)Bd(1-c)/S}{B(1-d) + (n(t)-1)Bd(1-c)/S + n_m(t)Bd_m(1-c)/S}\\
    &= \frac{S(1-d) + (n(t)-1)d(1-c)}{S(1-d) + (n(t)-1)d(1-c) + n_m(t)d_m(1-c)}.
    \end{aligned}
    $$

    Here $B(1-d)$ is the number of non-dispersing offspring produced by the resident in that patch, $(n(t)-1)Bd(1-c)/S$ is the number of resident offspring dispersing to the patch from elsewhere, and $n_m(t)Bd_m(1-c)/S$ is the number of mutant offspring dispersing to the patch. 
    
    The probability that a mutant offspring replaces a mutant $p_{mm}$ is the same expression with $d$ and $d_m$ and $n$ and $n_m$ exchanged,
    
    $$
    \begin{aligned}
    p_{mm} &= \frac{S(1-d_m) + (n_m(t)-1)d_m(1-c)}{S(1-d_m) + (n_m(t)-1)d_m(1-c) + n(t)d(1-c)}.
    \end{aligned}
    $$
    
    The probability that a resident offspring wins an empty patch is, analogously,

    $$
    \begin{aligned}
    p_{re} &= \frac{n(t)Bd(1-c)}{n(t)Bd(1-c) + n_m(t)Bd_m(1-c)}\\
    &= \frac{n(t)d}{n(t)d + n_m(t)d_m}.
    \end{aligned}
    $$
    
    These probabilities allow us to calculate the expected number of resident and mutant individuals in the next generation,

    $$
    \begin{aligned}
    n(t+1) &= n(t) p_{rr} + n_m(t) (1 - p_{mm}) + (S - n(t) - n_m(t)) p_{re}\\
    n_m(t+1) &= n(t) (1-p_{rr}) + n_m(t) p_{mm} + (S - n(t) - n_m(t)) (1-p_{re}).
    \end{aligned}
    $$    

    Now consider an equilibrium with no mutants, $\hat{n}_m=0$. Setting $n(t+1) = n(t) = \hat{n}$ gives $\hat{n}=S$.

    We determine the stability of this equilibrium with the Jacobian. From the general analysis above we know how to calculate the two eigenvalues,

    $$
    \begin{aligned}
    \lambda_1 &= \frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)}\bigg|_{n_m=0,n=S} = 0\\
    \lambda_2 &= \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)}\bigg|_{n_m=0,n=S} = \frac{S(1-d_m) - d_m(1-c)}{S(1-d_m) + Sd(1-c) - d_m(1-c)} + \frac{Sd_m(1-c)}{S(1-d) + d(1-c)(S-1)}.
    \end{aligned}
    $$

    The first eigenvalue is less than 1 in absolute value, meaning the resident equilibrium is stable. 
    
    The second eigenvalue is our invasion fitness, $\lambda(d_m,d)=\lambda_2$. The mutant invades when this is greater than one. This simplifies a little in the limit of a large number of sites,

    $$
    \begin{aligned}
    \lim_{S\rightarrow\infty}\lambda(d_m,d) &> 1\\
    \frac{(1-d_m)}{(1-d_m) + d(1-c)} + \frac{d_m(1-c)}{(1-d) + d(1-c)} &> 1\\
    \frac{(1-d_m)}{(1-d_m) + d(1-c)} + \frac{d_m(1-c)}{(1-d) + d(1-c)} - 1 &> 0\\
    \frac{(1-c)(d_m-d)(1-cd-d_m)}{(1-cd)((1-d_m)+d(1-c))} &> 0.
    \end{aligned}
    $$

    The sign of the left-hand side is the sign of $(d_m-d)(1-cd-d_m)$. So the mutant will invade if both terms are positive ($d<d_m<1-cd$) or both terms are negative ($1-cd<d_m<d$). In either case, the mutant will invade when it has a trait value between $1-cd$ and $d$. This is plotted below for three different resident dispersal rates. 


<pre data-executable="true" data-language="python">
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

d,n,c,n_m,d_m=symbols('d,n,c,n_m,d_m') 
invfit = (1-d_m)/((1-d_m) + d*(1-c)) + (d_m*(1-c))/((1-d) + d*(1-c)) #invasion fitness

ds = [0.5,0.75,1]
xs = np.linspace(0,1,100)

fig, axs = plt.subplots(1,len(ds),figsize=(12,3))
for i,di in enumerate(ds):
    ax = axs[i]
    inv = invfit.subs(c,0.33).subs(d,di)
    ax.plot(xs, [inv.subs(d_m,x) for x in xs])
    ax.plot(xs, [1 for _ in xs], c='k')
    ax.set_ylim(0.6,1.2)
    ax.set_title(r'resident dispersal rate $d=%s$' %di, fontsize=10)
    ax.set_xlabel(r'mutant dispersal rate $d_m$')

axs[0].set_ylabel(r'invasion fitness $\lambda(d_m,d)$')

plt.show()
</pre>


    
![png](lecture-16_files/lecture-16_6_0.png)
    


<span id='section2'></span>
## 2. Evolutionarily singular strategies
<hr>

When the mutant trait value is very close to the resident trait value, we can approximate invasion fitness with a Taylor series around $z_m = z$,

$$
\begin{aligned}
\lambda(z_m,z) &\approx \lambda(z,z) + \frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=z} (z_m-z)\\
&= 1 + \frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=z} (z_m-z).
\end{aligned}
$$

This allows us to determine which direction evolution will proceed from the current resident value:

- if $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}>0$ then invasion when $z_m>z$
- if $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}<0$ then invasion when $z_m<z$

The direction of evolution by small steps is given by $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}$, which we call the **selection gradient**.

Potential evolutionary endpoints, also called **evolutionarily singular strategies**, are the resident trait values $z=\hat{z}$ where there is no directional selection

$$
\frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=\hat{z}, z=\hat{z}} = 0.
$$

!!! note "The evolution of dispersal"

    With infinite sites, the selection gradient is

    $$
    \frac{\partial \lim_{S\rightarrow\infty}\lambda}{\partial d_m}\Big|_{d_m=d} = \frac{(1-d-dc)(1-c)}{(1-cd)^2}
    $$

    The sign of this is the sign of 1-d-dc, meaning the dispersal rate will increase whenever d<1/(1+c) and decrease whenever 1/(1+c)<d. Setting the selection gradient to zero and solving for the evolutionarily singular strategy gives 
    
    $$\hat{d} = 1/(1+c).$$


<pre data-executable="true" data-language="python">
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

d,n,c,n_m=symbols('d,n,c,n_m') 
selgrad = (1-d-d*c)*(1-c)/(1-c*d)**2 #selection gradient
cval=0.33
xs = np.linspace(0,1,100)
sel = selgrad.subs(c,cval)
plt.plot(xs, [sel.subs(d,x) for x in xs])
plt.plot(xs, [0 for _ in xs], c='k')
plt.scatter([1/(1+c).subs(c,cval)],[0],c='k',s=50,zorder=2) #singular strategy
plt.xlabel(r'resident dispersal rate $d$')
plt.ylabel(r'selection gradient')
plt.show()
</pre>


    
![png](lecture-16_files/lecture-16_10_0.png)
    


<span id='section3'></span>
## 3. Summary
<hr>

- For a given stable resident strategy we can derive the invasion fitness of a mutant, $\lambda(z_m,z)$.
- In simpler models we can determine which mutants can invade from examining $\lambda(z_m,z)>1$ (discrete time) or $\lambda(z_m,z)>0$ (continuous time)
- More generally, we can determine the direction of evolution by small mutations from the selection gradient, $\left.\frac{\mathrm{d}\lambda}{\mathrm{d}z_m}\right|_{z_m=z}$
- Evolutionary singular strategies $\hat z$ are trait values where there is no selection, $\left.\frac{\mathrm{d}\lambda}{\mathrm{d}z_m}\right|_{z_m=\hat z}=0$
- The stability of singular strategies will be examined in the next lecture

Practice questions from the textbook: 12.1-2.
