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

# Lecture 18: Evolutionary invasion analysis

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
2. [Evolutionary singular strategies](#section2)
3. [Evolutionary stable strategies](#section3)
4. [Evolutionary convergence](#section4)
5. [Pairwise invasibility plots](#section5)
6. [Summary](#section6)

In the models we've discussed so far we've taken the parameters to be fixed. In reality, many of these parameters can evolve. For example, in our model of exponential growth in discrete time \[n(t+1)=n(t) R\] we took $R$ to be the same for all individuals for all time. But clearly any mutation causing a larger $R$ would increase in frequency, causing the value of $R$ to increase over time. In this lecture we'll explore how to determine the direction of evolution and the stability of evolutionary endpoints for more complex models using a technique called <b>evolutionary invasion analysis</b>.

The idea:
  - determine which parameters of our model an evolving trait affects
  - take the population to be fixed for some "resident" trait value
  - determine the equilibria and stability of the system with only the resident trait
  - derive an equation for the growth of a rare "mutant" allele that affects the trait
  - ask when the mutant will invade
  - look for potential evolutionary endpoints
  - determine the stability of those endpoints

<span id='section1'></span>
## 1. Invasion fitness
<hr>

Let's think about this analysis very generally (in discrete time).

Let the number of individuals with the resident trait value be $n$ and the number of individuals with the mutant trait value $n_m$. (And we'll assume asexual reproduction for simplicity, so that residents produce residents and mutant produce mutants.)

Let the potentially nonlinear dynamics of these two groups of individuals depend on their respective trait values, $z$ and $z_m$,

\[\begin{aligned}
n(t+1) &= n(t) R(n(t), n_m(t), z, z_m)\\
n_m(t+1) &= n_m(t) R_m(n(t), n_m(t), z, z_m)
\end{aligned}\]

The Jacobian of this system is

\[\begin{aligned}
\mathbf{J} &= 
\begin{pmatrix}
\frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n(t+1)}{\mathrm{d}n_m(t)} \\
\frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)} 
\end{pmatrix}
\end{aligned}\]

Now consider some non-zero resident equilibrium, $\hat{n}>0$, without the mutant, $\hat{n}_m=0$.

Assuming that the resident does not produce mutants, $\frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n(t)}\big|_{n_m=0}=0$,
the Jacobian evaluated at this equilibrium simplifies to

\[\begin{aligned}
\mathbf{J}\big|_{n_m=0,n=\hat{n}} &= 
\begin{pmatrix}
\frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)} & \frac{\mathrm{d}n(t+1)}{\mathrm{d}n_m(t)} \\
0 & \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)} 
\end{pmatrix}_{n_m=0,n=\hat{n}}
\end{aligned}\]

We can immediately see that the two eigenvalues of this upper triangular matrix are

\[\begin{aligned}
\lambda_1 &= \frac{\mathrm{d}n(t+1)}{\mathrm{d}n(t)}\bigg|_{n_m=0,n=\hat{n}}\\
\lambda_2 &= \frac{\mathrm{d}n_m(t+1)}{\mathrm{d}n_m(t)}\bigg|_{n_m=0,n=\hat{n}}\\
\end{aligned}\]

The first, $\lambda_1$, determines whether the resident equilibrium, $\hat{n}>0$, is stable in the absence of mutants. We'll take $0< \lambda_1 < 1$ as given.

The second, $\lambda_2$, determines whether the resident equilibrium is stable in the presence of a small number of mutants. We call $\lambda_2$ the \textbf{invasion fitness}, $\lambda(z_m,z)$.

In particular, the mutant will invade whenever $\lambda(z_m,z) > 1$.

In some simple cases we might be able to use the invasion criterium, $\lambda(z_m,z)>1$, to determine what values of $z_m$ (relative to $z$) can invade. In most cases, however, $\lambda(z_m,z)$ will be complex enough that this will not be possible, so we rely on a simple approximation.

!!! todo

    Give example where we can look at \(\lambda(z_m,z)>1\) directly.

<span id='section2'></span>
## 2. Evolutionary singular strategies
<hr>

When the mutant trait value is very close to the resident trait value, we can use a first order Taylor series approximation around $z_m = z$

\[\begin{aligned}
\lambda(z_m,z) &\approx \lambda(z,z) + \frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=z} (z_m-z)\\\pause
&= 1 + \frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=z} (z_m-z)
\end{aligned}\]

This allows us to determine which direction evolution will proceed:

\begin{itemize}
\item[] $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}>0 \implies$ invasion, $\lambda(z_m,z)>1$, when $z_m>z$\pause
\item[] $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}<0 \implies$ invasion, $\lambda(z_m,z)>1$, when $z_m<z$
\end{itemize}

The direction of evolution by small steps is given by $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}$, the selection gradient.

Potential evolutionary endpoints, also called <b>evolutionarily singular strategies</b>, are the resident trait values $z=\hat{z}$ where there is no directional selection

\[\frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=\hat{z}, z=\hat{z}} = 0\]


<span id='section3'></span>
## 3. Evolutionary stable strategies
<hr>

An evolutionarily singular strategy, $\hat{z}$, will only be an <b>evolutionarily stable strategy</b> (ESS), $z^*$, if it cannot be invaded.

In the evolution of dispersal example, we could show that $\hat{d} = \frac{1}{1+c}$ was \textit{globally} stable because $\lambda(d_m,d)|_{d=\hat{d}}<1$ for all $d_m$.

Global stability will be impossible to prove for more complex models, and so we often focus on \textit{local} stability, which requires that $\lambda(z_m,z)|_{z=\hat{z}}$ is concave at $z_m=\hat{z}$

\[\frac{\partial^2 \lambda}{\partial z_m^2}\bigg|_{z_m=\hat{z}, z=\hat{z}} < 0\]

To summarize, an evolutionarily stable strategy, $z^*$, satisfies both

\[\begin{aligned}
\frac{\partial \lambda}{\partial z_m}\bigg|_{z_m=z^*, z=z^*} &= 0\\
\frac{\partial^2 \lambda}{\partial z_m^2}\bigg|_{z_m=z^*, z=z^*} &< 0
\end{aligned}\]

i.e., $z^*$ is a (local) fitness maximum.

<span id='section4'></span>
## 4. Evolutionary convergence
<hr>

There is one more characteristic of evolutionarily singular strategies that we care about, and
that is whether evolution actually leads to that strategy or not. For evolution to move the trait value towards a singular strategy, $\hat{z}$, we need evolution to increase the trait value when it is less than $\hat{z}$ and decrease the trait value when it is greater than $\hat{z}$. In other words, we need the selection gradient $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z}$ to decrease as we move through $z=\hat{z}$

\[\frac{\mathrm{d}}{\mathrm{d} z}\left( \frac{\partial \lambda}{\partial z_m}\Big|_{z_m=z} \right)_{z=\hat{z}} < 0\]

Singular strategies that satisfy this criteria are called <b>convergence stable</b>. Interestingly, not all evolutionarily stable strategies are convergence stable and not all convergence stable singular strategies are evolutionarily stable! Evolutionarily stable strategies that are not convergence stable are called <b>Garden of Eden strategies</b>.Singular strategies that are convergence stable but not evolutionarily stable are called <b>evolutionary branching points</b>. The latter are of particular interest because the system evolves towards a state where multiple strategies can coexist, leading to diversification.

<span id='section5'></span>
## 5. Pairwise invasibility plots
<hr>

A helpful way to visualize the two types of stability at an evolutionarily singular strategy is called a <b>pairwise invasibility plot (PIP)</b>

!!! to do

    show PIP

<span id='section6'></span>
## 6. Summary
<hr>

There are four types of evolutionarily singular strategies $\hat{z}$, where $\frac{\partial \lambda}{\partial z_m}\Big|_{z_m=\hat{z},z=\hat{z}}=0$.

!!! to do

    show PIPs for each type
