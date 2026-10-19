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

# Tutorial 6: Stability in multivariate linear models

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

Records from the Hudson's Bay Company famously show ~10-year population size cycles in hare and lynx. Moran (1953) built a model to fit to the lynx data collected from the Mackenzie River district from 1821 to 1934. 

Let $Y(t)$ be the number of lynx caught in year $t$. Following Moran, we will work with the logarithm, $X(t) = \log_{10} Y(t)$. The expected deviation from the mean, $X(t) - m$ (where $m$ is the mean value of $X$ over years) is predicted to change over time like 

$$
x(t+1) = a x(t) + b x(t-1).
$$

Here, $a$ and $b$ are constants that Moran estimated from data but we'll consider more generally: $a$ describes how much next year's population size is determined by the current year and $b$ describes how much next year's population size is determined by the previous year (a time-lag). These parameters could be negative (eg, big population sizes in the past mean less resources for the future). Linear models with such time-lags are often called autoregressive. 

Consider $x(t+1)$ and $x(t)$ as our two variables, with equations

$$
\begin{align}
x(t+1) = a x(t) + b x_(t-1)\\
x(t) = x(t).
\end{align}
$$

This can be written in matrix form,

$$
\begin{pmatrix}
x(t+1)\\
x(t)
\end{pmatrix}
= \begin{pmatrix}
a & b\\
1 & 0
\end{pmatrix}
\begin{pmatrix}
x(t)\\
x(t-1)
\end{pmatrix}.
$$

If we define $\vec{x}(t)= 
\begin{pmatrix}
x(t)\\
x(t-1)
\end{pmatrix}$,
this matrix equation is simply

$$
\vec{x}(t+1) = \mathbf{M} \vec{x}(t).
$$

## Questions

Determine the equilibrium.

Calculate the eigenvalues of $\mathbf{M}$. 

Determine what needs to be true about the parameters for there to be oscillations/cycles. Interpret the condition biologically. 

Given there are cycles, determine when the equilibrium is stable.
