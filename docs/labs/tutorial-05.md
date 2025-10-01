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

# Tutorial 5: long-term dynamics of linear multivariate models

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

## Discrete time

We now know that the general solution of a system of linear recursion equations, $\vec{x}(t+1) = \mathbf{M}\vec{x}(t) + \vec{m}$, can be written 

$$
\vec{x}(t) = \mathbf{A}\mathbf{D}^t\mathbf{A}^{-1}(\vec{x}(0)-\hat{\vec{x}}) + \hat{\vec{x}},
$$ 

where $\hat{\vec{x}}=-(\mathbf{M} - \mathbf{I})^{-1}\vec{m}$ is the equilibrium, $\mathbf{A}$ is a matrix with right eigenvectors as columns, and $\mathbf{D}$ is a diagonal matrix with eigenvalues along the diagonal.

Note that as time proceeds, the eigenvalues in $\mathbf{D}$ get raised to higher and higher powers. As a result, the eigenvalue with the largest absolute value, which we call the **leading eigenvalue**, comes to dominate. To see this, let the leading eigenvalue be $\lambda_1$ and factor it out of $\mathbf{D}^t$,

$$
\mathbf{D}^t = \lambda_1^t 
\begin{pmatrix} 
1 & 0 & \cdots & 0\\
0 & (\lambda_2/\lambda_1)^t & \cdots & 0\\
\vdots & \vdots & \vdots & \vdots\\ 
0 & 0 & \cdots & (\lambda_n/\lambda_1)^t\\  
\end{pmatrix}.
$$

Since $|\lambda_i/\lambda_1|<1$ for all $i$, for large $t$ these all go to zero and we have

$$
\mathbf{D}^t \approx \tilde{\mathbf{D}}^t  \equiv \lambda_1^t 
\begin{pmatrix} 
1 & 0 & \cdots & 0\\
0 & 0 & \cdots & 0\\
\vdots & \vdots & \vdots & \vdots\\ 
0 & 0 & \cdots & 0\\  
\end{pmatrix}.
$$

We can therefore approximate $\vec{x}(t)$ after a sufficient amount of time as,

$$
\begin{aligned}
\tilde{\vec{x}}(t) &= \mathbf{A}\tilde{\mathbf{D}}^t\mathbf{A}^{-1}(\vec{x}(0)-\hat{\vec{x}}) + \hat{\vec{x}}\\
&= \lambda_1^t \vec{v}_1 \vec{u}_1 (\vec{x}(0)-\hat{\vec{x}}) + \hat{\vec{x}}
\end{aligned}
$$

where $\vec{v}_1$ is the right eigenvector associated with the leading eigenvalue (the first column of $\mathbf{A}$) and $\vec{u}_1$ is the left eigenvector associated with the leading eigenvalue, $\lambda_1$ (the first row of $\mathbf{A}^{-1}$).

!!! warning

    Finding the left eigenvectors via the inverse of $\mathbf{A}$ guarantees that the eigenvectors have been scaled such that $\vec{u}_i\vec{v}_i = 1$. If you instead find the left eigenvectors by solving $\vec{u}\mathbf{M}=\lambda\vec{u}$, make sure you scale $\vec{u}_1$ so that $\vec{u}_1\vec{v}_1 = 1$, eg, divide your unscaled left eigenvector, $\vec{u}_1'$, by $\vec{u}_1'\vec{v}_1$. Otherwise the long-term approximation will be off by a constant factor.

This convenient approximation tells us two things:

- the system will eventually converge to the equilibrium $\hat{\vec{x}}$ if and only if the leading eigenvalue has absolute value less than one, $|\lambda_1|<1$
- since $\vec{u}_1 (\vec{x}(0)-\hat{\vec{x}})$ is a scalar (row vector times column vector), $\vec{v}_1$ tells us the relative deviation of each $x_i$ from its equilibrium value in the long term, eg, if $\vec{v}_1= \begin{pmatrix} 2 \\ 1 \end{pmatrix}$ then $x_1$ is eventually twice as far from its equilibrium value than $x_1$ is

## Continuous time

In continuous time we the general solution is

$$
\vec{x}(t) = \mathbf{A}\exp(\mathbf{D}t)\mathbf{A}^{-1}(\vec{x}(0) - \hat{\vec{x}}) + \hat{\vec{x}},
$$

with $\hat{\vec{x}}=-\mathbf{M}^{-1}\vec{m}$.

As $t$ increases $\exp(\mathbf{D}t)$ becomes dominated by the eigenvalue with the largest (ie, most positive) value, which we call the leading eigenvalue. The long term dynamics can then be approximated 

$$
\vec{x}(t) \approx \exp(\lambda_1 t) \vec{v}_1 \vec{u}_1 (\vec{x}(0)-\hat{\vec{x}}) + \hat{\vec{x}}.
$$

Therefore

- the system will eventually converge to the equilibrium $\hat{\vec{x}}$ if and only if the leading eigenvalue is negative, $\lambda_1<0$
- $\vec{v}_1$ tells us the relative deviation of each $x_i$ from its equilibrium value in the long term

## Problem

Many  eukaryotic  genes  have  both  exons  and  introns,  where  only  exons code  for  protein  sequence.  Messenger RNAs (mRNAs)  transcribed  from  such  genes  initially  include the  introns,  which  must  be  spliced  out  before  the  mRNAs  can  be  translated into proteins.  Let $n_1$ be the number of unspliced “pre-mRNAs,” which are produced by the cell at rate $c$ and let $n_2$ be the number of spliced “processed mRNAs”.  If pre-mRNAs are  spliced at rate $a$ and processed mRNAs degrade at rate $d$, the numbers of pre-mRNAs and processed mRNAs change according to

$$
\begin{align}
\frac{\mathrm{d}n_1}{\mathrm{d}t} &= c - a n_1 \\ 
\frac{\mathrm{d}n_2}{\mathrm{d}t} &= a n_1 - d n_2 .
\end{align}
$$

Writing this in matrix form, $\mathrm{d}\vec{n}/\mathrm{d}t = \mathbf{M}\vec{n} + \vec{m}$, determine $\mathbf{M}$ and $\vec{m}$. Determine the eigenvalues (and associated eigenvectors if you have time). Is the equilibrium stable? 


<pre data-executable="true" data-language="python">

</pre>
