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

# Lecture 12: Stability (nonlinear multivariate)

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

1. [Continuous time](#section1)
2. [Discrete time](#section2)
5. [Summary](#section3)

Now that we know how to find equilibria in nonlinear multivariate models, let's see how to determine if they're stable.

<span id='section1'></span>
## 1. Continuous time
<hr>

In general, if we have $n$ interacting variables, $x_1, x_2, ..., x_n$, we can write any continuous time model like

$$
\begin{aligned}
\frac{\mathrm{d}x_1}{\mathrm{d}t} &=  f_1(x_1, x_2, ..., x_n)\\
\frac{\mathrm{d}x_2}{\mathrm{d}t} &=  f_2(x_1, x_2, ..., x_n)\\
&\vdots\\
\frac{\mathrm{d}x_n}{\mathrm{d}t} &=  f_n(x_1, x_2, ..., x_n).
\end{aligned}
$$

If we then want to find the equilibria, $\hat{x}_1,\hat{x}_2, ..., \hat{x}_n$, we set all these equations to 0 and solve for one variable at a time. Let's assume we've already done this and now want to know when an equilibrium is stable.

In the linear multivariate case we determined stability with the eigenvalues of a matrix $\textbf{M}$ of parameters. Unfortunately we can not generally write a system of nonlinear equations in matrix form with a matrix composed only of parameters. In order to use what we've learned about eigenvalues and eigenvectors we're first going to have to **linearize** the system so that the corresponding matrices do not contain variables.

As we saw in nonlinear univariate models, one useful way to linearize a system is to measure the system relative to equilibrium, $\epsilon = x - \hat{x}$. Then assuming that the deviation from equilibrium, $\epsilon$, is small, we used a Taylor series expansion to approximate the nonlinear system with a linear system. To do that with multivariate models we'll need to take a Taylor series expansion of a multivariate function.

??? note "Taylor series expansion of a multivariate function"

    Taking the series of $f$ around $x_1=a_1$, $x_2=a_2$, ..., $x_n=a_n$ gives

    $$
    \begin{aligned}
    f(x_1, x_2, ..., x_n) &= f(a_1, a_2, ..., a_n)\\ 
    &+ \sum_{i=1}^{n} \left( \frac{\partial f}{\partial x_i} \bigg|_{x_1=a_1, x_2=a_2, ..., x_n=a_n} \right) (x_i - a_i)\\
    &+ \sum_{i=1}^{n}\sum_{j=1}^n \left( \frac{\partial f}{\partial x_i \partial x_j} \bigg|_{x_1=a_1, x_2=a_2, ..., x_n=a_n} \right) (x_i - a_i)(x_j - a_j)\\
    &+ \cdots
    \end{aligned}
    $$
    
    where $\frac{\partial f}{\partial x_i}$ is the "partial derivative" of $f$ with respect to $x_i$, meaning that we treat all the other variables as constants when taking the derivative. Considering only the first two terms gives a linear approximation.

So let $\epsilon_i = x_i - \hat{x}_i$ be the deviation of variable $x_i$ from its equilibrium value, $\hat{x}_i$, and write a system of equations describing the change in the deviations for all of our variables,

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &=  f_1(x_1, x_2, ..., x_n)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &=  f_2(x_1, x_2, ..., x_n)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &=  f_n(x_1, x_2, ..., x_n).  
\end{aligned}
$$

We can take a Taylor series around $x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n$ to get a linear approximation of our system near the equilibrium,

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx  f_1(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx  f_2(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  f_n(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i).
\end{aligned}
$$

And then note that all $f_i(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n)=0$ by definition of a equilibrium, leaving

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i).
\end{aligned}
$$

Each of the partial derivatives $\frac{\partial f_i}{\partial x_j}$ is evaluated at the equilibrium, so these are constants. And $x_i - \hat{x}_i = \epsilon_i$. So we now have a linear system,

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i,
\end{aligned}
$$

which we can write in matrix form

$$
\begin{aligned}
\begin{pmatrix} \frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} \\ \frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} \\ \vdots \\ \frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} \end{pmatrix}
=
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{pmatrix}_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n}
\begin{pmatrix} \epsilon_1 \\ \epsilon_2 \\ \vdots \\ \epsilon_n \end{pmatrix}.
\end{aligned}
$$

The matrix of first derivatives is called the **Jacobian**,

$$
\mathbf{J} = 
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{pmatrix}.
$$

If any of the equations, $f_i$, are nonlinear then $\textbf{J}$ depends on the variables. To determine **local stability** we find the eigenvalues of the Jacobian evaluated at an equilibrium.

### Example: epidemiology

Let's return to our model of disease spread from the previous lecture. The rate of change in the number of susceptibles and infecteds is

$$\begin{aligned}
\frac{\mathrm{d}S}{\mathrm{d}t} &= \theta - \beta S I - d S + \gamma I \\
\frac{\mathrm{d}I}{\mathrm{d}t} &= \beta S I - (d + v) I - \gamma I. 
\end{aligned}$$

We found two equilibria, the diease-free equilibrium,

$$\begin{aligned}
\hat{S} &= \theta/d \\
\hat{I} &= 0,
\end{aligned}$$

which is always valid, and the endemic equilirbium,

$$\begin{aligned}
\hat{S} &= (d + v + \gamma)/\beta \\
\hat{I} &= \frac{\theta - d(d + v + \gamma)/\beta}{d+v},
\end{aligned}$$

which is valid when $\beta\theta/d - (d + v + \gamma)>0$.

Now when are these stable? 

We first calculate the Jacobian,

$$\begin{aligned}
\mathbf{J} 
&= 
\begin{pmatrix}
\frac{\partial}{\partial S}\left(\frac{\mathrm{d}S}{\mathrm{d}t}\right) & \frac{\partial}{\partial I}\left(\frac{\mathrm{d}S}{\mathrm{d}t}\right) \\
\frac{\partial}{\partial S}\left(\frac{\mathrm{d}I}{\mathrm{d}t}\right) & \frac{\partial}{\partial I}\left(\frac{\mathrm{d}I}{\mathrm{d}t}\right)
\end{pmatrix}\\
&=
\begin{pmatrix}
-d-\beta I & -\beta S+\gamma \\
\beta I & \beta S-(d+v+\gamma)
\end{pmatrix}.
\end{aligned}$$

We can now determine the local stability of an equilibrium by evaluating the Jacobian at that equilibrium and calculating the eigenvalues. 

Let's do that first for the simpler disease-free equilibrium. Plugging $\hat{S}$ and $\hat{I}$ into the Jacobian gives 

$$\begin{aligned}
\mathbf{J}_\mathrm{disease-free} 
&= 
\begin{pmatrix}
-d & -\beta \theta/d+\gamma \\
0 & \beta \theta/d-(d+v+\gamma)
\end{pmatrix}.
\end{aligned}$$

This is an upper triangular matrix, so the eigenvalues are just the diagonal elements, $\lambda = -d, \beta\theta/d-(d+v+\gamma)$. Because all the parameters are rates they are all non-negative, and therefore the only eigenvalue that can have a positive real part (and therefore cause instability) is $\lambda=\beta\theta/d-(d+v+\gamma)$. The equilibrium is unstable when this is positive, $\beta\theta/d-(d+v+\gamma)>0$. Because this equilibrium has no infected individuals, instability in this case means the infected individuals will increase in number from rare -- ie, the disease can invade. 

We can rearrange the instability condition to get a little more intuition. The disease can invade whenever

$$\begin{aligned}
\beta\theta/d - (d+v+\gamma)& > 0 \\
\beta\theta/d &> d+v+\gamma \\
\frac{\beta\theta/d}{d+v+\gamma} &> 1.
\end{aligned}$$

The numerator is $\beta$ times the number of susceptibles at the disease-free equilibrium, $\hat{S}=\theta/d$. This is the rate that a rare disease infects new individuals. The denominator is the rate at which the disease is removed from the population. Therefore a rare disease that infects faster than it is removed can spread. This ratio, in our case $\frac{\beta\theta/d}{d+v+\gamma}$, turns out to be the expected number of new infections per infection when the disease is rare. This is termed $R_0$, which is a very key epidemiological quantity (you may remember estimates of $R_0$ in the news from a certain recent virus...).

Now for the endemic equilibrium. Plugging these values into the Jacobian and simplifying gives

$$\begin{aligned}
\mathbf{J}_\mathrm{endemic} 
&= 
\begin{pmatrix}
-\frac{\beta \theta - d \gamma}{d+v} & -(d+v) \\
\frac{\beta \theta - d (d+v+\gamma)}{d+v} & 0
\end{pmatrix}
\end{aligned}$$

Here, instead of calculating the eigenvalues explicitly, we can use the **Routh-Hurwitz stability criteria** for a 2x2 matrix, a negative trace and positive determinant. The determinant is $\beta \theta - d (d+v+\gamma)$. For this to be positive we need $\beta \theta/d > (d+v+\gamma)$, which was the instability condition on the disease-free equilibrium ($R_0>1$). The trace is $-\frac{\beta \theta - d \gamma}{d+v}$. For this to be negative we need $\beta \theta/d > \gamma$, which is guaranteed if the determinant is positive. So in conclusion, the endemic equilibrium is valid and stable whenever the disease can invade, $R_0>1$.

<span id='section2'></span>
## 2. Discrete time
<hr>

!!! note

    The short version of this section is that we can do the same thing in discrete time -- local stability is determined by the eigenvalues of the Jacobian, where the functions in that Jacobian are now our recursions, $x_i(t+1) = f_i(x_1(t), x_2(t), ..., x_n(t))$.

We can do something very similar for nonlinear multivariate models in discrete time

$$
\begin{aligned}
x_1(t+1) &=  f_1(x_1(t), x_2(t), ..., x_n(t))\\
x_2(t+1) &=  f_2(x_1(t), x_2(t), ..., x_n(t))\\
&\vdots\\
x_n(t+1) &=  f_n(x_1(t), x_2(t), ..., x_n(t)).
\end{aligned}
$$

Now the equilibria are found by setting all $x_i(t+1) = x_i(t) = \hat{x}_i$ and solving for the $\hat{x}_i$ one at a time. Let's assume we've done this.

To linearize the system around an equilibrium we again measure the system in terms of deviation from the equilibrium, $\epsilon_i(t) = x_i(t) - \hat{x}_i$, giving

$$
\begin{aligned}
\epsilon_1(t+1) &=  f_1(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1\\
\epsilon_2(t+1) &=  f_2(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1\\
&\vdots\\
\epsilon_n(t+1) &=  f_n(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1.
\end{aligned}
$$

Then taking the Taylor series of each $f_i$ around $x_1(t) = \hat{x}_1, ..., x_n(t) = \hat{x}_n$ we can approximate our system near the equilibrium as

$$
\begin{aligned}
  \epsilon_1(t+1) &=  f_1(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n)  + \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_1\\
  \epsilon_2(t+1) &=  f_2(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_2\\
  &\vdots\\
  \epsilon_n(t+1) &=  f_n(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_n.
  \end{aligned}
$$

Noting that $f_i(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) = \hat{x}_i$ and $x_i(t) - \hat{x}_i = \epsilon_i(t)$ we have a linear system,

$$
\begin{aligned}
  \epsilon_1(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t)\\
  \epsilon_2(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t)\\
  &\vdots\\
  \epsilon_n(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t),
  \end{aligned}
$$

which can be written in matrix form

$$
\begin{aligned}
\begin{pmatrix} \epsilon_1(t+1) \\ \epsilon_2(t+1) \\ \vdots \\ \epsilon_n(t+1) \end{pmatrix}
=
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{pmatrix}_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n}
\begin{pmatrix} \epsilon_1(t) \\ \epsilon_2(t) \\ \vdots \\ \epsilon_n(t) \end{pmatrix}.
\end{aligned}
$$

As in continuous time, the dynamics near an equilibrium are described by the **Jacobian** matrix,

$$
\mathbf{J} = 
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}.
\end{pmatrix}
$$

We assess the **local stability** of an equilibrium by evaluating the Jacobian at that equilibrium and finding the eigenvalues.

!!! note "to-do"

    Give density-dependent natural selection example. For now see pages 320-322 in the textbook.


<span id='section3'></span>
## 3. Summary
<hr>

We can determine the stability of nonlinear multivariate models with the eigenvalues of the Jacobian evaluated at an equilibrium. The recipe is
 
- Find the equilibrium of interest, $\hat{x}_1, \hat{x}_2, ... \hat{x}_n$
- Calculate the Jacobian, $\mathbf{J}$
- Evaluate the Jacobian at the equilibrium of interest, $\mathbf{J}_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n}$
- Calculate the characteristic polynomial $|\mathbf{J}_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} - \lambda\mathbf{I}|$ 
- Set the characteristic polynomial equal to 0 and solve for the $n$ eigenvalues, $\lambda$
  
!!! info "Stability reminder"

    **Continuous time**

    - if all eigenvalues have a **negative real part** the equilibrium is stable
    - if any eigenvalue has a non-zero imaginary part there will be cycling

    **Discrete time**

    - if all eigenvalues have an **absolute value less than one** the equilibrium is stable
    - if any eigenvalue has a non-zero imaginary part there will be cycling

Practice questions from the textbook: 8.1, 8.4, 8.5, 8.6, 8.7, 8.8, 8.12.


<pre data-executable="true" data-language="python">

</pre>
