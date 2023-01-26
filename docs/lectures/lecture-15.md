<script type="text/x-thebe-config">
  {
      requestKernel: true,
      mountActivateWidget: true,
      mountStatusWidget: true,
      binderOptions: {
      repo: "mmosmond/executable-cells",
      ref: "main",
      binderUrl: "https://gke.mybinder.org",
      },
  }
</script>
<script src="https://unpkg.com/thebe@latest/lib/index.js"></script>

# Lecture 15: Equilibria and stability (nonlinear multivariate)

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

Now that we've covered linear multivariate models, we turn our attention to the most common type of model: one with multiple *interacting* variables.

For example, a model of the number of susceptible, $S$, and infected, $I$, individuals often includes the interaction between these two variables, $SI$, describing the rate at which these two classes of individuals meet one another (and potentially spread disease). 

Similarly, models of predator, $P$, and prey, $N$, abundance often include terms like $NP$ describing the rate at which predators encounter prey (and potentially eat them).

Let's first see how to deal with these models in general and then apply those techniques to specific circumstances like those mentioned above (in the next two lectures).

<span id='section1'></span>
## 1. Continuous time
<hr>

In general, if we have $n$ interacting variables, $x_1, x_2, ..., x_n$, we can write any continuous time model like

$$
\begin{aligned}
\frac{\mathrm{d}x_1}{\mathrm{d}t} &=  f_1(x_1, x_2, ..., x_n)\\
\frac{\mathrm{d}x_2}{\mathrm{d}t} &=  f_2(x_1, x_2, ..., x_n)\\
&\vdots\\
\frac{\mathrm{d}x_n}{\mathrm{d}t} &=  f_n(x_1, x_2, ..., x_n)  
\end{aligned}
$$

If we then want to find the equilibria, $\hat{x}_1,\hat{x}_2, ..., \hat{x}_n$, we set all these equations to 0 and solve for one variable at a time (note that solving for the equilibrium is not always possible in nonlinear models!).

Now note that we can no longer write this system of equations in matrix form with a matrix composed only of parameters.
 
In order to use what we've learned about eigenvalues and eigenvectors we're first going to have to **linearize** the system so that the corresponding matrices do not contain variables.

As we saw in nonlinear univariate models, one useful way to linearize a system is to measure the system relative to equilibrium, $\epsilon = n - \hat{n}$.

Then assuming that the deviation from equilibrium, $\epsilon$, is small, we used a Taylor series expansion to approximate the nonlinear system with a linear system.

To do that with multivariate models we'll need to know how to take a Taylor series expansion of multivariate functions

!!! note "Taylor series expansion of a multivariate function"

    Taking the series of $f$ around $x_1=a_1$, $x_2=a_2$, ..., $x_n=a_n$ gives

    $$
    \begin{aligned}
    f(x_1, x_2, ..., x_n) &= f(a_1, a_2, ..., a_n)\\ 
    &+ \sum_{i=1}^{n} \left( \frac{\partial f}{\partial x_i} \bigg|_{x_1=a_1, x_2=a_2, ..., x_n=a_n} \right) (x_i - a_i)\\
    &+ \sum_{i=1}^{n}\sum_{j=1}^n \left( \frac{\partial f}{\partial x_i \partial x_j} \bigg|_{x_1=a_1, x_2=a_2, ..., x_n=a_n} \right) (x_i - a_i)(x_j - a_j)\\
    &+ \cdots
    \end{aligned}
    $$
    
    where $\frac{\partial f}{\partial x_i}$ is the "partial derivative" of $f$ with respect to $x_i$, meaning that we treat all the other variables as constants when taking the derivative. 

Then when the difference between each variable and its value, $x_i-a_i$, is small enough we can ignore all the terms with a $(x_i-a_i)(x_j-a_j)$, and we are left with a **linear** approximation of $f$.

So let $\epsilon_i = x_i - \hat{x}_i$ be the deviation of variable $x_i$ from its equilibrium value, $\hat{x}_i$.
  
Then we can write a system of equations describing the change in the deviations for all of our variables

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &=  f_1(x_1, x_2, ..., x_n)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &=  f_2(x_1, x_2, ..., x_n)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &=  f_n(x_1, x_2, ..., x_n)  
\end{aligned}
$$

And then we can take a Taylor series around $x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n$ to get a linear approximation of our system near the equilibrium

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx  f_1(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx  f_2(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  f_n(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)
\end{aligned}
$$

And then note that all $f_i(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n)=0$ by definition of a equilibrium, leaving

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i - \hat{x}_i)
\end{aligned}
$$

Each of the partial derivatives $\frac{\partial f_i}{\partial x_j}$ is evaluated at the equilibrium, so these are constants. And $x_i - \hat{x}_i = \epsilon_i$. So we now have a linear system

$$
\begin{aligned}
\frac{\mathrm{d}\epsilon_1}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i\\
\frac{\mathrm{d}\epsilon_2}{\mathrm{d}t} &\approx \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i\\
&\vdots\\
\frac{\mathrm{d}\epsilon_n}{\mathrm{d}t} &\approx  \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i
\end{aligned}
$$

We can now write our approximate system around the equilibrium in matrix form

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
\begin{pmatrix} \epsilon_1 \\ \epsilon_2 \\ \vdots \\ \epsilon_n \end{pmatrix}
\end{aligned}
$$

This is a special matrix called the **Jacobian**

$$
\mathbf{J} = 
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{pmatrix}
$$

And now that we have a linear system around an equilibrium, we can assess its **local stability** just as we did with linear multivariate models (see [Summary](#section3)).

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
x_n(t+1) &=  f_n(x_1(t), x_2(t), ..., x_n(t))\\
\end{aligned}
$$

Now the equilibria are found by setting all $x_i(t+1) = x_i(t) = \hat{x}_i$ and solving for the $\hat{x}_i$ one at a time.

To linearize the system around an equilibrium we again measure the system in terms of deviation from the equilibrium, $\epsilon_i(t) = x_i(t) - \hat{x}_i$, giving

$$
\begin{aligned}
\epsilon_1(t+1) &=  f_1(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1\\
\epsilon_2(t+1) &=  f_2(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1\\
&\vdots\\
\epsilon_n(t+1) &=  f_n(x_1(t), x_2(t), ..., x_n(t)) - \hat{x}_1\\
\end{aligned}
$$

Then taking the Taylor series of each $f_i$ around $x_1(t) = \hat{x}_1, ..., x_n(t) = \hat{x}_n$ we can approximate our system near the equilibrium as

$$
\begin{aligned}
  \epsilon_1(t+1) &=  f_1(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n)  + \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_1\\
  \epsilon_2(t+1) &=  f_2(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_2\\
  &\vdots\\
  \epsilon_n(t+1) &=  f_n(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) + \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) (x_i(t) - \hat{x}_i) - \hat{x}_n\\
  \end{aligned}
$$

And noting that $f_i(\hat{x}_1, \hat{x}_2, ..., \hat{x}_n) = \hat{x}_i$ and $x_i(t) - \hat{x}_i = \epsilon_i(t)$ we have

$$
\begin{aligned}
  \epsilon_1(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_1}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t)\\
  \epsilon_2(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_2}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t)\\
  &\vdots\\
  \epsilon_n(t+1) &= \sum_{i=1}^{n} \left( \frac{\partial f_n}{\partial x_i} \bigg|_{x_1=\hat{x}_1, x_2=\hat{x}_2, ..., x_n=\hat{x}_n} \right) \epsilon_i(t)\\
  \end{aligned}
$$

We can therefore write our approximation around the equilibrium in matrix form

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
\begin{pmatrix} \epsilon_1(t) \\ \epsilon_2(t) \\ \vdots \\ \epsilon_n(t) \end{pmatrix}
\end{aligned}
$$

As in continuous time, the dynamics are described by the **Jacobian** matrix

$$
\mathbf{J} = 
\begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{pmatrix}
$$

We therefore assess the **local stability** of an equilibrium by evaluating the Jacobian at that equilibrium and finding the eigenvalues (see [Summary](#section3)).

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

    - the leading eigenvalue is the one with the **largest real part**
    - if the leading eigenvalue has a **negative real part** the equilibrium is stable
    - if any eigenvalue has a non-zero complex part there will be cycling

    **Discrete time**

    - the leading eigenvalue is the one with the **largest absolute value**
    - if the leading eigenvalue has an **absolute value less than one** the equilibrium is stable
    - if any eigenvalue has a non-zero complex part there will be cycling
