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

# Lecture 21: Probability II (demographic stochasticity)

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

1. [Poisson random variable](#section1)
2. [Demographic stochasticity](#section2)
3. [Extinction](#section3)
4. [Establishment](#section4)

!!! note "credits"

    This lecture was created by PhD student Puneeth Deraje as part of a course development TA position -- thanks Puneeth! If you are following along with the text, this lecture does not follow along as closely.

In [Lecture 18](lecture-18.md) we learned how to model stochasticity in population genetics (genetic drift). In this lecture we'll learn how to model stochasticity in population dynamics (**demographic stochasticity**). 

To do so we'll work with the simplest model of population dynamics possible, exponential growth (see [Lecture 3](lecture-03.md)). In discrete time this is

$$n_{t+1} = R n_t$$

where $n_t$ is the number of individuals at time $t$ and $R$ is the reproductive factor. In this deterministic model every individual had a reproductive factor of exactly $R$ at every time step. In reality there will be stochasticity in $R$ across individuals and time. We can account for that by replacing $R$ with a random variable. 

<span id='section1'></span>
## 1. Poisson random variable
<hr>

Recall the binomial random variable from the previous lecture, $X \sim \mathrm{Bin}(n,p)$. The probability this random variable takes on value $k$ is then 

$$\Pr(X = k) = {n \choose k} p^k (1-p)^{n-k}$$

We could model the reproductive factor as a binomial random variable. For example, perhaps each individual at time $t$ produces $n$ offspring before dying, and each offspring survives to become an adult with probability $p$. Using this model requires estimates of two parameters, $n$ and $p$. But there is a simpler way.

Let the mean number of offspring produced be $\lambda=np$. Rearranging this in terms of $p$, we can write the binomial as $X\sim \mathrm{Bin}(n,\lambda/n)$. Note that the mean is always $\lambda$, regardless of $n$. Now imagine that individuals in the population we are modelling tend to have a very large number of offspring, very few of which survive. We can approximate this in the extreme by taking $n\rightarrow \infty$. Let $Y$ be this random variable, $Y = \lim_{n \rightarrow \infty} \mathrm{Bin}(n, \lambda/n))$. The distribution of $Y$ is then 

$$
\begin{aligned}
\Pr(Y=k) &= \lim_{n \rightarrow \infty} {n \choose k} \left(\frac{\lambda}{n}\right)^k \left(1-\frac{\lambda}{n}\right)^{n-k} \\
&= \lim_{n \rightarrow \infty} \frac{n(n-1)(n-2)...(n-k+1)}{k!} \frac{\lambda^k}{n^k} \left(1-\frac{\lambda}{n}\right)^{n-k}\\
&= \frac{\lambda^k}{k!}\lim_{n \rightarrow \infty} \frac{n(n-1)(n-2)...(n-k+1)}{n^k} \left(1-\frac{\lambda}{n}\right)^{n-k}\\
&= \frac{\lambda^k}{k!}\lim_{n \rightarrow \infty} 1\left(1-\frac{1}{n}\right)\left(1-\frac{2}{n}\right)...\left(1-\frac{k-1}{n}\right) \left(1-\frac{\lambda}{n}\right)^{n-k}\\
&= \frac{\lambda^k}{k!}\lim_{n \rightarrow \infty} \left(1-\frac{\lambda}{n}\right)^{n-k}\\
&= \frac{\lambda^k e^{-\lambda}}{k!}
\end{aligned}
$$

We call $Y$ a **Poisson random variable** with mean $\lambda$ and denoted it by $Y\sim\mathrm{Poi}(\lambda)$. We can now model the reproductive factor as a random variable with only one parameter, $\lambda$. 

Note that the simpler Poisson distribution is a good approximation of the binomial distribution even for fairly moderate values of $n$, as seen in the plot below.


<pre data-executable="true" data-language="python">
import sympy
import numpy as np
import matplotlib.pyplot as plt

def binomial(n,p,k):
    return sympy.binomial(n,k) * p**k * (1-p)**(n-k)

def poisson(lam,k):
    return lam**k * np.exp(-lam)/sympy.factorial(k)

n = 100
p = 0.1
lam = n*p
ks = range(n+1)

fig, ax = plt.subplots()

ax.plot(ks, [binomial(n,p,k) for k in ks], label='binomial')
ax.plot(ks, [poisson(lam,k) for k in ks], label='Poisson')

ax.set_xlabel('number of successes')
ax.set_ylabel('probability')
ax.legend()
plt.show()
</pre>


    
![png](lecture-21_files/lecture-21_5_0.png)
    


Two key properties of a Poisson random variable are

1) the variance is equal to the mean

$$ 
\begin{aligned} 
\mathrm{Var}(Y) &= \lim_{n \rightarrow \infty} \mathrm{Var}\left(\mathrm{Bin}\left(n,\frac{\lambda}{n}\right)\right) \\
&= \lim_{n \rightarrow \infty} n \frac{\lambda}{n} \left(1-\frac{\lambda}{n}\right) \\
&= \lambda  \lim_{n \rightarrow \infty} \left(1-\frac{\lambda}{n}\right) \\
&= \lambda
\end{aligned}
$$

2) if $Y_1$ and $Y_2$ are two independent Poisson random variables with means $\lambda_1$ and $\lambda_2$, then their sum is distributed as a Poisson with mean $\lambda_1+\lambda_2$

$$Y_1 + Y_2 \sim \mathrm{Poi}(\lambda_1 + \lambda_2)$$

<span id='section2'></span>
## 2. Demographic stochasticity
<hr>

We can now model demographic stochasticity with the Poisson distribution. 

Assume each individual in the population produces a Poisson number of surviving offspring with mean $\lambda$, independent of all other individuals in the populations, and then dies. Write the number of offspring for individual $i$ as $X_i\sim\mathrm{Poi}(\lambda)$. Let the current number of individuals be $n_t$. Then the population size in the next generation, $n_{t+1}$, is distributed like the sum of $n_t$ independent and identical Poisson's

$$\begin{aligned}
n_{t+1} &\sim \sum_{i=1}^{n_t}\mathrm{Poi}(\lambda)\\
&= \mathrm{Poi}\left(\sum_{i=1}^{n_t} \lambda\right) \\
&= \mathrm{Poi}(\lambda n_t) \\
\end{aligned}$$

This implies that the expected population size in the next generation is $\lambda n_t$, as is the variance. In Lab 11 we'll simulate this to get a better sense of the resulting dynamics.

<span id='section3'></span>
## 3. Extinction
<hr>

Let us now look at one property of this model -- the probability of extinction.

Consider a given individual at time $t$. Let $\eta$ be the probability this individual does not have descendants in the long-term, i.e., that its lineage goes extinct.

To find $\eta$ we note that the probability that this lineage goes extinct, $\eta$, is the probability that this individual has $k$ surviving offspring (for all values of $k$) and all of those offspring lineages go extinct (with probability $\eta^k$). Given the probability of having $k$ surviving offspring is $\lambda^k e^{-\lambda}/k!$, this implies

$$
\begin{aligned}
\eta 
&= \sum_{k=0}^{\infty} \frac{e^{-\lambda}\lambda^k}{k!} \eta^k \\
&= e^{-\lambda} \sum_{k=0}^{\infty} \frac{(\eta \lambda)^k}{k!} \\
&= e^{-\lambda}e^{\lambda \eta} \\
&= e^{-\lambda (1-\eta)}
\end{aligned}
$$

One solution is $\eta=1$, certain extinction. But there can be a second biologically valid solution (between 0 and 1), meaning that extinction is not certain. This occurs when the mean number of surviving offspring is greater than one, $\lambda>1$, as shown in plot below (the solutions are where the curve intersects the 1:1 line).


<pre data-executable="true" data-language="python">
xs = np.linspace(0,1,100)
fig,ax=plt.subplots()

lam = 2

ax.plot(xs, xs)
ax.plot(xs, [np.exp(-lam*(1-x)) for x in xs])

ax.set_xlabel(r'$\eta$')
ax.set_ylabel(r'$e^{-\lambda (1-\eta)}$')

plt.show()
</pre>


    
![png](lecture-21_files/lecture-21_11_0.png)
    


Above we calculated the probability a single lineage goes extinct, $\eta$. From this, the probability that the entire population goes extinct is the probability that all $n_t$ lineages go extinct, $\eta^{n_t}$. 

The two major conclusions from this are

1) when the mean number of surviving offspring per parent is $\lambda=1$ the deterministic model predicts a constant population size but the stochastic model says that extinction is certain.

2) when the mean number of surviving offspring per parent is $\lambda>1$ the deterministic model predicts exponential growth but the stochastic model says there is still some non-zero probability of extinction.

<span id='section4'></span>
## 4. Establishment
<hr>

Before moving on, the above result about extinction is mathematically identical to a classic result in population genetics concerning the establishment of a beneficial allele.

Consider a population at equilibrium such that the mean number of offspring per parent is 1. And now consider a beneficial allele that tends to have more offspring, $\lambda>1$. Then we know from above that there is some non-zero probability this lineage does not go extinct, $\eta<1$. Writing the above equation about extinction in terms of the **establishment probability**, $p=1-\eta$, we have

$$\begin{aligned}
\eta &= e^{-\lambda (1-\eta)}\\
1 - p &= e^{-\lambda p}\\
p &= 1 - e^{-\lambda p}
\end{aligned}$$

Now assume that the beneficial allele increases the number of offspring only slightly, so that $\lambda=1+s$ with $s$ small. We can then solve for $p$ explicitly using a Taylor series expansion of $e^{-\lambda p}=e^{-(1+s) p}$ around $s=0$

$$ 
\begin{aligned}
p &= 1 - e^{-(1+s)p} \\
p &\approx 1 - \left(1 - (1+s)p + \frac{(1+s)^2p^2}{2}\right) \\
p &\approx (1+s)p - \frac{(1+s)^2p^2}{2} \\
1 &\approx 1 + s - \frac{(1+s)^2p}{2} \\
0 &\approx s - \frac{(1+s)^2p}{2} \\
p &\approx \frac{2s}{(1+s)^2}\\
p &\approx 2s
\end{aligned}
$$

This says that the probability a weakly beneficial allele establishes (i.e., is not lost by genetic drift) is roughly twice its selective advantage.
