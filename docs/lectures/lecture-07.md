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

# Lecture 7: Equilibria

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

1. [Equilibria](#section1)
2. [Exponential growth](#section2)
3. [Logistic growth](#section3)
4. [Haploid selection](#section4)
5. [Diploid selection](#section5)
6. [Summary](#section6)

<span id='section1'></span>
## 1. Equilibria
<hr>

An **equilibrium** is any state of a system which tends to persist unchanged over time.

For *discrete-time* models, the equilibria are defined as those values of the variables where no changes occur from one time step to the next. 

For example, those values of allele frequency $p_t$ where

$$
\begin{aligned}
\Delta p &= 0\\
p_{t+1} - p_t &= 0\\
p_{t+1} &= p_t
\end{aligned}
$$

Similarly, for *continuous-time* models, the equilibria are defined as those values of the variables for which the rate of change in the variables equals zero. 

For example, those values of allele frequency $p$ where

$$
\frac{\mathrm{d}p}{\mathrm{d}t} = 0
$$


What are the equilibria for the following models?

| Model | Discrete time | Continous time |
| ----- | ------------- | -------------- |
| Exponential growth | $n_{t+1} = R n_t$ | $\frac{\mathrm{d}n}{\mathrm{d}t} = r n$ |
| Logistic growth | $n_{t+1} = (1 + r(1 - \frac{n_t}{K}))n_t$ | $\frac{\mathrm{d}n}{\mathrm{d}t} = r(1 - \frac{n}{K})n$ |
| Haploid selection | $p_{t+1} = \frac{W_A p_t}{W_A p_t + W_a q_t}$ | $\frac{\mathrm{d}p}{\mathrm{d}t} = s p(1-p)$ |
| Diploid selection | $p_{t+1} = \frac{p_t^2W_{AA} + p_t q_tW_{Aa}}{p_t^2W_{AA} + 2 p_t q_tW_{Aa} +  q_t^2W_{aa}}$ | Not derived |

<span id='section2'></span>
## 2. Exponential growth
<hr>


Here, we will solve for the equilibria in both the discrete- and continuous-time exponential-growth models.

### Discrete time

$$
n_{t+1} = Rn_t
$$

Set $n_{t+1} = n_t = \hat n$ and solve for $\hat{n}$

$$
\begin{aligned}
\hat n &= R\hat n\\
\hat n &= 0
\end{aligned}
$$

### Continuous time

$$
\frac{\mathrm{d}n}{\mathrm{d}t} = r n
$$

Set $\mathrm{d}n/\mathrm{d}t = 0$ and $n = \hat n$ and solve for $\hat n$

$$
\begin{aligned}
0 &= r \hat{n}\\
\hat n &= 0
\end{aligned}
$$

So the only equilibrium in both discrete- and continuous-time exponential growth is extinction, $\hat{n}=0$.

!!! note "Special case of parameters"

    Notice above that $R=1$ and $r=0$ also satisfy the conditions for an equilibrium. These are called **special cases of parameters**. Here this refers to the case where individuals perfectly replace themselves so that the population remains constant from *any* starting value of $n$.

<span id='section3'></span>

## 3. Logistic growth
<hr>

Here, we will solve for the equilibria in both the discrete- and continuous-time logistic-growth models.

### Discrete time

$$
n_{t+1} = \left(1 + r\left(1 - \frac{n_t}{K}\right)\right)n_t
$$

As above, we substitute $n_{t+1} = n_t = \hat n$ and want to solve for $\hat{n}$.

$$
\hat n = \left(1 + r\left(1 - \frac{\hat n}{K}\right)\right)\hat{n} 
$$

Notice that one equilibrium is $\hat n = 0$. However, this isn't the only equilibrium because dividing both sides by $\hat n$ results in

$$
\begin{aligned}
1 &= 1 + r\left(1 - \frac{\hat n}{K}\right)\\
0 &= r\left(1 - \frac{\hat n}{K}\right)
\end{aligned}
$$

Here we have a special case of parameters, $r=0$, or

$$
\begin{aligned}
0 &= 1 - \frac{\hat n}{K}\\
\hat n &= K
\end{aligned}
$$

There are therefore two equilibria: extinction, $\hat{n}=0$, or carrying capacity, $\hat{n}=K$.

### Continuous time

$$
\frac{\mathrm{d}n}{\mathrm{d}t} = r \left(1 - \frac{n}{K}\right)n
$$

We set $\mathrm{d}n/\mathrm{d}t=0$ and $n=\hat{n}$

$$
0 =  r \left(1 - \frac{\hat n}{K}\right)\hat{n}
$$

which is the same equation we had above in discrete-time, so the equilibria ($\hat n = 0,K$) and the special case of parameters ($r = 0$) are also the same.

<span id='section4'></span>
## 4. Haploid selection
<hr>

Here, we will solve for the equilibria in both the discrete- and continuous-time haploid-selection models.

### Discrete time

$$
p_{t+1} = \frac{p_tW_A}{p_tW_A + q_tW_a}
$$

Replace $p_{t+1}$ and $p_t$ with $\hat p$ and replace $q_t$ with $\hat q$ and solve for $\hat p$ and $\hat q$

$$
\begin{aligned}
\hat{p} &= \frac{\hat p W_A}{\hat p W_A + \hat q W_a}\\
\end{aligned}
$$

We first see that $\hat{p}=0$ is an equilibrium. But there is more, since dividing by $\hat p$ gives

$$
\begin{aligned}
1 &= \frac{W_A}{\hat p W_A + \hat q W_a}\\
\hat p W_A + \hat q W_a &= W_A\\
\hat q W_a &= (1-\hat p) W_A
\end{aligned}
$$

At this point we use $q=1-p$ to write this in terms of $p$ only

$$
(1-\hat p) W_a = (1-\hat p) W_A
$$

So $\hat p =1$ is another equilibrium. 

And finally, dividing by $(1-\hat p)$ gives a special case of parameters, $W_A=W_a$.

To summarize, the allele frequency will not change from one generation to the next in our discrete-time haploid-selection model when

- $\hat p = 0 \Longrightarrow$ the population is "fixed" for the $a$ allele
- $\hat p = 1 \Longrightarrow$ the population is fixed for the $A$ allele
- $W_A = W_a \Longrightarrow$ the two alleles have equal fitness ("neutrality")

### Continuous time

$$
\frac{\mathrm{d}p}{\mathrm{d}t} = sp(1-p)
$$

In the continuous-time model, we set the derivative equal to zero and $p=\hat{p}$

$$
\begin{aligned}
0 &= s\hat p(1 -\hat p)
\end{aligned}
$$

And we again find the same equilibria ($\hat p=0,1$) and special case of parameters ($s=0$, i.e., neutrality).

<span id='section5'></span>
## 5. Diploid selection
<hr>

### Discrete time

Here, we will solve for the equilibria in the discrete-time diploid-selection model

$$
p_{t+1} = \frac{p_t^2 W_{AA} + p_t  q_t W_{Aa}}{p_t^2 W_{AA} + 2 p_t q_t W_{Aa} + q_t^2 W_{aa}}
$$

We replace $p_{t+1}$ and $p_t$ with $\hat p$ and $q_t$ with $\hat{q}$ and solve for $\hat p$ and $\hat q$

$$
\begin{aligned}
\hat{p} &= \frac{\hat{p}^2 W_{AA} + \hat{p} \hat{q} W_{Aa}}{\hat{p}^2 W_{AA} + 2 \hat{p} \hat{q} W_{Aa} + \hat{q}^2 W_{aa}}\\
\hat{p} &= \frac{\hat{p}(\hat p W_{AA} + \hat{q} W_{Aa})}{\hat{p}^2 W_{AA} + 2 \hat{p} \hat{q} W_{Aa} + \hat{q}^2 W_{aa}}
\end{aligned}
$$

We see that $\hat{p}=0$ is one equilibrium. Moving on, dividing by $\hat p$ gives

$$
\begin{aligned}
1 &= \frac{\hat p W_{AA} + \hat{q} W_{Aa}}{\hat{p}^2 W_{AA} + 2 \hat{p} \hat{q} W_{Aa} + \hat{q}^2 W_{aa}}\\
\hat{p}^2 W_{AA} + 2 \hat{p} \hat{q} W_{Aa} + \hat{q}^2 W_{aa} &= \hat p W_{AA} + \hat{q} W_{Aa}\\
0 &= (\hat{p} - \hat{p}^2) W_{AA} + (\hat{q} - 2 \hat{p} \hat{q}) W_{Aa} - \hat{q}^2 W_{aa}\\
0 &= \hat{p}(1 - \hat{p}) W_{AA} + \hat{q}(1 - 2 \hat{p}) W_{Aa} - \hat{q}^2 W_{aa}\\
0 &= \hat{q}(\hat{p} W_{AA} + (1 - 2 \hat{p}) W_{Aa} - \hat{q} W_{aa})
\end{aligned}
$$

And so $\hat{q}=0\implies\hat{p}=1$ is another equilibrium. Dividing by $\hat{q}$ and putting everything in terms of $p$ we have

$$
\begin{aligned}
0 &= \hat{p} W_{AA} + (1 - 2 \hat{p}) W_{Aa} - \hat{q} W_{aa}\\
0 &= \hat{p} W_{AA} + (1 - 2 \hat{p}) W_{Aa} - (1 - \hat{p}) W_{aa}\\
0 &= \hat{p}(W_{AA} -2W_{Aa} + W_{aa}) + W_{Aa} - W_{aa}\\
W_{aa} - W_{Aa} &= \hat p(W_{AA} -2W_{Aa} + W_{aa})\\
\frac{W_{Aa} - W_{aa}}{2W_{Aa} - W_{AA} - W_{aa}} &= \hat p\\
\end{aligned}
$$

We therefore have *three* equilibria under diploid selection: $\hat{p}=0,\frac{W_{Aa} - W_{aa}}{2W_{Aa} - W_{AA} - W_{aa}},1$.

Since a frequency is bounded between 0 and 1, we must have $0 \leq p \leq 1$. We therefore call $\hat{p}=0$ and $\hat{p}=1$ **boundary equilibria**.
  
These bounds also imply the third equilibrium is only **biologically valid** when 

$$
0 \leq \frac{W_{Aa} - W_{aa}}{2 W_{Aa} -W_{AA} - W_{aa}} \leq 1
$$

When $W_{Aa} = W_{aa}$ this equilibrium reduces to $\hat{p}=0$ and when $W_{Aa} = W_{AA}$ this reduces to $\hat{p}=1$ (check this for yourself).

The third equilibrium will be an **internal equilibrium**, representing a population with both $A$ and $a$ alleles, when

$$
0 < \frac{W_{Aa} - W_{aa}}{2 W_{Aa} -W_{AA} - W_{aa}} < 1
$$

The equilibrium is positive when the numerator and denominator have the same sign (i.e., are both positive or both negative).

Let's split this into two "cases". Case A will have a positive numerator, $W_{Aa} > W_{aa}$, and Case B will have a negative numerator, $W_{Aa} < W_{aa}$.

So, in Case A, the equilibrium is positive when the denominator is positive, $2 W_{Aa} - W_{AA} - W_{aa} > 0$.

While in case B the equilibrium is positive when the denominator is negative, $2 W_{Aa} - W_{AA} - W_{aa} < 0$.

Now we can rearrange the equilibrium to show that it is less than 1 when

$$
\begin{aligned}
\frac{W_{Aa} - W_{aa}}{2 W_{Aa} -W_{AA} - W_{aa}} &< 1\\
\frac{W_{Aa} - W_{aa}}{2 W_{Aa} -W_{AA} - W_{aa}} - 1 &< 0\\
\frac{W_{Aa} - W_{aa} - (2 W_{Aa} -W_{AA} - W_{aa})}{2 W_{Aa} -W_{AA} - W_{aa}} &< 0\\
\frac{W_{AA} - W_{Aa}}{2 W_{Aa} -W_{AA} - W_{aa}} &< 0\\
\frac{W_{Aa} - W_{AA}}{2 W_{Aa} -W_{AA} - W_{aa}} &> 0\\
\end{aligned}
$$
  
Again we need the numerator and denominator to have the same sign for this inequality to hold.
  
In case A, where we've said that denominator is positive, this means we also need the numerator to be positive, $W_{Aa} > W_{AA}$.

While in case B we said that the denominator is negative, so we also need the numerator to be negative, $W_{Aa} < W_{AA}$.
  
Putting this all together, there is a biologically-relevant internal equilibrium when either

- Case A: $W_{Aa} > W_{aa}$ and $W_{Aa} > W_{AA}$ (which ensures $2 W_{Aa} - W_{AA} - W_{aa} > 0$; go ahead and check!)
- Case B: $W_{Aa} < W_{aa}$ and $W_{Aa} < W_{AA}$  (which ensures $2 W_{Aa} - W_{AA} - W_{aa} < 0$)

Case A therefore represents "heterozygote advantage", $W_{AA} < W_{Aa} > W_{aa}$, while Case B represents "heterozygote disadvantage", $W_{AA} > W_{Aa} < W_{aa}$. 

<span id='section6'></span>
## 6. Summary
<hr>

In summary, the equilibria for the models we have looked at are:

| Model | Discrete-time equilibria | Continuous-time equilibria |
| ----- | ------------------------ | -------------------------- |
| Exponential growth | $\hat n = 0$ | $\hat n = 0$ |
| Logistic growth | $\hat n = 0, \hat n = K$ | $\hat n = 0, \hat n = K$ |
| Haploid selection | $\hat p = 0, \hat p = 1$ | $\hat p = 0, \hat p = 1$
| Diploid selection | $\hat p = 0, \hat p = 1, \hat p = \frac{W_{Aa} - W_{aa}}{2W_{Aa} - W_{AA} - W_{AA}}$ | Not derived |

Make sure that you understand how to determine equilibria in discrete- and continuous-time and can derive the equilibria of the models above on your own.
