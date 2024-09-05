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

# Lecture 1: Introduction

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

1. [Why use mathematical models in ecology and evolution?](#section1)
2. [Syllabus](#section2)

<span id='section1'></span>
## 1. Why use mathematical models in ecology and evolution?
<hr>

Mathematics permeates ecology and evolution, from simple back-of-the-envelope calculations to the development of sophisticated mathematical models. This is because mathematics is a unique tool that can take what we already know (or assume to be true) and rigorously lead us to the logical conclusions.

To see this (while also introducing you to the kind of models I work with), let's look at two examples.

### Example 1: HIV

See sections 1.2-1.4 of the text for more info.

Human immunodeficiency virus (HIV) is, as the name suggests, a virus that infects humans and causes acquired immunodeficiency syndrome (AIDS). This is a terrible diesase that has caused over 20 million deaths worldwide. It is transmitted via bodily fluids. Once inside the body HIV particles (virions) attach to a protein called CD4 on the cell membrane of helper T cells (part of our immune system) and others. Once attached, the virus inserts its RNA into the host cell, which is reverse transcribed into DNA and becomes part of the host genome. The host cell can remain in the 'latently infected' state for some time. When the viral DNA is eventually transcribed by the host cell, starting an 'active infection', hundreds of new virions are produced, often killing the host cell. These new virions then go on to infect other CD4+ cells. A large decline in the number of helper T cells is what causes AIDS. This is because helper T cells bind to viruses and secrete chemical signals to stimulate the immune system. So without helper T cells the immune system is very comprimised.

<center>

![](lecture-01-img/hiv.png){ width=50% }

</center>

<center><sup>HIV and T cells. Source: bio.libretexts.org</sup></center>

Based on this, once infected with HIV (and without treatment) we expect the number of virions to rapidly increase and the number of helper T cells to decline. This is generally what is observed. However, the number of virions then tends to steeply decline. Why might this be?

Here are two hypotheses:

1. the immune system recognizes HIV and suppresses it
2. the decline in helper T cells prevents HIV from replicating

To decide whether the second hypothesis is valid [Phillips (1996)](https://www.science.org/doi/abs/10.1126/science.271.5248.497) built a mathematical model describing the rate of change in the number of CD4+ cells and free virions. 

See the note below for the details of the model if you are interested (but don't worry if you don't understand much of this yet!).

??? note "Details of the Philips (1996) model"

    See Box 2.4 in the text for more details.

    Philips built a dynamical model of the infection process using 4 differential equations:

    $\frac{dR}{dt} = \Gamma \tau - \mu R - \beta V R$

    $\frac{dL}{dt} = p \beta V R - \mu L - \alpha L$

    $\frac{dE}{dt} = (1 - p) \beta V R + \alpha L - \delta E$

    $\frac{dV}{dt} = \pi E - \sigma V$

    These equations describe the rate of change in the number of susceptible CD4+ cells (R), latently infected cells (L), actively infected cells (E), and free virions (V). These are the 4 "variables" of the model, as their values change over time. The remainder of the symbols represent "parameters", whose values do not change (they are constants). The meaning of the parameters and the values used by Phillips (1996) are shown in the table below.

    | Symbol | Description | Value (units/day) |
    | ------ | ----------- | ----------------- |
    | $\Gamma$ | Rate that CD4+ cells are produced | 1.36 |
    | $\tau$ | Proportion of CD4+ cells that are susceptible to attack | 0.2 |
    | $\mu$ | HIV-independent death rate of susceptible CD4+ cells | 1.36 x 10^-3 |
    | $\beta$ | Rate of CD4+ cell infection per HIV virion | 0.00027 | 
    | $p$ | Proportion of newly infected cells becoming latently infected | 0.1 | 
    | $\alpha$ | Activation rate of latently infected cells | 3.6 x 10^-2 |
    | $\delta$ | Death rate of actively infected cells | 0.33 |
    | $\pi$ | Rate of production of virions by actively infected cells | 100 |
    | $\sigma$ | Removal rate of cell-free virus | 2 |
    
    Given these parameter descriptions, can you "read" the differential equations above? For example, what does $\Gamma \tau$ represent, biologically?

    Using these equations and parameter values Philips (1996) asked how the number of CD4+ cells (R + L + E) and free virions (V) changed over time following infection. 

To see how Philips' model behaves, activate the kernel at the top of the page and run the code below


<pre data-executable="true" data-language="python">
import numpy as np #numerical tools
import matplotlib.pyplot as plt #plotting tools

# define a function to iterate through the model
def philips_model(R=200, L=0, E=0, V=4e-5, days=120, steps=120): #default initial conditions and timesteps
    #choose parameter values
    gamma, mu, tau, beta, p, alpha, sigma, delta, pi = [1.36, 1.36e-3, 0.2, 0.00027, 0.1, 3.6e-2, 2, 0.33, 100] 
    #update variables
    record = []
    ts = np.linspace(0, days, steps) #times we record population sizes at
    for t in ts:  #looping over timesteps      
        dRdt = gamma * tau - mu * R - beta * V * R
        R += dRdt #uninfected T cells
        dLdt = p * beta * V * R - mu * L - alpha * L
        L += dLdt #latently infected T cells
        dEdt = (1-p) * beta * V * R + alpha * L - delta * E
        E += dEdt #actively infected T cells
        dVdt = pi * E - sigma * V
        V += dVdt #free virions
        record += ([[R + L + E, V]]) #total T cells, virions
    return ts, np.array(record)

# iterate and record population sizes
ts, ns = philips_model() #times and population sizes

# initialize plot
fig, left_ax = plt.subplots()
right_ax = left_ax.twinx()
fig.set_size_inches(8,4)

# plot T cell data on left axis
left_ax.plot(ts, ns[:,0], label='R');
left_ax.set_ylabel('CD4+ cells (R + E + L)', fontsize=14)
left_ax.legend(frameon=False, bbox_to_anchor=(0.99, 0.99))
left_ax.set_xlabel('Days from infection', fontsize=14)

# plot virion data on right axis
right_ax.plot(ts, ns[:,1], label='V', color=plt.cm.tab10(1));
right_ax.set_ylabel('Free virions V', fontsize=14)
right_ax.legend(frameon=False, bbox_to_anchor=(0.99, 0.85))

# format plot and show it
plt.yscale('log')
fig.tight_layout()
plt.show()
</pre>


    
![png](lecture-01_files/lecture-01_5_0.png)
    


Compare to Figure 1.3 and Figure 1.4 in the text. We see the initial increase in virions (orange) and decline in CD4+ cells (blue), followed by a decline virions (note the log scale on the right axis -- this is a big decline, from over 1000 to about 10). Because this model does not include an immune response against the virions but still exhibits the decline in virions, we conclude that the second hypothesis is valid, that it is theoretically plausible that the decline in virions is due to a lack of CD4+ cells to infect. A few years later this hypothesis was empirically tested and validated -- a nice example of theory guiding science.

Feel free to play around with the code above, changing parameter values or even the structure of the model. Do the dynamics change as you expected?

### Example 2: Extreme Events

A second example is a model that I helped Dr. Kelsey Lyberger (then a PhD student at UC Davis) with in [Lyberger et al 2021](https://www.biorxiv.org/content/10.1101/2020.04.02.014951v2).

<center>

![](lecture-01-img/kelsey.png){ width=50% }

</center>

<center><sup>Kelsey Lyberger, doing <i>Daphnia</i> fieldwork I suppose.</sup></center>

Kelsey was interested in how populations respond to extreme climatic events, like lizards to hurricanes. It has long been clear that such events can impact the size of a population, e.g., by causing extra mortality, and may in fact put populations at risk of extinction. More recently it has become apparent that extreme events can also impose strong natural selection, and that populations can quickly adapt to the new environment. Some examples include:

- [Hurricanes select on lizard limbs and toe pads](https://www.pnas.org/doi/10.1073/pnas.2000801117)
- [Ice-storms select on sparrow body size](https://www.jstor.org/stable/2406980)
- [Droughts select on Darwin finch beaks](https://www.science.org/doi/10.1126/science.aad8786)
- [Droughts select on flowering time in *Brassica*](https://nph.onlinelibrary.wiley.com/doi/10.1111/j.1469-8137.2010.03603.x)

<center>

![](lecture-01-img/lizard.png)
    
</center>

<center><sup>Lizard being blown off perch by a leaf blower. Source: colindonihue.com</sup></center>

Now, how should such rapid adaptive evolution impact population size? This is the question Kelsey set out to answer with a mathematical model.

??? note "Details of Kelsey's model"

    Kelsey assumed each individual has a quantitative genetic trait, such as lizard limb length, that is determined by many alleles of small effect plus some environmental noise. Fitness is assumed to be a bell-shaped function of the difference between the optimum trait value, $\theta$, and an individual's trait value, $z$, which we write as $W(\theta - z)$. 
    
    The change in population size, $N$, from one generation to rate the next is the current population size times mean fitness, $\overline{W}(\theta - \bar{z})$, the latter depending on the population mean trait value, $\bar{z}$. This gives

    $\Delta N = N \overline{W}(\theta - \bar{z}).$
    
    The change in the mean trait value from one generation to the next is roughly the product of genetic variance in the trait, $V_g$, and the strength of selection. The strength of selection is defined as the derivative of the natural logarithm of population mean fitness with respect to the mean trait value, $\mathrm{d}\ln\overline{W}(\theta - \bar{z})/\mathrm{d}\bar{z}$. This gives
    
    $\Delta \bar{z} = V_g \frac{\mathrm{d}}{\mathrm{d}\bar{z}}\ln\overline{W}(\theta - \bar{z}).$
    
    Together, these coupled recursion equations, $\Delta N$ and $\Delta \bar{z}$, can be used to describe how evolution affects population size under an extreme event, which is modeled as a sudden but temporary change in the optimum phenotype, $\theta$.

Below is a stochastic simulation much like that used by Kelsey. With an activated kernel, run the code below to create a plot very similar to Figure 1 in Lyberger et al. (this may take a minute).


<pre data-executable="true" data-language="python">
import numpy as np
import matplotlib.pyplot as plt

def lyberger_model(Vo=0.75, Ve=0, event_duration=1, seed=0, other_parameters=[120, 500, 1, 2, 100, 0, 2.5]):
    
    # unpack parameters
    generations, K, w, lmbda, event_time, initial_theta_t, dtheta_t = other_parameters
    
    # initialize population
    members = np.random.normal(initial_theta_t, Vo**0.5, K) #K individuals with trait values distributed around the optimum
    
    # run simulations
    np.random.seed(seed)
    population_size, mean_breeding_value = [], []
    for g in range(generations):

        # optimum trait value
        if g in np.arange(event_time, event_time + event_duration):
            theta_t = initial_theta_t + dtheta_t #extreme event
        else:
            theta_t = initial_theta_t #normal

        # viability selection
        prob_survival = np.array([np.exp(-(theta_t - z + np.random.normal(0, Ve))**2 / (2*w**2)) for z in members])
        survived = np.array([True if p > np.random.uniform(0,1) else False for p in prob_survival])
        if len(survived) == 0:
            break
        members = members[survived]

        # random mating
        offspring = []
        for m in np.random.choice(members, len(members)):
            if len(offspring) > K: #if more than K offspring already then choose K and stop matings
                offspring = np.random.choice(offspring, K)
                break
            else:
                n_off = np.random.poisson(lmbda)
                mate = np.random.choice(members)
                offspring += [(m + mate)/2 for _ in range(n_off)] #give offspring mean parental value

        # sample new trait values for offspring
        offspring = np.array(offspring)
        offspring = np.random.normal(offspring, Vo) #add segregation variance

        # record statistics
        population_size.append(len(offspring))
        mean_breeding_value.append(np.mean(offspring))

        # update for next generation
        members = offspring
    
    return (
        np.arange(0, generations-event_time+1), #times
        np.array(population_size[event_time-1:]), #population size
        np.array(mean_breeding_value[event_time-1:]) #mean trait value
    )

# initialize plot
fig, ax = plt.subplots(2, sharex=True)
fig.set_size_inches(8,6)

# length of extreme event
event_duration = 1

# run 10 simulations per segregation (Vo) and environment variance (Ve) parameter combination
for Vo, Ve, c, lab in [[0.75, 0, 'black', 'with evolution'], [0, 1, 'red', 'without evolution']]:
    # plot simulations
    simulations = np.array([lyberger_model(Vo=Vo, Ve=Ve, event_duration=event_duration, seed=s) for s in range(10)])
    ax[0].plot(simulations[:,0].T, simulations[:,1].T, color=c, alpha=0.3);
    ax[1].plot(simulations[:,0].T, simulations[:,2].T, color=c, alpha=0.3);
    
    # hack together only one instance of the legend
    ax[0].plot([np.min(simulations[:,0].T)], [np.min(simulations[:,1].T)], alpha=1,
                  label = lab, color=c)
    ax[1].plot([np.min(simulations[:,0].T)], [np.min(simulations[:,2].T)], alpha=1,
                  label = lab, color=c)

# add environmental event duration
ax[0].fill_between([0,event_duration], y1=500, alpha=0.2)
ax[1].fill_between([0,event_duration], y1=-0.2, y2=2, alpha=0.2)

# add labels
ax[0].set_ylabel('Population size', fontsize=12)
ax[1].set_ylabel('Mean trait value', fontsize=12)
ax[1].set_xlabel('Generation', fontsize=14)

# add legend
plt.legend(frameon=False)
plt.show()
</pre>


    
![png](lecture-01_files/lecture-01_13_0.png)
    


The key result, that you can see in the plot above, is that when extreme events are short, adaptive evolution (black lines) can paradoxically *reduce* population size (relative to the red lines, where there is no evolution). The reason for this is that, while during the extreme event (shaded section) evolution is adaptive, once the extreme event ends the population finds itself maladapted to the original environment. Adaptive evolution can therefore hamper population persistence, and this is an important thing to keep in mind when documenting rapid adaptive evolution in response to extreme events -- it is not necessarily a good thing for the species (or our conservation goals).

<span id='section2'></span>
## 2. Syllabus
<hr>

OK, now that we've gone over some motivating examples of modeling in ecology and evolution, let's take a look at how we're going to learn to become modelers in this course. Point your browser over to the [syllabus](../../syllabus/general_info) and read each of the pages there.
