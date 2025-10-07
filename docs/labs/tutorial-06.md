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

# Tutorial 6: Linear multivariate practice

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

Metastasis is a process by which cancer cells spread throughout the body. Sometimes cancer cells move via the bloodstream and become lodged in the capillaries of different organs. Some of these new cells then move across the capillary wall, where they initiate new tumours. We will construct a model for the dynamics of cancer cells lodged in the capillaries of an organ, $C$, and the number of cancer cells that have actually invaded that organ, $I$. Suppose that cells are lost from the capillaries by dislodgement or death at per capita rate $\delta_1$ and that they invade the organ from the capillaries at a per capita rate $\beta$. Once cells are in the organ they die at rate $\delta_2$, and the cancer cells replicate at a per capita rate $\rho$. This gives the differential equations

$$
\begin{aligned}
\frac{\mathrm{d}C}{\mathrm{d}t} &= -\delta_1 C - \beta C\\
\frac{\mathrm{d}I}{\mathrm{d}t} &= \beta C -\delta_2 I + \rho I.
\end{aligned}
$$

## Questions

Write this in matrix form, $\frac{\mathrm{d}\vec{n}}{\mathrm{d}t} = \mathbf{M}\vec{n}$, specifying what $\mathbf{M}$ and $\vec{n}$ are.

What is the equilibrium of this model?

Calculate the eigenvalues of $\mathbf{M}$. When is the equilibrium stable? (Based on the description of the model we can assume all parameters are positive.) Explain what this means biologically.

Assuming the system is unstable, calculate the right eigenvector associated with the leading of $\mathbf{M}$. Where are most of the cells in the long term?
