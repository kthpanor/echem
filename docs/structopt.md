(sec:structopt-methods)=
# Optimization Methods

The first class of special points on the PES that we will discuss are local energy minima. These correspond to equilibrium molecular structures and are characterized by a vanishing first order energy derivative combined with a Hessian matrix which has only positive eigenvalues. The procedure to determine a local minimum (i.e. finding the coordinates that minimize the energy) is called a structure optimization.  

```{figure} /img/pes/PES.svg
:scale: 100%
```
In practical terms, the ingredients to perform a geometry optimization include: (1) the initial molecular coordinates, (2) a choice of coordinate system, (3) the energy at a specific geometry $E(\boldsymbol{\xi})$, (4) the gradient $\nabla E(\boldsymbol{\xi})$, (5) the Hessian, and (6) a procedure to update the coordinates and Hessian and move on the potential energy surface towards lower energy. 

Having addressed [the issue of coordinate system](coord), the remaining question is what procedure to use to move along the potential energy surface and arrive at a local energy minimum. There are several iterative methods to do this, some of which need only information on the energy gradient (e.g. gradient-descent, conjugate gradient), while others take into account also the Hessian (Newton--Raphson, quasi-Newton). For a detailed review of minimization techniques, see {cite}`Snyman2005`.

## Gradient descent

The simplest optimization procedure is to repeatedly take a step in the direction opposite to the local gradient:
%
\begin{equation}
\mathbf{x}_{i+1} = \mathbf{x}_i - k_i\nabla E(\mathbf{x}_i) \label{eq:grad-descent} \,, 
\end{equation}
where by $\mathbf{x}_{i+1}$ we denote the new coordinates (in a generic coordinate system -- either Cartesian or internal coordinates), $\mathbf{x}_i$ are the coordinates at the previous step $i$, $\nabla E$ is the energy gradient and $k_i$ is the step size. The step size can either be kept constant, or adjusted at each iteration, e.g. by the line search procedure.

The gradient-descent method is simple to implement and is guaranteed to converge, but has the disadvantage that it requires many steps and becomes slow when close to the minimum where the gradient is small. It always converges to a local minimum, given enough steps.  

## Conjugate gradient
An improved method over the gradient-descent approach is to use the "gradient history" (steps $i$ and $i-1$) to determine the coordinates at step $i+1$:
%
\begin{equation}
    \mathbf{x}_{i+1} = \mathbf{x}_i - k_i \mathbf{h}_i\,,\label{eq:conjugate_gradient}
\end{equation}
with $\mathbf{h}_i = \nabla E (\mathbf{x}_i)+\gamma_i\mathbf{h}_{i-1}$. The function $\gamma_i$ contains gradient information from steps $i$ and $i-1$ and can be defined in different ways. For example, in the Fletcher-Reeves conjugate gradient method: 
%
\begin{equation}
    \gamma_i = \frac{|\nabla E(\mathbf{x}_i)|^2}{|\nabla E(\mathbf{x}_{i-1})|^2}\,.
\end{equation}

## Newton--Raphson and quasi-Newton
The next step in the hierarchy of minimization methods is to use both the first and second order energy derivatives (i.e. gradient $ \nabla E$ and Hessian $\mathbf{H}$) in determining the next step in conformation space. This is based on a quadratic approximation for the local shape of the PES:
%
\begin{equation}
    E(\mathbf{x}+\Delta\mathbf{x}) \approx E(\mathbf{x}) + \nabla E(\mathbf{x})\Delta\mathbf{x} + \frac{1}{2}\Delta\mathbf{x}^\mathrm{T}\mathbf{H}\Delta\mathbf{x} \label{eq:quadratic_approx_PES} \,,
\end{equation}
%
where $\Delta \mathbf{x}$ is the Newton step used to update the coordinates: 
%
\begin{flalign}
    \Delta \mathbf{x} &= -\mathbf{H}^{-1}\nabla E(\mathbf{x})\label{eq:Newton_step}\,,\\
    \mathbf{x}_\mathrm{i+1} &= \mathbf{x}_\mathrm{i}+\Delta \mathbf{x}\,.
\end{flalign}
 %
 When redundant internal coordinates are used, it is important to ensure that the displacements are only performed in the non-redundant part of the internal coordinate space. This is achieved by applying a projector $\mathbf{P}=\mathbf{G}^{-}\mathbf{G}$ to the gradient and Hessian before constructing the Newton step {cite}`orcamanual`:
 %
 \begin{flalign}
 \tilde{\mathbf{g}}_q &= \mathbf{P}\nabla E(\mathbf{q}) \, , \\
 \tilde{\mathbf{H}}_q &= \mathbf{P}\mathbf{H}_q\mathbf{P}+\alpha(1-\mathbf{P}) \,,
 \end{flalign}
%
where $\alpha$ is an arbitrary large value (e.g. 1000).
%and is used to set the redundant part of 
%\textcolor{red}{I haven't perfectly understood what $\alpha$ does yet. The exact quote from the ORCA manual is "The second term for H sets the matrix elements of the redundant part of the internal coordinate space to very large values ($\alpha$ = 1000)"}

As evident from Eq.~\eqref{eq:Newton_step}, the second order derivatives of the energy with respect to nuclear displacements are required to generate the Newton step. The direct computation of these derivatives (which compose the Hessian matrix) is quite expensive, but good \textit{approximations} for the Hessian can be constructed using  the gradient history. For example, the Broyden-Fletcher-Goldfarb-Shanno (BFGS, used by \code{geomeTRIC}) approach uses the relation:
%
\begin{flalign}
\mathbf{H}_{i+1} &= \mathbf{H}_i + \frac{\mathbf{g}^{\textcolor{white}{.}}_i\mathbf{g}_i^\mathrm{T}}{\mathbf{g}_i^\mathrm{T}\mathbf{s}^{\textcolor{white}{.}}_i} - \frac{\mathbf{H}_i\mathbf{s}_i\left(\mathbf{H}_i\mathbf{s}_i\right)_{\textcolor{white}{i}}^\mathrm{T}}{\mathbf{s}_i^\mathrm{T}\mathbf{H}_i^{\textcolor{white}{.}}\mathbf{s}_i^{\textcolor{white}{.}}} \, ,\label{eq:BFGS}
\end{flalign}
with,
\begin{flalign}
\mathbf{g}_i &= \nabla E(\mathbf{x}_{i+1}) - \nabla E(\mathbf{x}_i) \, \\
\mathbf{s}_i &= \mathbf{x}_{i+1} - \mathbf{x}_{i} \,,
\end{flalign}
to update the Hessian at step $i+1$ using the Hessian at step $i$ and information about the gradient at the current and previous step.

Quasi-Newton methods, i.e., methods that use approximate Hessians, achieve a very quick convergence at the same computational cost as the gradient-descent method.  

```{figure} /img/pes/geom_opt_flowchart.svg
:scale: 100%
```
