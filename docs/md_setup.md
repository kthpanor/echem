<!-- #region -->
# Ensembles and time integration

Workflow of a simulation.
- Figure


## Choice of ensemble

When simulating a system using MD, the ensemble of that system must be considered, i.e., the states in which it is possible to find the system and how probable these are.


### Microcanonical ensemble ($\textit{NVE}$)

Of the ensembles, the microcanonical is the simplest to simulate. Also known as the $\textit{NVE}$ ensemble, it represents a completely closed system, in which the number of particles, $N$, volume, $V$, and total energy, $E$, are kept constant. In simulation terms, this just means that the integration is left on its own from the starting position, with energy being converted back and forth between potential and kinetic form while the total sum remains constant. The main problem for $\textit{NVE}$ simulations is that it is hard to maintain energy conservation when using numerical integration methods. As can be seen in the examples in the previous section, the accumulation of numerical errors causes a drift in energy proportional to the time step.


### Canonical ensemble ($\textit{NVT}$)

An ensemble that is of more practical use for realistic simulations is the canonical, or $\textit{NVT}$, ensemble. In this case, the temperature is kept constant instead of the energy, emulating a system in thermal equilibrium with a heat bath. For the simulation, this means two things: First, the average temperature of the system should be constant. Note that this does not mean that the kinetic energy, the time average of which the temperature is proportional to, is constant at all times, just that it varies around the correct value. Second, how it varies should be determined by the fact that the velocities of the individual particles follow the Maxwell--Boltzmann distribution. The task of maintaining these two criteria in an MD simulation is done by a thermostat.


#### Thermostat

The most basic thermostat is a simple rescaling of the system velocities by a factor
%
\begin{equation*}
\lambda = \sqrt{\frac{T_0}{T}},
\end{equation*}
%
where $T$ is the current temperature and $T_0$ is the desired one. While this fulfils the criterion that the temperature should be constant, since the total kinetic energy is kept constant between steps, this also means that the second criterion is not fulfilled. An improved version of this method comes in the form of the Berendsen thermostat, [cite berendsen1984] which also includes a coupling parameter, $\tau$, determining the rate of heat transfer between the heat bath and the system. The rescaling factor then becomes
%
\begin{equation*}
\lambda = \sqrt{1 + \frac{\Delta t}{\tau} \left( \frac{T_0}{T} - 1 \right)}.
\end{equation*}
The coupling parameter needs to be chosen with care, as a $\tau$ that is equal to the time step, $\Delta t$, just reproduces the velocity rescaling thermostat and a too high $\tau$ makes the coupling too weak, with the limit $\tau \rightarrow \infty$ instead describing the microcanonical ensemble. So while the Berendsen thermostat does not describe a true canonical ensemble, there are coupling parameters between $\Delta t$ and infinity that give a suitably close approximation for larger systems. As the Berendsen thermostat is quite quick to converge to the correct temperature, it is often used for an initial equilibration, after which a true canonical thermostat is used. 
Examples of such thermostats are found in the Andersen [cite andersen1980] and Nosé--Hoover [cite nose1984a,nose1984b,hoover1985] thermostats. The Andersen thermostat takes a stochastic approach, replacing the velocities of a selection of atoms at each step by new velocities taken from the Maxwell--Boltzmann distribution, simulating random collisions. The frequency at which velocities are replaced becomes the coupling parameter of the thermostat, determining the speed at which it converges to the correct temperature. The Nosé--Hoover thermostat, on the other hand, introduces an artificial particle, with a mass and velocity, representing the heat bath. The extended system, containing both the real and artificial systems, is described by a microcanonical ensemble, but due to the coupling between the two, the real system becomes canonical.


### Isothermal-isobaric ensemble $(\textit{NPT})$

In the isothermal-isobaric, or $\textit{NPT}$, ensemble, the pressure of the system is kept constant instead of the volume. This means that in addition to a thermostat regulating the temperature, a barostat is required to do the same for the pressure.


#### Barostat

For the three mentioned thermostats, there are corresponding barostats. In the case of the Berendsen barostat, the method is much the same as for the thermostat. Instead of velocities, however, it is the size of the periodic box containing the system and the coordinates within it that are scaled. If the pressure is too low, the box size is decreased and all atoms within it are brought closer to each other, and for a pressure that is too high the opposite is done. The Andersen and Nosé--Hoover barostats adopt a similar approach to the Nosé--Hoover thermostat, creating a coupling to an artificial system. This system acts as a piston, with artificial mass determining the strength of the coupling, trying to compress the real system.


## Time integration

As described in previous chapters, the potential energy of a system can be obtained for a given set of nuclear positions either through QM methods or through more approximate means such as molecular dynamics (MD). Given this knowledge, it is now possible to make the atoms move. Following the Born--Oppenheimer approximation, the behaviour of the nuclei should be described using MD, but such calculations are too resource demanding for any larger systems. As such, this work deals only with classical MD, in which the nuclei move according to Newton's equations.


### Euler integration

The first thing that is needed is the acceleration of each nucleus, obtained from the derivative of the potential energy function with respect to the nuclear coordinates. For the force field energy, a simple analytical function of $3N$ coordinates, this is easily done, but for the QM energy things are more complicated, possibly requiring numerical differentiation. Based on these derivatives, the acceleration for atom $j$ can easily be found using Newton's second law:
%
\begin{equation*}
\label{eq:acceleration}
\mathbf{a}_j = \frac{\mathbf{F}_j}{m_j} = - \frac{1}{m_j} \mbox{\boldmath$\nabla$}_j E,
\end{equation*}
%
where $m_j$ is the mass of atom $j$ and $\mathbf{F}_j$ are the forces acting on it.
If the atomic coordinates at a given time $t = t_i$ are $\mathbf{R}(t_i)$, the positions after a small time step, $\Delta t$, can be expressed as a Taylor expansion around $t = t_i$:
%
\begin{equation*}
\mathbf{R}(t_i + \Delta t)
= \mathbf{R}(t_i) + \frac{\partial\mathbf{R}}{\partial t} \Delta t + \frac{1}{2} \frac{\partial^2\mathbf{R}}{\partial t^2} \Delta t^2 + \frac{1}{6} \frac{\partial^3\mathbf{R}}{\partial t^3} \Delta t^3 + \cdots,
\end{equation*}
%
where each partial derivative is evaluated at $t = t_i$. Given a fixed set of time steps, each of length $\Delta t$, a series of positions are obtained:
%
\begin{equation*}
\mathbf{R}_{i+1}
= \mathbf{R}_i + \frac{\partial\mathbf{R}_i}{\partial t} \Delta t + \frac{1}{2} \frac{\partial^2\mathbf{R}_i}{\partial t^2} \Delta t^2 + \frac{1}{6} \frac{\partial ^3\mathbf{R}_i}{\partial t^3} \Delta t^3 + \cdots 
= \mathbf{R}_i + \mathbf{v}_i \Delta t + \frac{1}{2} \mathbf{a}_i \Delta t^2 + \frac{1}{6} \mathbf{j}_i \Delta t^3 + \mathcal{O}(\Delta t^4),
\end{equation*}
%
where $\mathbf{v}_i$, $\mathbf{a}_i$ and $\mathbf{j}_i$ are the velocities, acceleration and jerk of time step $i$, respectively. Truncating this to the first order gives
%
\begin{equation*}
\mathbf{R}_{i+1} = \mathbf{R}_i + \mathbf{v}_i \Delta t + \mathcal{O}(\Delta t^2),
\end{equation*} 
%
which can be differentiated with respect to time to obtain the velocities for step $i+1$:
%
\begin{align*}
\frac{\partial\mathbf{R}_{i+1}}{\partial t}&= \frac{\partial\mathbf{R}_i}{\partial t} + \frac{\partial\mathbf{v}_i }{\partial t}\Delta t + \mathcal{O}(\Delta t^2),\\
\mathbf{v}_{i+1}& = \mathbf{v}_i + \mathbf{a}_i \Delta t + \mathcal{O}(\Delta t^2).
\end{align*}
These two equations, together with the expression for the acceleration found in Equation \eqref{eq:acceleration} make up the Euler method. For each step, the acceleration is obtained from Equation \eqref{eq:acceleration}, then the velocities are calculated with Equation \eqref{eq:eulerv} and finally a new set of coordinates with Equation \eqref{eq:eulerr} before the process is repeated, step by step. This is an extremely simple method, with a local truncation error in order of $\Delta t^2$ for each step and a global error of first order for the trajectory as a whole. This can be seen from the fact that a trajectory of length $T$ requires $M = T / \Delta t$ steps. For each step, an error of order $\Delta t^2$ is added, leading to a total error of $T/ \Delta t \hspace{1mm} \mathcal{O}(\Delta t^2) = \mathcal{O}(\Delta t)$, i.e.\ first order. A simple example of integration using the Euler method is shown in Figure \ref{fig:euler}.
```{figure} ../img/md/MD_euler_even.svg
width: 300px
name: fig_euler
Simulation of a harmonic oscillator using Euler integration. The integration time step, $\Delta t$ is given as a fraction of the period of the oscillator, $T$.
```

### Verlet integration

A better choice of integrator can be found in the Verlet algorithm.[cite verlet1967] In this algorithm, the coordinates of the previous time step, $\mathbf{R}_{i-1}$, are required, which can be obtained by replacing $\Delta t$ with $-\Delta t$ in Equation \eqref{eq:md1}:
%
\begin{equation*}
\label{eq:md3}
\mathbf{R}_{i-1} = \mathbf{R}_i - \mathbf{v}_i \Delta t + \frac{1}{2} \mathbf{a}_i \Delta t^2 - \frac{1}{6} \mathbf{j}_i \Delta t^3 + \mathcal{O}(\Delta t^4).
\end{equation*}

If Equations \eqref{eq:md1} and \eqref{eq:md3} are added together, the velocity and jerk terms cancel and the expression for $\mathbf{R}_{i+1}$ can be written as
%
\begin{equation*}
\label{eq:md4}
\mathbf{R}_{i+1} = (2\mathbf{R}_i - \mathbf{R}_{i-1}) + \mathbf{a}_i \Delta t^2 + \mathcal{O}(\Delta t^4).
\end{equation*}

Neglecting higher order terms and using the acceleration from Equation \eqref{eq:acceleration}, this expression makes up the Verlet algorithm, with a local truncation error in the order of $\Delta t^4$. The global error, however, is of second order,[cite jensen2006] due to the presence of two previous coordinates in the expression. Figure \ref{fig:verlet} shows the improvement over the Euler method.


```{figure} ../img/md/MD_verlet_even.svg
width: 300px
name: fig_MD_verlet_even
Simulation of a harmonic oscillator using Verlet integration. The integration time step, $\Delta t$ is given as a fraction of the period of the oscillator, $T$.
```


### Velocity Verlet

In any simulation that deals with temperature, the kinetic energy of the system becomes a factor, most commonly dealt with through the atomic velocities. As the Verlet algorithm neither calculates nor uses these velocities, measurements and alterations of the temperature become difficult. For this reason, the velocity Verlet algorithm[cite swope1982] was created. This uses the Taylor expansion of Equation \eqref{eq:md1} up to the second order for the coordinates:
%
\begin{equation*}
\mathbf{R}_{i+1} = \mathbf{R}_i + \mathbf{v}_i \Delta t + \frac{1}{2} \mathbf{a}_i \Delta t^2 + \mathcal{O}(\Delta t^3).
\end{equation*}

Differentiating this expression with respect to time gives the expression for the velocity of the next step:
%
\begin{align*}
\frac{\partial\mathbf{R}_{i+1}}{\partial t} &= \frac{\partial \mathbf{R}_i}{\partial t} + \frac{\partial\mathbf{v}_i}{\partial t} \Delta t + \frac{1}{2} \frac{\partial \mathbf{a}_i}{\partial t} \Delta t^2 + \mathcal{O}(\Delta t^3),\\
\mathbf{v}_{i+1} &= \mathbf{v}_i + \mathbf{a}_i \Delta t + \frac{1}{2} \mathbf{j}_i \Delta t^2 + \mathcal{O}(\Delta t^3).
\end{align*}

An expression for $\mathbf{j}_i$ can be found by differentiating once more and truncating after the second term:
%
\begin{align*}
\frac{\partial \mathbf{v}_{i+1}}{\partial t} &= \frac{\partial\mathbf{v}_i}{\partial t} + \frac{\partial\mathbf{a}_i}{\partial t} \Delta t + \mathcal{O}(\Delta t^2),\\
\mathbf{a}_{i+1} &= \mathbf{a}_i + \mathbf{j} \Delta t + \mathcal{O}(\Delta t^2).
\end{align*}

This can be rearranged into
%
\begin{equation*}
\mathbf{j}_i \Delta t = \mathbf{a}_{i+1} - \mathbf{a}_i + \mathcal{O}(\Delta t^2),
\end{equation*}
%
which can be inserted in Equation \eqref{eq:md6}, resulting in
%
\begin{equation*}
\label{eq:md7}
\mathbf{v}_{i+1} = \mathbf{v}_i + \frac{1}{2} (\mathbf{a}_i + \mathbf{a}_{i+1})\Delta t + \mathcal{O}(\Delta t^3).
\end{equation*}

Equations \eqref{eq:acceleration}, \eqref{eq:md5} and \eqref{eq:md7} make up the velocity Verlet algorithm. First, Equation \eqref{eq:md5} is used to calculated the coordinates of the next step, which are used in Equation \eqref{eq:acceleration} to calculate the acceleration of that step. Finally, the velocities are calculated using Equation \eqref{eq:md7}, completing the step.

As is evident from these equations, the local truncation error of the velocity Verlet algorithm is in the order of $\mathcal{O}(\Delta t^3)$, both for coordinates and velocities, but the global error is of second order, same as for the standard Verlet algorithm. This can be seen in Figure \ref{fig:vverlet}, where the total error of the velocity Verlet algorithm is comparable to that of the Verlet algorithm in Figure \ref{fig:verlet}.


```{figure} ../img/md/MD_vverlet_even.svg
width: 300px
name: fig_MD_vverlet_even
Simulation of a harmonic oscillator using velocity Verlet integration. The integration time step, $\Delta t$ is given as a fraction of the period of the oscillator, $T$.
```

---
<!-- #endregion -->