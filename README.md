# PDE
These codes show various methods of solving different PDEs; Let's start with the transport equation:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;&plus;&space;v&space;\frac{\partial&space;u}{\partial&space;x}&space;=&space;0" title="http://latex.codecogs.com/svg.latex?\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = 0" />

if we used:

<img src="http://latex.codecogs.com/svg.latex?\\\frac{\partial&space;u}{\partial&space;t}&space;=&space;\frac{u^{n&plus;1}_j&space;-u^n_j}{\Delta&space;t}\\\\\frac{\partial&space;u}{\partial&space;x}&space;=&space;\frac{u^n_{j&plus;1}&space;-u^n_{j-1}}{2&space;\Delta&space;x}\\&space;" title="http://latex.codecogs.com/svg.latex?\\\frac{\partial u}{\partial t} = \frac{u^{n+1}_j -u^n_j}{\Delta t}\\\\\frac{\partial u}{\partial x} = \frac{u^n_{j+1} -u^n_{j-1}}{2 \Delta x}\\ " />

we can see that the apmlitude of solution diverges for any choice of \Deltat and \Deltax; we can therefore use the Lax method where the value of the solution at time n and position j is replaced with the average of the values ​​of the solution at position j + 1 and j-1:

<img src="http://latex.codecogs.com/svg.latex?u^{n&plus;1}_j&space;=&space;\frac{u^n_{j&plus;1}&plus;u^n_{j-1}}{2}&space;-&space;\frac{v&space;\Delta&space;t}{2&space;\Delta&space;x}&space;(u^n_{j&plus;1}&space;-u^n_{j-1})&space;" title="http://latex.codecogs.com/svg.latex?u^{n+1}_j = \frac{u^n_{j+1}+u^n_{j-1}}{2} - \frac{v \Delta t}{2 \Delta x} (u^n_{j+1} -u^n_{j-1}) " />

which has as a condition of stability: <img src="http://latex.codecogs.com/svg.latex?&space;\frac{v&space;\Delta&space;t}{\Delta&space;x}<&space;1&space;" title="http://latex.codecogs.com/svg.latex? \frac{v \Delta t}{\Delta x}< 1 " />

An example is in the code trasp_lax.f

The Lax method to prevent divergences introduces a numerical diffusion in fact it is the Forward Time Centered Space of a diffusion equation. To reduce the diffusion we can use the Lax – Wendroff method which is of the second order both in time and in space:

<img src="http://latex.codecogs.com/svg.latex?\\\text{the&space;lax&space;step:}&space;\\\\u^{n&plus;1/2}_{j&plus;1/2}&space;=&space;\frac{1}{2}(u^n_{j&plus;1}&space;&plus;&space;u^n_j)&space;-&space;\frac{v&space;\Delta&space;t}{2\Delta&space;x}(u^n_{j&plus;1}&space;-&space;u^n_j)\\\\u^{n&plus;1/2}_{j-1/2}&space;=&space;\frac{1}{2}(u^n_j&space;&plus;&space;u^n_{j-1})&space;-&space;\frac{v&space;\Delta&space;t}{2\Delta&space;x}(u^n_j&space;-&space;u^n_{j-1})\\\\\text{the&space;leap&space;frog&space;step:}&space;\\\\u^{n&plus;1}_j&space;=&space;u^n_j&space;-&space;&space;&space;\frac{v&space;\Delta&space;t}{\Delta&space;x}&space;(u^{n&plus;1/2}_{j&plus;1/2}&space;-&space;u^{n&plus;1/2}_{j-1/2})" title="http://latex.codecogs.com/svg.latex?\\\text{the lax step:} \\\\u^{n+1/2}_{j+1/2} = \frac{1}{2}(u^n_{j+1} + u^n_j) - \frac{v \Delta t}{2\Delta x}(u^n_{j+1} - u^n_j)\\\\u^{n+1/2}_{j-1/2} = \frac{1}{2}(u^n_j + u^n_{j-1}) - \frac{v \Delta t}{2\Delta x}(u^n_j - u^n_{j-1})\\\\\text{the leap frog step:} \\\\u^{n+1}_j = u^n_j - \frac{v \Delta t}{\Delta x} (u^{n+1/2}_{j+1/2} - u^{n+1/2}_{j-1/2})" />

if desired, the method can be generalized to an equation such as:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u(x,t)}{\partial&space;t}&plus;\frac{\partial&space;f(u(x,t))}{\partial&space;x}=0" title="http://latex.codecogs.com/svg.latex?\frac{\partial u(x,t)}{\partial t}+\frac{\partial f(u(x,t))}{\partial x}=0" />

An example is in the code trasp_lw.f