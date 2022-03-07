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

An example is in the code trasp_lw.f

If desired, the method can be generalized to an equation such as:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u(x,t)}{\partial&space;t}&plus;\frac{\partial&space;f(u(x,t))}{\partial&space;x}=0" title="http://latex.codecogs.com/svg.latex?\frac{\partial u(x,t)}{\partial t}+\frac{\partial f(u(x,t))}{\partial x}=0" />

For example in the code traspnl_lw.f is solved:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;\frac{\partial&space;f(u)}{\partial&space;x}&space;\hspace{5&space;mm}&space;\text{with}&space;\hspace{2.5&space;mm}&space;f(u)=&space;v&space;u^2&space;\Rightarrow&space;\frac{\partial&space;u}{\partial&space;t}&space;=&space;2v\frac{\partial&space;u}{\partial&space;x}" title="http://latex.codecogs.com/svg.latex?\frac{\partial u}{\partial t} = \frac{\partial f(u)}{\partial x} \hspace{5 mm} \text{with} \hspace{2.5 mm} f(u)= v u^2 \Rightarrow \frac{\partial u}{\partial t} = 2v\frac{\partial u}{\partial x}" />

the plot.py code can be used to analyze the results of fortran codes.

Another interesting equation to deal with is the heat equation, which is a diffusion equation:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;D&space;\frac{\partial^2&space;u}{\partial&space;x^2}" title="http://latex.codecogs.com/svg.latex?\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}" />

this time we can use the Forward Time Centered Space that produce:

<img src="http://latex.codecogs.com/svg.latex?u_{i}^{n&space;&plus;&space;1}&space;=&space;&space;u_{i}^{n}&space;&plus;&space;\frac{D&space;\Delta&space;t}{\Delta&space;x^2}(u_{i&space;&plus;&space;1}^{n}&space;-&space;2&space;u_{i}^{n}&space;&plus;&space;u_{i&space;-&space;1}^{n})" title="http://latex.codecogs.com/svg.latex?u_{i}^{n + 1} = u_{i}^{n} + \frac{D \Delta t}{\Delta x^2}(u_{i + 1}^{n} - 2 u_{i}^{n} + u_{i - 1}^{n})" />

the stability analysis shows that the scheme is stable if:  <img src="http://latex.codecogs.com/svg.latex?\frac{D&space;\Delta&space;t}{\Delta&space;x^2}&space;<&space;\frac{1}{2}" title="http://latex.codecogs.com/svg.latex?\frac{D \Delta t}{\Delta x^2} < \frac{1}{2}" />

we can also adopt an implicit scheme:

<img src="http://latex.codecogs.com/svg.latex?\\u_{i}^{n&space;&plus;&space;1}&space;=&space;u_{i}^{n}&space;&plus;&space;\frac{&space;D&space;\Delta&space;t}{\Delta&space;x^2}&space;(u_{i&space;&plus;&space;1}^{n&plus;1}&space;-&space;2&space;u_{i}^{n&plus;1}&space;&plus;&space;u_{i&space;-&space;1}^{n&plus;1})\hspace{5&space;mm}&space;r&space;=&space;&space;\frac{&space;D&space;\Delta&space;t}{\Delta&space;x^2}&space;\\\\-ru^{n&plus;1}_{i-1}&space;&plus;&space;(1&plus;2r)u^{n&plus;1}_i&space;-&space;ru^{n&plus;1}_{i&plus;1}&space;=&space;u^n_j&space;\\\\\text{we&space;must&space;therefore&space;solve:}\\\\\begin{pmatrix}1&space;&plus;2r&space;&&space;-r&space;&&space;0&space;\\-r&space;&&space;&space;1&space;&plus;2r&space;&&space;\ddots&space;\\0&space;&&space;\ddots&space;&&space;\ddots\end{pmatrix}u^{n&plus;1}&space;=&space;u^n&space;" title="http://latex.codecogs.com/svg.latex?\\u_{i}^{n + 1} = u_{i}^{n} + \frac{ D \Delta t}{\Delta x^2} (u_{i + 1}^{n+1} - 2 u_{i}^{n+1} + u_{i - 1}^{n+1})\hspace{5 mm} r = \frac{ D \Delta t}{\Delta x^2} \\\\-ru^{n+1}_{i-1} + (1+2r)u^{n+1}_i - ru^{n+1}_{i+1} = u^n_j \\\\\text{we must therefore solve:}\\\\\begin{pmatrix}1 +2r & -r & 0 \\-r & 1 +2r & \ddots \\0 & \ddots & \ddots\end{pmatrix}u^{n+1} = u^n " />

Examples of explicit and implicit schema are found respectively in: calore1D_exp.py and calore1D_imp.py.
Also in calore1D_imp.f is implemented the implicit scheme and to show solution can be used calore.py