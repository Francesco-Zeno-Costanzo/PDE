# PDE
These codes show various methods of solving different PDEs; Let's start with the transport equation:

## Transport equation ##

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

It is therefore possible to see the formation and breakage of the wave front, (and of the code when the front is vertical).The plot.py code can be used to analyze the results of fortran codes.

## Heat equation ##

Another interesting equation to deal with is the heat equation, which is a diffusion equation:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;D&space;\frac{\partial^2&space;u}{\partial&space;x^2}" title="http://latex.codecogs.com/svg.latex?\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}" />

this time we can use the Forward Time Centered Space that produce:

<img src="http://latex.codecogs.com/svg.latex?u_{i}^{n&space;&plus;&space;1}&space;=&space;&space;u_{i}^{n}&space;&plus;&space;\frac{D&space;\Delta&space;t}{\Delta&space;x^2}(u_{i&space;&plus;&space;1}^{n}&space;-&space;2&space;u_{i}^{n}&space;&plus;&space;u_{i&space;-&space;1}^{n})" title="http://latex.codecogs.com/svg.latex?u_{i}^{n + 1} = u_{i}^{n} + \frac{D \Delta t}{\Delta x^2}(u_{i + 1}^{n} - 2 u_{i}^{n} + u_{i - 1}^{n})" />

the stability analysis shows that the scheme is stable if:  <img src="http://latex.codecogs.com/svg.latex?\frac{D&space;\Delta&space;t}{\Delta&space;x^2}&space;<&space;\frac{1}{2}" title="http://latex.codecogs.com/svg.latex?\frac{D \Delta t}{\Delta x^2} < \frac{1}{2}" />

we can also adopt an implicit scheme:

<img src="http://latex.codecogs.com/svg.latex?\\u_{i}^{n&space;&plus;&space;1}&space;=&space;u_{i}^{n}&space;&plus;&space;\frac{&space;D&space;\Delta&space;t}{\Delta&space;x^2}&space;(u_{i&space;&plus;&space;1}^{n&plus;1}&space;-&space;2&space;u_{i}^{n&plus;1}&space;&plus;&space;u_{i&space;-&space;1}^{n&plus;1})\hspace{5&space;mm}&space;r&space;=&space;&space;\frac{&space;D&space;\Delta&space;t}{\Delta&space;x^2}&space;\\\\-ru^{n&plus;1}_{i-1}&space;&plus;&space;(1&plus;2r)u^{n&plus;1}_i&space;-&space;ru^{n&plus;1}_{i&plus;1}&space;=&space;u^n_j&space;\\\\\text{we&space;must&space;therefore&space;solve:}\\\\\begin{pmatrix}1&space;&plus;2r&space;&&space;-r&space;&&space;0&space;\\-r&space;&&space;&space;1&space;&plus;2r&space;&&space;\ddots&space;\\0&space;&&space;\ddots&space;&&space;\ddots\end{pmatrix}u^{n&plus;1}&space;=&space;u^n&space;" title="http://latex.codecogs.com/svg.latex?\\u_{i}^{n + 1} = u_{i}^{n} + \frac{ D \Delta t}{\Delta x^2} (u_{i + 1}^{n+1} - 2 u_{i}^{n+1} + u_{i - 1}^{n+1})\hspace{5 mm} r = \frac{ D \Delta t}{\Delta x^2} \\\\-ru^{n+1}_{i-1} + (1+2r)u^{n+1}_i - ru^{n+1}_{i+1} = u^n_j \\\\\text{we must therefore solve:}\\\\\begin{pmatrix}1 +2r & -r & 0 \\-r & 1 +2r & \ddots \\0 & \ddots & \ddots\end{pmatrix}u^{n+1} = u^n " />

Examples of explicit and implicit schema are found respectively in: calore1D_exp.py and calore1D_imp.py.
Also in calore1D_imp.f is implemented the implicit scheme and to show solution can be used calore.py.

It is possible use both implicit and explicit method via the crank_nicolson method:

<img src="https://latex.codecogs.com/svg.image?\\\frac{u_i^{n&plus;1}&space;-&space;u_i^n}{\Delta&space;t}&space;=&space;\frac{D}{2&space;(\Delta&space;x)^2}(&space;(u_{i&plus;1}^{n&plus;1}&space;-&space;2&space;u_i^{n&plus;1}&space;&plus;&space;u_{i-1}^{n&plus;1})&space;&plus;&space;&space;(u_{i&plus;1}^n&space;-&space;2&space;u_i^n&space;&plus;&space;u_{i-1}^n))\\\\\text{define:&space;}r&space;=&space;\frac{a&space;\Delta&space;t}{2(\Delta&space;x)^2}\\\\-r&space;u_{i&plus;1}^{n&plus;1}&space;&plus;&space;(1&space;&plus;&space;2r)&space;u_i^{n&plus;1}&space;-&space;r&space;u_{i-1}^{n&plus;1}&space;=&space;r&space;u_{i&plus;1}^n&space;&plus;&space;(1&space;-&space;2r)&space;u_i^n&space;&plus;&space;r&space;u_{i-1}^n" title="https://latex.codecogs.com/svg.image?\\\frac{u_i^{n+1} - u_i^n}{\Delta t} = \frac{D}{2 (\Delta x)^2}( (u_{i+1}^{n+1} - 2 u_i^{n+1} + u_{i-1}^{n+1}) + (u_{i+1}^n - 2 u_i^n + u_{i-1}^n))\\\\\text{define: }r = \frac{a \Delta t}{2(\Delta x)^2}\\\\-r u_{i+1}^{n+1} + (1 + 2r) u_i^{n+1} - r u_{i-1}^{n+1} = r u_{i+1}^n + (1 - 2r) u_i^n + r u_{i-1}^n" />

## Burger equation ##
Wanting to combine transport and diffusion, we obtain the burger equation:

<img src="http://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;&plus;&space;u&space;\frac{\partial&space;u}{\partial&space;x}&space;=&space;\nu&space;\frac{\partial^2&space;u&space;}{\partial&space;x^2}" title="http://latex.codecogs.com/svg.latex?\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u }{\partial x^2}" />

In the code burger1D_FTCS this equation is solved with the scheme:

<img src="http://latex.codecogs.com/svg.latex?u^{n&plus;1}_j&space;=&space;u^n_j&space;-&space;\frac{dt}{2dx}u^n_j(u^n_{j&plus;1}&space;-&space;u^n_{j-1})&space;&plus;&space;\frac{\nu&space;dt}{dx^2}(u^n_{j&plus;1}&space;-2u^n_j&space;&plus;&space;u^n_{j-1})" title="http://latex.codecogs.com/svg.latex?u^{n+1}_j = u^n_j - \frac{dt}{2dx}u^n_j(u^n_{j+1} - u^n_{j-1}) + \frac{\nu dt}{dx^2}(u^n_{j+1} -2u^n_j + u^n_{j-1})" />

Another possible thing is to pass a Fourier transform in space in order to make the pde become an ode: then we calculate the spatial dierivates in transform make the inverse and then we evolve over time.

<img src="https://latex.codecogs.com/svg.image?u&space;\xrightarrow{FFT}&space;\hat{u};&space;\hspace{5&space;mm}&space;2\pi&space;k&space;\hat{u}&space;\xrightarrow{IFFT}&space;\frac{\partial&space;u}{\partial&space;x}&space;" title="https://latex.codecogs.com/svg.image?u \xrightarrow{FFT} \hat{u}; \hspace{5 mm} 2\pi k \hat{u} \xrightarrow{IFFT} \frac{\partial u}{\partial x} " />



## Wave equation ##

"
On the other hand, even if we cannot see beauty in particular measured results, we can already claim to see a certain beauty in the equations which describe general physical laws. For example, in the wave equation, there's something nice about the regularity of the appearance of the x, y, z, and t. And this nice symmetry in appearance of the x, y, z, and t suggests to the mind still a  greater beauty which has to do with the four d1mensions, the possibility that space has four-dimensional symmetry, the possibility of analyzing that and the developments of the special theory of relativity. So there is plenty of intellectual beauty associated with the equations.
"

Feynman R. Feynman's Lectures On Physics Volume 2, chapter 20.3.

Let's start with the one-dimensional equation:

<img src="https://latex.codecogs.com/svg.image?\frac{\partial^2&space;u}{\partial&space;x^2}&space;-&space;v^2&space;\frac{\partial^2&space;u}{\partial&space;t^2}&space;=&space;0" title="https://latex.codecogs.com/svg.image?\frac{\partial^2 u}{\partial x^2} - v^2 \frac{\partial^2 u}{\partial t^2} = 0" />

we can use the FTCS in the first instance as shown in the code onde1D_FTCS.py

<img src="https://latex.codecogs.com/svg.image?\frac{u^{n&plus;1}_j&space;-&space;2u^n_j&space;&plus;&space;u^{n-1}_j}{\Delta&space;t^2}&space;-&space;v^2\frac{u^n_{j&plus;1}&space;-&space;2u^n_j&space;&plus;&space;u^n_{j-1}}{\Delta&space;x^2}=&space;0" title="https://latex.codecogs.com/svg.image?\frac{u^{n+1}_j - 2u^n_j + u^{n-1}_j}{\Delta t^2} - v^2\frac{u^n_{j+1} - 2u^n_j + u^n_{j-1}}{\Delta x^2}= 0" />

that requies: <img src="https://latex.codecogs.com/svg.image?\frac{v^2&space;\Delta&space;t^2}{\Delta&space;x^2}&space;<&space;1" title="https://latex.codecogs.com/svg.image?\frac{v^2 \Delta t^2}{\Delta x^2} < 1" />

But with some manipulation is possible to write:

<img src="https://latex.codecogs.com/svg.image?\\r&space;=&space;v&space;\frac{\partial&space;u}{\partial&space;x}&space;\hspace{7&space;mm}&space;s&space;=&space;\frac{\partial&space;u}{\partial&space;t}\\\\\frac{\partial^2&space;u}{\partial&space;x^2}&space;=&space;v^2&space;\frac{\partial^2&space;u}{\partial&space;t^2}&space;\Rightarrow&space;\begin{cases}\frac{\partial&space;r}{\partial&space;t}&space;=&space;v&space;\frac{\partial&space;s}{\partial&space;x}\\\frac{\partial&space;s}{\partial&space;t}&space;=&space;v&space;\frac{\partial&space;r}{\partial&space;x}\end{cases}&space;" title="https://latex.codecogs.com/svg.image?\\r = v \frac{\partial u}{\partial x} \hspace{7 mm} s = \frac{\partial u}{\partial t}\\\\\frac{\partial^2 u}{\partial x^2} = v^2 \frac{\partial^2 u}{\partial t^2} \Rightarrow \begin{cases}\frac{\partial r}{\partial t} = v \frac{\partial s}{\partial x}\\\frac{\partial s}{\partial t} = v \frac{\partial r}{\partial x}\end{cases} " />

which are two coupled transport equations that we can solve with lax wendroff as seen above. With some reworking:

<img src="https://latex.codecogs.com/svg.image?\\r^{n&plus;1}_j&space;=&space;r^n_j&space;&plus;&space;\alpha&space;\Biggl[&space;\frac{1}{2}(s^n_{j&plus;1}-s^n_{j-1})&space;&plus;&space;\frac{\alpha}{2}(r^n_{j&plus;1}&space;-&space;2r^n_j&space;+&space;r^n_{j-1})\Biggr]\\\\s^{n&plus;1}_j&space;=&space;s^n_j&space;&plus;&space;\alpha&space;\Biggl[&space;\frac{1}{2}(r^n_{j&plus;1}-r^n_{j-1})&space;&plus;&space;\frac{\alpha}{2}(s^n_{j&plus;1}&space;-&space;2s^n_j&space;+&space;s^n_{j-1})\Biggr]\\\\u^{n&plus;1}_j&space;=&space;u^n_j&space;&plus;&space;\frac{\Delta&space;t}{2}(s^n_{j&plus;1}&plus;s^n_j)" title="https://latex.codecogs.com/svg.image?\\r^{n+1}_j = r^n_j + \alpha \Biggl[ \frac{1}{2}(s^n_{j+1}-s^n_{j-1}) + \frac{\alpha}{2}(r^n_{j+1} - 2r^n_j - r^n_{j-1})\Biggr]\\\\s^{n+1}_j = s^n_j + \alpha \Biggl[ \frac{1}{2}(r^n_{j+1}-r^n_{j-1}) + \frac{\alpha}{2}(s^n_{j+1} - 2s^n_j - s^n_{j-1})\Biggr]\\\\u^{n+1}_j = u^n_j + \frac{\Delta t}{2}(s^n_{j+1}+s^n_j)" />

with: <img src="https://latex.codecogs.com/svg.image?\alpha&space;=&space;\frac{v&space;\Delta&space;t}{\Delta&space;x}" title="https://latex.codecogs.com/svg.image?\alpha = \frac{v \Delta t}{\Delta x}" />

## Laplace equation ##

From electrostatics to gravity we frequently see this equation:(would be the poisson equation, the laplace equation is the associated homogeneous)

<img src="https://latex.codecogs.com/svg.image?\nabla^2&space;\phi&space;=&space;f" title="https://latex.codecogs.com/svg.image?\nabla^2 \phi = f" />

in two dimensions we can discretize with the usual rule and obtain:

<img src="https://latex.codecogs.com/svg.image?\frac{\phi_{j&plus;1,i}&space;-&space;2\phi_{j,i}&space;&plus;&space;\phi_{j-1,i}}{\Delta&space;x^2}&space;&plus;&space;\frac{\phi_{j,i&plus;1}&space;-&space;2\phi_{j,i}&space;&plus;&space;\phi_{j,i-1}}{\Delta&space;y^2}&space;=&space;f" title="https://latex.codecogs.com/svg.image?\frac{\phi_{j+1,i} - 2\phi_{j,i} + \phi_{j-1,i}}{\Delta x^2} + \frac{\phi_{j,i+1} - 2\phi_{j,i} + \phi_{j,i-1}}{\Delta y^2} = f" />

if we choose the spatial spacings on equal x and y we can rearrange the problem in a linear system like:

<img src="https://latex.codecogs.com/svg.image?\begin{pmatrix}A&space;&&space;I&space;&&space;0&space;\\I&space;&&space;A&space;&&space;\ddots&space;\\0&space;&&space;\ddots&space;&&space;\ddots\end{pmatrix}\phi&space;=&space;f\hspace{10&space;mm}\text{with:}\hspace{4&space;mm}&space;A=\begin{pmatrix}&space;-4&space;&&space;1&space;&&space;0&space;\\1&space;&&space;-4&space;&&space;\ddots&space;\\0&space;&&space;\ddots&space;&&space;\ddots\end{pmatrix};&space;\hspace{10&space;mm}I=\begin{pmatrix}&space;1&space;&&space;0&space;&&space;0&space;\\0&space;&&space;1&space;&&space;0&space;\\0&space;&&space;0&space;&&space;\ddots\end{pmatrix}&space;&space;" title="https://latex.codecogs.com/svg.image?\begin{pmatrix}A & I & 0 \\I & A & \ddots \\0 & \ddots & \ddots\end{pmatrix}\phi = f\hspace{10 mm}\text{with:}\hspace{4 mm} A=\begin{pmatrix} -4 & 1 & 0 \\1 & -4 & \ddots \\0 & \ddots & \ddots\end{pmatrix}; \hspace{10 mm}I=\begin{pmatrix} 1 & 0 & 0 \\0 & 1 & 0 \\0 & 0 & \ddots\end{pmatrix} " />

i.e. a block triadiagonal matrix where the blocks on the diagonal are tridiagonal matrices NxN and those on the other diagonals are the identity NxN.
in the laplace2D.py code the resolution occurs through the numpy library, while in the laplace2D.f code the succesive over relaxation algorithm is implemented to solve the system. As always there is also the python code to display the result of the fortran program, in this case: laplaceplot.py
