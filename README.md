# PDE
These codes solve, with the finite difference method, the equations of heat and waves (in one and two spatial dimensions). For the heat there is a heated bar which cools to a uniform temperature; for the waves, in one dimension there is an elastic rope with fixed ends, in two dimensions it is the mambrana of a rectangular drum.
There is also an animation of the solution.

In these codes the solution method consists simply in approximating the derivatives with the incremental ratio and adding everything


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\frac{df(x,&space;y)}{dx}&space;=&space;\frac{f(x_{k&plus;1},&space;y_k)-f(x_k,&space;y_k)}{h}&space;&plus;&space;\mathcal{O}(h)\\&space;\frac{d^2f(x,&space;y)}{dx^2}&space;=&space;\frac{f(x_{k&plus;1},&space;y_k)-2f(x_k,&space;y_k)&plus;f(x_{k-1},&space;y_k)}{h^2}&space;&plus;&space;\mathcal{O}(h^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\frac{df(x,&space;y)}{dx}&space;=&space;\frac{f(x_{k&plus;1},&space;y_k)-f(x_k,&space;y_k)}{h}&space;&plus;&space;\mathcal{O}(h)\\&space;\frac{d^2f(x,&space;y)}{dx^2}&space;=&space;\frac{f(x_{k&plus;1},&space;y_k)-2f(x_k,&space;y_k)&plus;f(x_{k-1},&space;y_k)}{h^2}&space;&plus;&space;\mathcal{O}(h^2)" title="\\ \frac{df(x, y)}{dx} = \frac{f(x_{k+1}, y_k)-f(x_k, y_k)}{h} + \mathcal{O}(h)\\ \frac{d^2f(x, y)}{dx^2} = \frac{f(x_{k+1}, y_k)-2f(x_k, y_k)+f(x_{k-1}, y_k)}{h^2} + \mathcal{O}(h^2)" /></a>


![](eq_onda_2d.gif)