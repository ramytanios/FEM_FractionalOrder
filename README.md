# Option Pricing under Feller-Lévy models. 

* This repository contains the C++ codes for my work as a research assistant at the Seminar for Applied Mathematics at ETH Zürich. 
* The topic of my project was: Option Pricing under Feller-Lévy models. 
* The codes contain the Finite Element Method implementation for option pricing. 
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">

## Overview
The arbitrage-free price of financial products with payoff <img src="https://render.githubusercontent.com/render/math?math=g"> at time <img src="https://render.githubusercontent.com/render/math?math=t \in [0,T]">, is given by
<img src="https://render.githubusercontent.com/render/math?math=V(t,x) = \mathbb{E}[\e^{-rT}g(X_T) | X_t=x]">.

The price process follow <img src="https://render.githubusercontent.com/render/math?math=dX_t = b(X_{t^-})dt + a(X_{t^-})^{1/\alpha}d L_t^{\alpha}, \ \ \ X_0=x, ">, where  <img src="https://render.githubusercontent.com/render/math?math=(L_t^{\alpha})_{t\ge0}"> is an alpha-stable Lévy process.

Using the Feynman-Kac theorem, the fractional partial differential equation governing the price of the option is given by <img src="https://render.githubusercontent.com/render/math?math=v_t(t,x) - b(x)\partial_xv(t,x) + a(x)(-\partial_{xx})^{\alpha/2}v(t,x) + rv(t,x) = 0, & \text{in}\ \mathbb{R}_{>0}\times\mathbb{R}_{\ge0}">, with suitable initial conditions and r the risk free interest rate.

Finally, the finite element method is applied to the above equation to solve for the option price process. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{E}&space;&plus;&space;r" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{E}&space;&plus;&space;r" title="\mathbb{E} + r" /></a>
