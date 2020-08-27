# Option Pricing under Feller-Lévy models. 

* This repository contains the C++ codes for my work as a research assistant at the Seminar for Applied Mathematics at ETH Zürich. 
* The topic of my project was: Option Pricing under Feller-Lévy models. 
* The codes contain the Finite Element Method implementation for option pricing. 

## Overview
The arbitrage-free price of financial products with payoff <a href="https://www.codecogs.com/eqnedit.php?latex=g" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g" title="g" /></a> at time <a href="https://www.codecogs.com/eqnedit.php?latex=t&space;\in&space;[0,T]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t&space;\in&space;[0,T]" title="t \in [0,T]" /></a>, is given by
<a href="https://www.codecogs.com/eqnedit.php?latex=V(t,x)&space;=&space;\mathbb{E}[e^{-rT}g(X_T)&space;|&space;X_t=x]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V(t,x)&space;=&space;\mathbb{E}[e^{-rT}g(X_T)&space;|&space;X_t=x]" title="V(t,x) = \mathbb{E}[e^{-rT}g(X_T) | X_t=x]" /></a>.

The stock price process follow <a href="https://www.codecogs.com/eqnedit.php?latex=dX_t&space;=&space;b(X_{t^-})dt&space;&plus;&space;a(X_{t^-})^{1/\alpha}d&space;L_t^{\alpha},&space;\&space;\&space;\&space;X_0=x," target="_blank"><img src="https://latex.codecogs.com/gif.latex?dX_t&space;=&space;b(X_{t^-})dt&space;&plus;&space;a(X_{t^-})^{1/\alpha}d&space;L_t^{\alpha},&space;\&space;\&space;\&space;X_0=x," title="dX_t = b(X_{t^-})dt + a(X_{t^-})^{1/\alpha}d L_t^{\alpha}, \ \ \ X_0=x," /></a>, where  <a href="https://www.codecogs.com/eqnedit.php?latex=(L_t^{\alpha})_{t\ge0}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(L_t^{\alpha})_{t\ge0}" title="(L_t^{\alpha})_{t\ge0}" /></a> is an alpha-stable Lévy process, and hence the stock price is a jump process.

Using the Feynman-Kac theorem, the fractional partial differential equation governing the price of the option is given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation*}&space;\label{pricing_eq}&space;\left\{&space;\begin{array}{@{}ll@{}}&space;v_t(t,x)&space;-&space;b(x)\partial_xv(t,x)&space;&plus;&space;a(x)(-\partial_{xx})^{\alpha/2}v(t,x)&space;&plus;&space;rv(t,x)&space;=&space;0,&space;&&space;\text{in}\&space;\mathbb{R}_{>0}\times\mathbb{R}_{\ge0},&space;\\&space;v(0,x)&space;=&space;g(x),&space;&&space;\forall&space;x\in\mathbb{R}_{\ge0}.&space;\end{array}\right.&space;\end{equation*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation*}&space;\label{pricing_eq}&space;\left\{&space;\begin{array}{@{}ll@{}}&space;v_t(t,x)&space;-&space;b(x)\partial_xv(t,x)&space;&plus;&space;a(x)(-\partial_{xx})^{\alpha/2}v(t,x)&space;&plus;&space;rv(t,x)&space;=&space;0,&space;&&space;\text{in}\&space;\mathbb{R}_{>0}\times\mathbb{R}_{\ge0},&space;\\&space;v(0,x)&space;=&space;g(x),&space;&&space;\forall&space;x\in\mathbb{R}_{\ge0}.&space;\end{array}\right.&space;\end{equation*}" title="\begin{equation*} \label{pricing_eq} \left\{ \begin{array}{@{}ll@{}} v_t(t,x) - b(x)\partial_xv(t,x) + a(x)(-\partial_{xx})^{\alpha/2}v(t,x) + rv(t,x) = 0, & \text{in}\ \mathbb{R}_{>0}\times\mathbb{R}_{\ge0}, \\ v(0,x) = g(x), & \forall x\in\mathbb{R}_{\ge0}. \end{array}\right. \end{equation*}" /></a>

where r the risk free interest rate.

Finally, the finite element method is applied to the above equation to solve for the option price process. 

