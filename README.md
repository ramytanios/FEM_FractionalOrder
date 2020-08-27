# Option Pricing under Feller-Lévy models. 

* This repository contains the C++ codes for my work as a research assistant at the Seminar for Applied Mathematics at ETH Zürich. 
* The topic of my project was: Option Pricing under Feller-Lévy models. 
* The codes contain the Finite Element Method implementation for option pricing. 
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">

## Overview
The arbitrage-free price of financial products with payoff <img src="https://render.githubusercontent.com/render/math?math=g"> at time <img src="https://render.githubusercontent.com/render/math?math=t \in [0,T]">, is given by
<img src="https://render.githubusercontent.com/render/math?math=V(t,x) = \mathbb{E}[\e^{-rT}g(X_T) | X_t=x]">.

The price process follow <img src="https://render.githubusercontent.com/render/math?math=dX_t = b(X_{t^-})d t + a(X_{t^-})^{1/\alpha}d L_t^{\alpha}, \ \ \ X_0=x, ">, where  <img src="https://render.githubusercontent.com/render/math?math=(L_t^{\alpha})_{t\ge0}"> is an alpha-stable Lévy process.

