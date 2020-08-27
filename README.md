<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/katex.min.css" integrity="sha384-yFRtMMDnQtDRO8rLpMIKrtPCD5jdktao2TV19YiZYWMDkUR5GQZR/NOVTdquEx1j" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/katex.min.js" integrity="sha384-9Nhn55MVVN0/4OFx7EE5kpFBPsEMZxKTCnA+4fqDmg12eCTqGi6+BB2LjY8brQxJ" crossorigin="anonymous"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.2/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous" onload="renderMathInElement(document.body);"></script>

# Option Pricing under Feller-Lévy models. 

* This repository contains the C++ codes for my work as a research assistant at the Seminar for Applied Mathematics at ETH Zürich. 
* The topic of my project was: Option Pricing under Feller-Lévy models. 
* The codes contain the Finite Element Method implementation for option pricing. 
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">

## Overview
The arbitrage-free price of financial products with payoff <img src="https://render.githubusercontent.com/render/math?math=g"> at time <img src="https://render.githubusercontent.com/render/math?math=t \in [0,T]">, is given by
<img src="https://render.githubusercontent.com/render/math?math=V(t,x) = \mathbb{E}[\myexp^{-rT}g(X_T) | X_t=x]">

