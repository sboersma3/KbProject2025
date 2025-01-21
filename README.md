# KB project 2025

For questions reach out to: sjoerd.boersma@wur.nl


## Model 
This directory contains a model that ca be used to generate simulation data. The model is defined in state-space format (continuous time and discrete time):

$$
\begin{aligned}
\frac{\text{d}x(t)}{\text{d}t} &= f\left( x(t),u(t),d(t),w(t),p \right), \qquad  &y(t) = g\left( x(t),u(t),d(t),v(t),p \right) \\
x(k+1) &= f_\text{d}\left( x(k),u(k),d(k),w(k),p \right), \qquad &y(k) = g\left( x(k),u(k),d(k),v(k),p \right)
\end{aligned}
$$

with

symbol|meaning
---|---
$x(t) \in \mathbb{R}^6$|state variable
$y(t) \in \mathbb{R}^{13}$|measurement
$u(t) \in \mathbb{R}^3$|control signal
$d(t) \in \mathbb{R}^4$|weather disturbance
$w(t)\in \mathbb{R}^6$|process noise 
$v(t)\in \mathbb{R}^{13}$|measurement noise
$p$|model parameter
$f(\cdot),f_\text{d}(\cdot)$|nonlinear function (model)
$t \in \mathbb{R},k \in \mathbb{Z}$|continuous, discrete time

The data that is generated with the model can be used to test hypothesis, before using real data. The meaning of the variables is defined below.



## Optimal Control


## Reinforcement Learning

