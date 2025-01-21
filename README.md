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

$x(t)$|$y(t)$|$u(t)$|$d(t)$
---|---|---|---
enthalpy content of ice buffer | enthalpy content of ice buffer
brine temperature | brine temperature
costs of AC electricity consumption | costs of AC electricity consumption
penalty for AC power uptake above a user-defined set limit Pel_AC_max | penalty for AC power uptake above a user-defined set limit Pel_AC_max
penalty on Qib_discharge > 0.1*Pth_comp | penalty on Qib_discharge > 0.1*Pth_comp
penalty op Tbrine > zuigdrukregelwaarde + 4 | penalty op Tbrine > zuigdrukregelwaarde + 4
- | Ticebank, temperatur of water/ice in ice bank temperature 
- | h_icebank, height of liquid level in ice buffer
- | m_ice/ModelParameters.ib.pm(1), ice fraction
- | Pth_comp, thermal compressor capacity
- | Q_ib_charge, heat flow rate from water/ice to  brine in ice bank
- | Q_ib_discharge, heat flow rate from liquid refrigerant to  water/ice in ice bank
- | Pel_AC, electric power uptake from AC power grid
 
 
## Optimal Control


## Reinforcement Learning

