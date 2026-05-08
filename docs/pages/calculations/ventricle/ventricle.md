---
layout: default
title: Ventricle
parent: Calculations
has_children: false
nav_order: 20
---

# Ventricle

## Dimensions

FiberVent assumes the left-ventricle is a hemi-ellipse with internal radius $r$ and height $z$.

The chamber volume is given by

$$
\begin{equation}
V = \frac{2}{3} \pi r^2 z
\end{equation}
$$

where

$$
\begin{equation}
z = z_{scale} (\frac{l}{l_0})^{z_{exp}}
\end{equation}
$$

and

+ $z_{scale}$ is a scaling factor
+ $l$ and $l_0$ are the current and reference half-sarcomere lengths respectively
+ $z_{exp}$ is an exponent

If $z_{scale}$ is greater than 1, the ventricle is longer than it is wide.

The $(\frac{l}{l_0})^{z_{exp}}$ term allows the ventricle to shorten as it contracts mimicking some of the geometrical advantages provided by the helical array of sheet angles found in normal hearts.

## Pressure

FiberVent calculates the intraventricular pressure $P$ using the thick-wall approximation of Laplace's law.

$$
\begin{equation}
P = \frac {\sigma T \left(2 + \frac{T}{r} \right) } { r }
\end{equation}
$$

 where:
 + $\sigma$ is wall stress
 + $T$ is the wall thickness
 + $r$ is the internal radius.


