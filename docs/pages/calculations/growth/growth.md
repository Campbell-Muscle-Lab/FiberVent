---
layout: default
title: Growth
parent: Calculations
has_children: false
nav_order: 2
---

# Growth

## Overview

FiberVent's hemi-elliptical ventricle can grow in two modes: concentric (wall thickening / thinning) and eccentric (chamber dilation / constriction).


## Concentric growth

Concentric growth mimics myocytes add / removing sarcomeres and mitochondria in parallel by changeing the wall thickness $T$. The relevant differential equation is

$$
\begin{equation}
\frac{dT(t)}{dt} = G_m T(t) g_c(t)
\end{equation}
$$

where:
+ $G_m$ is the master growth rate
+ $g_c$ is the concentric growth signal

In turn, $g_c$ is defined by:

$$
\begin{equation}
\frac{dg_c(t)}{dt} = g_{cp} s_c(t) + g_{cd}\frac{ds_c(t)}{dt}
\end{equation}
$$


where
+ $g_{cp}$ is scaling factor defining proportional feedback
+ $g_{cd}$ is scaling factor defining derivative feedback'

and

$$
\begin{equation}
s_c(t) = \frac{[ATP] - ATP_{set}}{ATP_{set}}
\end{equation}
$$

where:

+ $ATP_{set}$ is a homeostatic set-point


## Eccentric growth

Eccentric growth mimics myocytes adding / removing sarcomeres in series by changing the number of half-sarcomeres around the ventricle's circumference $n$. The relevant differential equation is

$$
\begin{equation}
\frac{dn(t)}{dt} = G_m n(t) g_e(t)
\end{equation}
$$

where:
+ $G_m$ is the master growth rate
+ $g_e$ is the ecceentric growth signal

In turn, $g_e$ is defined by:

$$
\begin{equation}
\frac{dg_e(t)}{dt} = g_{ep} s_e(t)
\end{equation}
$$

and

$$
\begin{equation}
s_e(t) = \frac{F_{titin}(t) - F_{titin\_set}}{F_{titin\_set}}
\end{equation}
$$

where:

+ $F_{titin}$ is the stress in titin molecules
+ $F_{titin_set}$ is a homestatic set-point
