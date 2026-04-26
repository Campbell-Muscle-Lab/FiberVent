---
layout: default
title: Growth
parent: Modules
has_children: false
nav_order: 2
---

# Growth

## Overview

FiberVent's hemi-elliptical ventricle can grow in two modes: concentric (wall thickening / thinning) and eccentric (chamber dilation / constriction).

These two behaviors are defined in the growth object in the [model file](../../files/model/model_file.html).

Here is an example.

```
"growth":
{
    "master_rate": 1,
    "control":
    [
        {
            "type": "concentric",
            "level": "muscle",
            "signal": "muscle_ATP_concentration",
            "set_point": 0.007,
            "prop_gain": -0.01,
            "deriv_gain": -0.3,
            "control_period_s": 5
        },
        {
            "type": "eccentric",
            "level": "FiberSim_half_sarcomere",
            "signal": "fs_stress_titin",
            "set_point": 105,
            "prop_gain": 0.01,
            "control_period_s": 10,
            "max_rate": 0.01
        }
    ]
}
```

## Concentric growth

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

 Concentric growth is defined by the differential equation 

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
