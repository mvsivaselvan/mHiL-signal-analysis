# Time varying parameter estimation

This directory contains code to solve the 
time-varying parameter estimation/optimal control problem,

$$
\begin{gathered}
\min_{x(\cdot),u(\cdot),p(\cdot)} \int_0^T  (C x(t) - y_m(t))^\top Q (C x(t) - y_m(t))
                                          + (u(t) - u_m(t))^\top R (u(t) - u_m(t))
										  + (p(t) - p_\text{nom})^\top \rho (p(t) - p_\text{nom}) dt \\
\text{subject to:}\quad \dot{x}(t) = A(p(t))x(t) + B(p(t)) u(t)
\end{gathered}
$$