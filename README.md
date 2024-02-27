# Moire_continuum_model
Based on continuum model to construct the Hamiltonian of twisted bilayer (multi-layers) system


## Introduction
Continuum model transforms the real space basis to k space basis. 


## Dirac Hamiltonian
For simple graphene, the Hamiltonian is given as: 

$$
H(\boldsymbol{q})=\hbar v_F\left[\begin{array}{cc}
0 & \pm q_x-i q_y \\
\pm q_x+i q_y & 0
\end{array}\right]= \pm \hbar v_F \vec{\sigma} \cdot \vec{q}
$$




where $v_F=\frac{3dt}{2\hbar},\vec{\sigma} $ is the Pauli matrix vector. $\vec{q}$ is the $k$ vector relative to the $K\ or\ K'$ point.

## Rotated Dirac Hamiltonian

After rotation of angle $\theta$, the hamiltonian is:

$$
H_\theta(q)=-v (q_xcos(\theta_k-\theta),q_ysin(\theta_k-\theta))\cdot \begin{bmatrix}\sigma_x\\\sigma_y\end{bmatrix}=-v q\left[\begin{array}{cc}
0 & e^{i\left(\theta_{\mathbf{k}}-\theta\right)} \\
e^{-i\left(\theta_{\mathbf{k}}-\theta\right)} & 0
\end{array}\right]
$$

The hoppings between different bases are constrained by the momentum conservation.


## TBG Hamiltonian
### Selection rule
The hoppings of intralayer atoms can be easily determined. The problem is to determine the interlayer one $t_\perp^{\alpha\beta}(\vec r)$

$$
T_{12}^{\alpha \beta}\left(\mathbf{q}_1, \mathbf{q}_2\right)=\frac{1}{A_{\text {u.c. }}} \sum_{\mathbf{G}_1, \mathbf{G}_2} e^{i \mathbf{G}_1 \cdot \boldsymbol{\tau}_{1, \alpha}} t_{12}^{\alpha \beta}\left(\mathrm{K}_1+\mathbf{q}_1+\mathbf{G}_1\right) e^{-i \mathbf{G}_2 \cdot \boldsymbol{\tau}_{2, \beta}} \delta_{\mathrm{K}_1+\mathbf{q}_1+\mathbf{G}_1, \mathrm{~K}_2+\mathbf{q}_2+\mathbf{G}_2}
$$
for $|q_1|,|1_2|\ll |K|$, $
\mathbf{G}_{\ell}=\mathbf{g}_{\ell, 1}, \mathbf{g}_{\ell, 2}, \mathbf{g}_{\ell, 3} $  with $ \mathbf{g}_{\ell, 1}=\mathbf{0}, \mathbf{g}_{\ell, 2}=\mathbf{b}_{\ell, 2} \text { and } \mathbf{g}_{\ell, 3}=-\mathbf{b}_{\ell, 1}
$

In Macdonald(2011),

$$
T^{\alpha \beta}(\mathbf{r})=w \sum_{j=1}^3 \exp \left(-i \mathbf{q}_j \cdot \mathbf{r}\right) T_j^{\alpha \beta}
$$

$$ 
\begin{gathered}
T_1=\mathrm{T}_{\mathbf{q}_{\mathrm{b}}}=\frac{t_{\perp}(|\mathrm{K}|)}{A_{u . c .}}\left[\begin{array}{ll}
1 & 1 \\
1 & 1
\end{array}\right], \\
T_2=\mathrm{T}_{\mathbf{q}_{t r}}=\frac{t_{\perp}(|\mathrm{K}|)}{A_{u . c .}} e^{-i \mathbf{g}_{1,2} \cdot \tau_0}\left[\begin{array}{cc}
e^{i \phi} & 1 \\
e^{-i \phi} & e^{i \phi}
\end{array}\right], \\
T_3=\mathrm{T}_{\mathbf{q}_{t l}}=\frac{t_{\perp}(|\mathrm{K}|)}{A_{\text {u.c. }}} e^{-i \mathbf{g}_{1,3} \cdot \tau_0}\left[\begin{array}{cc}
e^{-i \phi} & 1 \\
e^{i \phi} & e^{-i \phi}
\end{array}\right],
\end{gathered}
$$
where $\phi=\frac{2}{3}\pi$

For the k space Hamiltonian, only following items are non-zero:
$$q_2-q_1=q_b\\
q_2-q_1=q_{tr}\\
q_2-q_1=q_{tl}
$$
So,
$$
T_{21}^{\alpha \beta}\left(\mathbf{q}_2, \mathbf{q}_1\right)=\left(T_{\mathbf{q}_b}^{\beta \alpha}\right)^* \delta_{\mathbf{q}_2-\mathbf{q}_1, \mathbf{q}_b}+\left(T_{\mathbf{q}_{t r}}^{\beta \alpha}\right)^* \delta_{\mathbf{q}_2-\mathbf{q}_1, \mathbf{q}_{t r}}+\left(T_{\mathbf{q}_{t l}}^{\beta \alpha}\right)^* \delta_{\mathbf{q}_2-\mathbf{q}_1, \mathbf{q}_{t l}}
$$

### Two example TBG Hamiltonians
With Hamiltonian written in 4 basis: $
\left|1, \mathrm{~K}_1+\mathbf{q}, \alpha\right\rangle,\left|2, \mathrm{~K}_2+\mathbf{q}+\mathbf{q}_b, \alpha\right\rangle,\left|2, \mathrm{~K}_2+\mathbf{q}+\mathbf{q}_{t t^{\prime}} \alpha\right\rangle,\left|2, \mathrm{~K}_2+\mathbf{q}+\mathbf{q}_{t p^{\prime}} \alpha\right\rangle
$

Based on the selection rule, the Hamiltonian is:
$$
H_{4 }^{\mathrm{K}}(\mathbf{q})=\left[\begin{array}{cccc}
H_1^{\mathrm{K}}(\mathbf{q}) & T_{\mathbf{q}_b} & T_{\mathbf{q}_{t r}} & T_{\mathbf{q}_{t l}} \\
T_{\mathbf{q}_b}^{\dagger} & H_2^{\mathrm{K}}\left(\mathbf{q}+\mathbf{q}_b\right) & 0 & 0 \\
T_{\mathbf{q}_{t r}}^{\dagger} & 0 & H_2^{\mathrm{K}}\left(\mathbf{q}+\mathbf{q}_{t r}\right) & 0 \\
T_{\mathbf{q}_{t l}}^{\dagger} & 0 & 0 & H_2^{\mathrm{K}}\left(\mathbf{q}+\mathbf{q}_{t l}\right)
\end{array}\right]
$$


With Hamiltonian written in 10 basis:$
\begin{gathered}
|1,(0,0)\rangle,|2,(0,0)\rangle,|2,(0,1)\rangle|2,(-1,0)\rangle, 
|1,(0,-1)\rangle,|1,(1,0)\rangle,|1,(0,1)\rangle,|1,(1,1)\rangle,|1,(-1,0)\rangle,\\|1,(-1,-1)\rangle,
\end{gathered}
$

$$
H_{10, \mathrm{tLLG}}^{\mathrm{K}}(\mathbf{q})=\left[\begin{array}{cccccccccc}
H_1^{\mathrm{K}} & T_{\mathbf{q}_b} & T_{\mathbf{q}_{t r}} & T_{\mathbf{q}_{t l}} & 0 & 0 & 0 & 0 & 0 & 0 \\
T_{\mathbf{q}_b}^{\dagger} & H_2^{\mathrm{K}} & 0 & 0 & T_{\mathbf{q}_{t r}}^{\dagger} & T_{\mathbf{q}_{t l}}^{\dagger} & 0 & 0 & 0 & 0 \\
T_{\mathbf{q}_{t r}}^{\dagger} & 0 & H_2^{\mathrm{K}} & 0 & 0 & 0 & T_{\mathbf{q}_b}^{\dagger} & T_{\mathbf{q}_{t l}}^{\dagger} & 0 & 0 \\
T_{\mathbf{q}_{t l}}^{\dagger} & 0 & 0 & H_2^{\mathrm{K}} & 0 & 0 & 0 & 0 & T_{\mathbf{q}_b}^{\dagger} & T_{\mathbf{q}_{t r}}^{\dagger} \\
0 & T_{\mathbf{q}_{t r}} & 0 & 0 & H_1^{\mathrm{K}} & 0 & 0 & 0 & 0 & 0 \\
0 & T_{\mathbf{q}_{t l}} & 0 & 0 & 0 & H_1^{\mathrm{K}} & 0 & 0 & 0 & 0 \\
0 & 0 & T_{\mathbf{q}_b} & 0 & 0 & 0 & H_1^{\mathrm{K}} & 0 & 0 & 0 \\
0 & 0 & T_{\mathbf{q}_{t l}} & 0 & 0 & 0 & 0 & H_1^{\mathrm{K}} & 0 & 0 \\
0 & 0 & 0 & T_{\mathbf{q}_b} & 0 & 0 & 0 & 0 & H_1^{\mathrm{K}} & 0 \\
0 & 0 & 0 & T_{\mathbf{q}_{t r}} & 0 & 0 & 0 & 0 & 0 & H_1^{\mathrm{K}}
\end{array}\right]
$$

## Band plot
Interlayer hopping $110meV$ is empirical value.






