""" Task 1.4 and 1.5 """
#%% implementing two-field variational formulation
import time
st = time.time()

import numpy as np
import matplotlib.pyplot as plt

""" mesh generation """
n = 8
num_GP = 2
GA = 80*(10**9)*3.63*(10**(-4))
EI = 200*(10**9)*171*(10**(-8))
L = 5
Load = -1
l_ele = L/n
J = l_ele/2

""" global matrices """
K = np.zeros(((4*n+4),(4*n+4)))
F = np.zeros((4*n+4,1))

""" Gauss coordinates and their weights """
xi_gp, w_gp = np.polynomial.legendre.leggauss(num_GP)

""" defining basis functions and their derivatives """
def N1(xi): return 0.5 * (1.0 - xi)
def N2(xi): return 0.5 * (1.0 + xi)
def dN1_dxi(_): return -0.5
def dN2_dxi(_): return  0.5

""" loop on elements """
for i_ele in range(n):
    k_ele = np.zeros((8,8))
    for xi, wi in zip(xi_gp, w_gp):
        n1, n2 = N1(xi), N2(xi)
        dn1dx = dN1_dxi(xi) / J
        dn2dx = dN2_dxi(xi) / J
        wprime_phi = np.array ([dn1dx, n1, 0, 0, dn2dx, n2, 0, 0])
        phiprime = np.array ([0, dn1dx, 0, 0, 0, dn2dx, 0 , 0])
        M = np.array ([0, 0, n1, 0, 0, 0, n2 , 0])
        Q = np.array ([0, 0, 0, n1, 0, 0, 0 , n2])
        #print(phiprime)
        #print(wprime_phi)

        k_ele += ((1 / EI) * np.outer(M, M) + (1 / GA) * np.outer(Q, Q)
                  - np.outer(wprime_phi, Q) - np.outer(Q, wprime_phi)
                  - np.outer(phiprime, M) - np.outer(M, phiprime)) * wi * J

    # assembling the global K
    start = 4 * i_ele
    K[start:start + 8, start:start + 8] += k_ele

""" inserting boundary conditions """
F[4 * n, 0] = -Load
for dof in (0, 1):
    K[dof,:] = 0.0
    K[:,dof] = 0.0
    K[dof,dof] = 1.0
    F[dof,0] = 0.0

""" solving """
coefficients = np.linalg.solve(K, F)
w_nodes  = coefficients[0::4, 0]
phi_nodes = coefficients[1::4, 0]
M_nodes  = coefficients[2::4, 0]
V_nodes  = coefficients[3::4, 0]

x_nodes = np.linspace(0.0, L, n + 1)
plt.plot(x_nodes, w_nodes, marker="o")
plt.xlabel("x [m]")
plt.ylabel("w [m]")
plt.title("two-field Timoshenko beam with " + str(n) + " elements")
plt.grid(True, linestyle="-", alpha=0.6)
plt.tight_layout()
plt.show()

end = time.time()
print("Running time: " + str(end - st) + " s")
#%% conduct the converge study L2 norm
ref2 = 0
L2 = 0
import timo as tm

for i in range(n):
    xL = i*l_ele
    xR = xL + l_ele
    wi = w_nodes[i]
    wj = w_nodes[i+1]

    for xi, wi_gauss in zip(xi_gp, w_gp):
        # our x points
        xC = xL + (xi+1)*J

        # displacement at GP using one-field and Euler bernoulli
        w_one = N1(xi)*wi + N2(xi)*wj
        w_ref = tm.w(xC, Load, L, EI)

        L2 += ((w_one-w_ref)**2)*wi_gauss*J
        ref2 += (w_ref**2)*wi_gauss*J

err_L2 = np.sqrt(L2)
print(f"L2 norm error (displacement): {err_L2:.6e} m")

plt.plot(x_nodes, w_nodes, marker="s", color="darkblue", label="displacement (two-field)")
plt.plot(x_nodes, tm.w(x_nodes, Load, L, EI), linestyle="--", color="red", label="displacement (Euler-Bernoulli)")
plt.xlabel("x [m]")
plt.ylabel("displacement [m]")
plt.title("two-field Timoshenko beam vs Euler-Bernoulli")
plt.grid(True, linestyle="-", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()