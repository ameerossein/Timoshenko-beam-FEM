""" task 1.2 and 1.3 """
#%% implementing one-field variational formulation
import time
st = time.time()

import numpy as np
import matplotlib.pyplot as plt

""" mesh generation """
n = 64
num_GP = 2
GA = 80*(10**9)*3.63*(10**(-4))
EI = 200*(10**9)*171*(10**(-8))
L = 5
Load = -1
l_ele = L/n
J = l_ele/2

""" global matrices """
K = np.zeros(((2*n+2),(2*n+2)))
F = np.zeros((2*n+2,1))

""" Gauss coordinates and their weights """
xi_gp, w_gp = np.polynomial.legendre.leggauss(num_GP)

""" defining basis functions and their derivatives """
def N1(xi): return 0.5 * (1.0 - xi)
def N2(xi): return 0.5 * (1.0 + xi)
def dN1_dxi(_): return -0.5
def dN2_dxi(_): return  0.5

""" loop on elements """
for i_ele in range(n):
    k_ele = np.zeros((4,4))
    for xi, wi in zip(xi_gp, w_gp):
        n1, n2 = N1(xi), N2(xi)
        dn1dx = dN1_dxi(xi) / J
        dn2dx = dN2_dxi(xi) / J
        wprime_phi = np.array ([dn1dx, n1, dn2dx, n2])
        phiprime = np.array ([0, dn1dx, 0, dn2dx])
        #print(phiprime)
        #print(wprime_phi)

        # inserting delta (w´+phi)GA(w´+phi)
        # inserting delta (phi´)EI(phi´)
        k_ele += (GA * np.outer(wprime_phi, wprime_phi) + EI * np.outer(phiprime, phiprime)) * wi * J

    # assembling the global K
    start = 2 * i_ele
    K[start:start + 4, start:start + 4] += k_ele

""" inserting boundary conditions """
F[2 * n, 0] = Load
for dof in (0,1):
    K[dof, :] = 0
    K[:, dof] = 0
    K[dof, dof] = 1
    F[dof, 0] = 0

""" solving """
coefficients = np.linalg.solve(K, F)
w_nodes  = coefficients[0::2, 0]
phi_nodes = coefficients[1::2, 0]

x_nodes = np.linspace(0.0, L, n + 1)
plt.plot(x_nodes, w_nodes, marker="o")
plt.xlabel("x [m]")
plt.ylabel("w [m]")
plt.title("one-field Timoshenko beam with " + str(n) + " elements")
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

print("L2 norm error (displacement): " + str(err_L2) + " m")

plt.plot(x_nodes, w_nodes, marker="o", label="displacement (one-field)")
plt.plot(x_nodes, tm.w(x_nodes, Load, L, EI), linestyle="--", label="displacement (Euler-Bernoulli)")
plt.xlabel("x [m]")
plt.ylabel("displacement [m]")
plt.title("one-field Timoshenko beam vs Euler-Bernoulli")
plt.grid(True, linestyle="-", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

