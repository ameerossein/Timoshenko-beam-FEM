""" Analytical Euler-Bernoulli beam """
import numpy as np
import matplotlib.pyplot as plt

""" properties and data """
EI = 200*(10**9)*171*(10**(-8))
L = 5
Load = -1

""" defining w and theta """
def w(x, Load, L, EI):
    return Load * x**2 * (3*L - x) / (6 * EI)

def theta(x, Load, L, EI):
    return Load * x * (2*L - x) / (2 * EI)

""" BCs """
w_tip = w(L, Load, L, EI)
theta_tip = theta(L, Load, L, EI)

""" calculate the fields """
x = np.linspace(0.0, L, 11)
displacement = w(x, Load, L, EI)
rotation = theta(x, Load, L, EI)
#print(displacement)

""" plot the displacement """
plt.plot(x, displacement)
plt.xlabel("x [m]")
plt.ylabel("Displacement w [m]")
plt.title("Analytical Euler–Bernoulli solution")
plt.grid(True, linestyle="-", alpha=0.5)
plt.show()

