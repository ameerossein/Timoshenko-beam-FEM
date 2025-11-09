""" task 1.6 """
#%% one-field vs two-field
print("Timoshenko one-field weak form captures both bending and shear deformations. Only C0 continuity required"
      " On the other hand, it has shear locking for slender beams meaning that the beam behaves too stiff."
      " Consequently, it does not perform really well with coarse mesh."
      " However, Timoshenko two-field (HR) is locking free and is very accurate even with coarse mesh"
      " But it requires more DOFs (4 per node) which results in expensive computation.")

#%% pros and cons of Euler-Bernoulli
print("Although Euler-Bernoulli theory is very simple and exact for slender beams, it cancels out shear deformation"
      " which makes it inaccurate for short or thick beams. Besides, it requires a C1 continuity for exact theory.")