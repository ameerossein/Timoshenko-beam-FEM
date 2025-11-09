Timoshenko Beam Finite Element Implementation — Homework 1

This repository contains the source code and report files for Homework 1 of the course
Nonlinear Finite Element Methods — WiSe 2025/26.

The task involves implementing one-field and two-field finite-element formulations for a cantilever Timoshenko beam and comparing results to the analytical Euler–Bernoulli beam solution.

📘 Problem Summary

A steel cantilever beam (L = 5 m) with IPE100 section is loaded by a point force at its free end.
The assignment requires:

Derivation of Euler–Bernoulli analytical solution

One-field Timoshenko variational formulation (standard FE)

Two-field Timoshenko variational formulation (Hellinger–Reissner)

Discretization using 2-node linear finite elements

Convergence studies (L2 norm of displacement)

Comparison of formulations

📂 Files
File	Description
timo.py	Analytical Euler–Bernoulli reference solution
onefield.py	One-field Timoshenko FE implementation + convergence study
twofield.py	Two-field Timoshenko FE implementation + convergence study
comparison.py	Summary statements comparing methods
homework.pdf	Problem description / assignment sheet
🚀 How to Run

Install dependencies (optional):

pip install numpy matplotlib


Run scripts:

python onefield.py
python twofield.py
python comparison.py

📈 Output

The scripts generate:

Displacement plots for both FE formulations

Convergence curve (L2 displacement error)

Comparison with Euler–Bernoulli analytical solution

✅ Key Observations (Summary)

One-field Timoshenko: simple, shear locking for slender beams

Two-field (Hellinger–Reissner): locking-free, more unknowns → more expensive

Euler–Bernoulli: accurate for slender beams, neglects shear deformation

🎓 Course

Nonlinear Finite Element Methods
WiSe 2025/26

This repository serves academic purposes and demonstrates FEM beam formulations and numerical convergence behavior.
