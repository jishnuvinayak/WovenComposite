# Woven Composite Material Implementation for Monolithic FE² Framework

## Overview

This repository contains the computational implementation developed as part of my Master’s thesis in Computational Materials Science at TU Bergakademie Freiberg.

The objective of this work was the integration of a woven composite material definition into an efficient monolithic FE² framework. The implementation enables the coupled simulation of yarn and matrix phases within a multi-scale finite element environment.

---

## Scientific and Computational Contributions

The work comprises:

* Implementation of a **UMAT-like material interface** to combine matrix and yarn material models within the FE² framework.
* Automatic switching between yarn and matrix material routines depending on the active macroscopic material point.
* Handling of **local yarn orientation mapping** at each integration point.
* Development of an **Abaqus plugin (RSG-based)** for assignment of woven composite material properties and orientation data.
* Numerical verification using single-element and single-material-point test cases.
* RVE simulations to evaluate the combined response of yarn–matrix systems.
* Benchmark validation against experimental results at the macroscopic level.
* Parallelization study investigating scalability of the monolithic FE² framework.

---

## Numerical Aspects

The implementation addresses:

* Coupled multi-scale finite element analysis using a monolithic FE² approach.
* Constitutive integration of matrix plasticity models based on return mapping algorithms.
* Investigation of convergence behavior under automatic time stepping.
* Analysis of algorithmic limitations arising from the use of continuum tangents instead of algorithmically consistent material tangents.
* Scalability assessment for realistic macrostructure simulations.

---

## Limitations and Observations

* The current implementation considers a ±90° stacking sequence.
* Matrix-dominated loading cases (e.g., ±45° stacking) exhibit convergence issues due to the return mapping formulation.
* The cutting-plane based algorithm does not provide an algorithmically consistent tangent operator, which affects convergence for larger time steps.
* Future improvements require refinement of the return mapping scheme and improved initialization strategies for the local Newton–Raphson iteration.

---

If you would like, I can now:

* Add a **Repository Structure** section,
* Add a short **How to Run** section,
* Or slightly tighten the wording to make it more “framework developer–oriented” for the Hereon reviewers.

