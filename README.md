
# Micromechanical Modeling of MEE Composites: The Role of Graphene Nanoplatelets

This repository contains the Abaqus script and MATLAB implementation for the micromechanical modeling of **Magneto-Electro-Elastic (MEE)** multiphase composites. The code investigates the functional role of **Graphene Nanoplatelets (GNPs)** in enhancing the coupled properties of Cobalt Ferrite/Barium Titanate polymer-matrix composites.
In the first step, the mechanical properties of the GNPs/PDMS composite are predicted using the finite element method implemented in Abaqus. Subsequently, spherical piezoelectric and piezomagnetic fillers are incorporated into the composite, and the coupled mechanical and electromagnetic properties are evaluated using the Mori–Tanaka homogenization method in MATLAB.
## 📄 Project Overview
This project utilizes the **Mori-Tanaka mean-field homogenization method** to predict the effective material properties of a complex multi-phase composite consisting of:
1.  **Matrix:** Polymer (e.g., PDMS)
2.  **Piezoelectric Phase:** Barium Titanate (BaTiO₃)
3.  **Piezomagnetic Phase:** Cobalt Ferrite (CoFe₂O₄)
4.  **Reinforcement:** Graphene Nanoplatelets (GNPs)

## 🚀 Key Features
* **Finite Element Simulation:** predict mechanical properties of GNP/PDMS nanocomposites using finite element method in Abaqus sotware.
* **Mori-Tanaka Homogenization:** Implementation of the multi-phase Mori-Tanaka scheme for predicting effective moduli.
* **Eshelby Tensor Calculation:** Accurate calculation of Eshelby tensors for spherical piezoelectric and piezomagnetic fillers.
* **Parametric Study:** Analysis of the effect of piezoelectric and piezomagnetic volume fractions, GNPs volume fraction, thickness, distribution types and effect of interphase layer on:
    * Elastic Stiffness ($C_{ij}$)
    * Piezoelectric Coefficients ($e_{ij}$)
    * Piezomagnetic Coefficients ($q_{ij}$)

## 🛠️ How to Run
1.  Clone this repository.
2.
 
4.  Run Dastgir_magnetomechanics.m the script to generate property plots (Stiffness, Piezo-coefficients, etc.) vs. Piezomagnetic Volume Fraction.
5.  Run Dastgir_electric_properties.m the script to generate property plots (Stiffness, Piezo-coefficients, etc.) vs. Piezoelectric Volume Fraction.
6.  The excel file named (result of GNP vol on magnetimechanics properties.xlsx) presents the outputs of mechanical and magnetic properties of composite as function of the volume fractions of piezomagnetic and GNP fillers. For furthur detail, please refer to Figure 6 and 7 of the article.
7.  The excel file named (result of GNP vol on electric properties.xlsx) presents the outputs of electric properties of composite as function of the volume fractions of piezoelectric and GNP fillers. For furthur detail, please refer to Figure 7 of the article.

## 📚 References
This code is based on the material properties and theoretical framework presented in:
* *"Micromechanical modeling of the functional role of graphene nanoplatelets in the coupled properties of magnetoelectric cobalt ferrite/barium titanate polymer-matrix composites"*
