# LP2FAD
MATLAB implementation of the master thesis work carried out at the Faculty of Aerosapce Engineering, TU Delft.

# Master Thesis Details
Stiffness Design of Laminated Composites: Efficient Conversion of Lamination Parameters into Stacking Sequences

# Contributors
Master Thesis Student: Rakshith Manikandan\
Supervisors: Dr.Ir. D.M.J. Peeters, and Dr.ir J.M.J.F. van Campen

# Motivation
Converting a set of Lamination Parameters (LPs) into Stacking sequences (SS) is usually performed using a Genetic Algorithm or Decision-tree based tools. However, as the number of layers and permissible ply orientations increases, these methods can become computationally expensive, posing a significant challenge. Hence, this thesis proposed partitioning the process of converting LPs into SS. First, the number of plies belonging to every ply orientation, or the Fibre Angle Distribution (FAD), will be designed. Then, using them as a basis, the SS will be designed. \
\
As a first step in this direction, this master's thesis focuses on efficiently designing FADs, for which a novel Fast Fourier Transform(FFT)-based method is introduced. The MATLAB implementation of this method can be found in this repository. Future updates of this code will focus on converting FADs into SS.

# Repository
### This repository consists of the following:
-The function LP2FAD.m: A novel-FFT based implementation used to swiftly convert In-Plane Lamination Parameters into a laminate's Fibre Angle Distribution.\
-The code LP2FAD_Example.m: An example code on how to use the function 'LP2FAD.m'.\
-The 2D benchmark: Dataset with 2D targets of V1A and V3A Lamination Parameters. They were uniformly sampled from the feasible In-Plane Lamination Parameter space (for V1A and V3A).\
-The 4D benchmark: Dataset with 4D targets of V1A,V2A,V3A and V4A Lamination Parameters. They were quasi-randomly sampled from the feasible In-Plane Lamination Parameter space (for V1A,V2A,V3A, and V4A).\
-The Bristol Dataset: Dataset with 4D targets of V1A,V2A,V3A and V4A Lamination Parameters. They belong to real laminates made of 200 layers (in total) and consist of the orientations [0,45,90]. Sourced from (https://github.com/noemiefedon/LAYLA/tree/main). \
-The 4D benchmark: Dataset with 4D targets of V1A,V2A,V3A and V4A Lamination Parameters. They were quasi-randomly sampled from the feasible In-Plane Lamination Parameter space (for V1A,V2A,V3A, and V4A).\

# License
LP2FAD © 2023 by Rakshith Manikandan is licensed under CC BY-SA 4.0.\
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/

# Validation
This code has been validated with benchmark datasets available in the literature and custom-made datasets for more thorough tests. Nevertheless, the possibility of errors to appear is not null. In such cases, please do not hesitate to contact me (rakshith.m1505@gmail.com) to report any problems/fixes/suggestions for improvement.
