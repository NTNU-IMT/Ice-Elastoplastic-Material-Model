# Ice-Elastoplastic-Material-Model
This repository includes VUMAT and VUSDFLD subroutines for simulating the elastoplastic behaviour of ice crushed in ice-structure collision simulations. They are both written in Fortran for the Abaqus/Explicit solver (v. 2019). The VUMAT is based on the material model initially proposed by Liu et al. in 2011. The VUSDFLD should be used in conjunction with the Abaqus built-in Crushable Foam (CF) plasticity model to simulate the ice failure. Both plasticity models (i.e., VUMAT and CF) are based on the Tsai-Wu yield surface. However, the VUMAT plasticity model is governed by an associated flow rule while the CF model uses a non-associated flow rule with essentially zero plastic Poisson’s ratio. This zero plastic Poisson’s ratio was found to cause serious problems in ice crushing simulations. The details and the constitutive laws are discussed in the paper by Mokhtari et al. (2022) available in this repository and published in Marine Structures: https://doi.org/10.1016/j.marstruc.2022.103233. This paper has won the Moan-Faltinsen best paper award.
In the VUMAT subroutine, there is a damage model that has not been used in the simulations discussed in the paper noted above. It means that you should input 0 for Gf1, Gf2 and alphad parameters when implementing the code unless you want to test/use the damage model. In the case of using the damage model, appropriate references should be made to the master thesis by Xintong Wang where the damage model is proposed: “Analysis of Iceberg-Structure Interaction During Impacts”.
The VUMAT and VUSDFLD scripts were developed under the supervision of Professors Jørgen Amdahl and Ekaterina Kim whom I greatly appreciate. The financial support of the Research Council of Norway through the Centers of Excellence funding scheme, project AMOS (Grant number 223254) and the Centers for Research-based Innovation funding scheme, project CASA (Grant number 237885), at the Norwegian University of Science and Technology are acknowledged. 
