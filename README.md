# Metamodeling and Control of Medical Digital Twins

This repository contains the files used in the manuscript entitled: Metamodeling and Control of Medical Digital Twins.<br>
Under review in *Science Advances*. Preprint at [arXiv](https://doi.org/10.48550/arXiv.2402.05750). <br>

The manuscript uses two agent-based models (ABMs) to illustrate the different methodologies for creating ODE approximations of ABMs. The two ABMs are the sheep-wolves-grass version of the [Wolf Sheep Predation model](https://ccl.northwestern.edu/netlogo/models/WolfSheepPredation) and a purposely built model of a toy metabolic pathway.

## Files used for approximating the Wolf Sheep Predation model
 - `Wolf Sheep Predation original.nlogo` **Original** code for the Wolf Sheep Predation model.
 - `Wolf Sheep Predation bigworld.nlogo` Netlogo file used to generate **dataset I**.
 - `Wolf Sheep Predation bigworld 2nd DataSet.nlogo` Netlogo file used to generate **dataset II**.
 -  `Wolf Sheep Predation bigworld_ConGrass2.nlogo` Netlogo file used to generate **dataset III** obtained by removing `2%` of grass per time-step.
 -  `Wolf Sheep Predation bigworld_ConSheep2.nlogo` Netlogo file used to generate **dataset IV** obtained by removing `2%` of sheep per time-step.
 -  `Wolf Sheep Predation bigworld_ConWolves1.5.nlogo` Netlogo file used to generate **dataset V** obtained by removing `1.5%` of wolves per time-step.
 -  'Dados01.mat' MAT-file containing dataset I time series. This dataset resulted from **averaging 100 simulations** of the `Wolf Sheep Predation bigworld.nlogo` model.
 -  `Dados02.mat` MAT-file containing dataset II time series. This dataset resulted from **averaging 100 simulations** of the `Wolf Sheep Predation bigworld 2nd DataSet.nlogo` model.
 -  `DadosConGrass2.mat` MAT-file containing dataset III time series. This dataset resulted from **averaging 100 simulations** of the `Wolf Sheep Predation bigworld_ConGrass2.nlogo` model.
 -  `DadosConGrass2s.mat` MAT-file containing standard deviation associated with dataset III time series. Used only for plotting.
 -  `DadosConSheep2.mat` MAT-file containing dataset IV time series. This dataset resulted from averaging 100 simulations of the `Wolf Sheep Predation bigworld_ConSheep2.nlogo` model.
 -  `DadosConSheep2s.mat` MAT-file containing standard deviation associated with dataset IV time series. Used only for plotting.
 -  `DadosConWolves1.5.mat` MAT-file containing dataset V time series. This dataset resulted from averaging 100 simulations of the `Wolf Sheep Predation bigworld_ConWolves1.5.nlogo` model.
 -  `DadosConWolves1.5s.mat` MAT-file containing standard deviation associated with dataset V time series. Used only for plotting.
 -  `SWG_Case1_Mech_I_II.m` MATLAB script containing the model and best-fit parameter-set for the **mechanistic** approximation (Case 1) fitted against **datasets I and II**.
 -  `SWG_Case1_Mech_I_V.m` MATLAB script containing the model and best-fit parameter-set for the **mechanistic** approximation (Case 1) fitted against **datasets I, II, III, IV, and V**.
 -  `SWG_Case2_GMA_I_II.m` MATLAB script containing the model and best-fit parameter-set for the **GMA** approximation (Case 2) fitted against **datasets I and II**.
 -  `SWG_Case2_GMA_I_V.m` MATLAB script containing the model and best-fit parameter-set for the **GMA** approximation (Case 2) fitted against **datasets I, II, III, IV, and V**.
 -  `SWG_Case3_Linear_I.m` MATLAB script containing the model and best-fit parameter-set for the **first-order polynomial** approximation (Case 3) fitted against **dataset I**.
 -  `SWG_Case3_Quad_I_II.m` MATLAB script containing the model and best-fit parameter-set for a **second-order polynomial** approximation (Case 3) fitted against **datasets I and II**.
 -  `SWG_Case3_Quad_I_V.m` MATLAB script containing the model and best-fit parameter-set for a **second-order polynomial** approximation (Case 3) fitted against **datasets I, II, III, IV, and V**.
 -  `SWG_Case4_Ssystem_I_II.m` MATLAB script containing the model and best-fit parameter-set for the **S-system** approximation (Case 4) fitted against **datasets I and II**.
 -  `SWG_Case4_Ssystem_I_V.m` MATLAB script containing the model and best-fit parameter-set for the **S-system** approximation (Case 4) fitted against **datasets I, II, III, IV, and V**.

  ## Files used for approximating the Metabolic Pathway model
 -  `Met_Pathway_dataset_80k_20k_20k_10_10_NoDil.m` MATLAB script for the Metabolic Pathway model operating under **Batch mode**.
 -  `Met_Pathway_dataset_80k_20k_20k_10_10_wDil_wFeed.m`  MATLAB script for the Metabolic Pathway model operating under **Continuous mode** with an inflow of `1.0`.
 -  `Met_Pathwayv2_S80k_P20k_Q20k_noDil.mat` MAT-file containing **dataset I** time series. This dataset resulted from a single simulation of the Metabolic Pathway model in **Batch mode**, obtained using `Met_Pathway_dataset_80k_20k_20k_10_10_NoDil.m`.
 -  `Met_Pathwayv2_S80k_P20k_Q20k_Dil0005In1.mat` MAT-file containing **dataset II** time series. This dataset resulted from a single simulation of the Metabolic Pathway model in **Continuous mode** with an inflow of `1.0`, obtained using `Met_Pathway_dataset_80k_20k_20k_10_10_wDil_wFeed.m`.
 -  `MetPw_TrainingDatasets.mat` MAT-file containing **datasets III, IV, and V** time series. These datasets were obtained by averaging 100 simulations. **Dataset III** - same conditions as dataset I (**Batch mode**).  **Dataset IV** - same conditions as dataset II (**Continuous mode** with an inflow of `1.0`).  **Dataset V** - similar conditions as dataset II, in **Continuous mode**, but with an inflow of **`0.2`**.
 -  `MetPw_Case1_Mech_I.m` MATLAB script containing the model and best-fit parameter-set for the **mechanistic** approximation (Case 1) fitted against **datasets I and II**.
 -  `MetPw_Case1_Mech_C.m` MATLAB script containing the model and best-fit parameter-set for the **mechanistic** approximation (Case 1) fitted against **datasets III, IV, and V**.
 -  `MetPw_Case2_GMA_I.m` MATLAB script containing the model and best-fit parameter-set for the **GMA** approximation (Case 2) fitted against datasets **I and II**.
 -  `MetPw_Case2_GMA_C.m` MATLAB script containing the model and best-fit parameter-set for the **GMA** approximation (Case 2) fitted against **datasets III, IV, and V**.
 -  `MetPw_Case3_Linear_I.m` MATLAB script containing the model and best-fit parameter-set for the **first-order polynomial** approximation (Case 3) fitted against **datasets I and II**.
 -  `MetPw_Case3_Quad_I.m` MATLAB script containing the model and best-fit parameter-set for a **second-order polynomial** approximation (Case 3) fitted against **datasets I and II**.
 -  `MetPw_Case4_Ssystem_I.m` MATLAB script containing the model and best-fit parameter-set for the **S-system** approximation (Case 4) fitted against **datasets I and II**.
 -  `MetPw_Case4_Ssystem_C.m` MATLAB script containing the model and best-fit parameter-set for the **S-system** approximation (Case 4) fitted against **datasets III, IV, and V**.
  
