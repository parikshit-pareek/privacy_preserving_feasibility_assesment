# privacy_preserving_feasibility_assesment

This file contain the code related to the paper: 

Parikshit Pareek and Hung Nguyen, ["Privacy-Preserving Feasibility Assessment for P2P Energy Trading and Storage Integration" Accepted for 2022 PESGM"](https://www.researchgate.net/publication/358660003_Privacy-Preserving_Feasibility_Assessment_for_P2P_Energy_Trading_and_Storage_Integration), [Preprint](https://www.researchgate.net/publication/358660003_Privacy-Preserving_Feasibility_Assessment_for_P2P_Energy_Trading_and_Storage_Integration)

Cite As: 
```
@article{pareek2022privacy,
  title={Privacy-Preserving Feasibility Assessment for P2P Energy Trading and Storage Integration,
  author={Pareek, Parikshit and Singh, Anshuman and Sampath, LP Mohasha Isuru and Gooi, HB and Nguyen, Hung D},
 booktitle={2022 IEEE Power Engineering Society General Meeting--Accepted},
  pages={1--5},
  year={2022},
  organization={IEEE}
}
```
In perticular, the code generates figure 2 of the manusript and provides the all the results in cell array **R**. 

## Details of Files: 
- `CFPF_QD_Kernel_DER_bus.m` : Function obtaining the CFPF Approximation using Quadratic Kernel.
- `input_dataset_Load.m` : Creating load data set for training and testing 
- `MCS_output.m`    : Monte-Carlo Simulation to obtain testing data points
- `rand_sample_x.m` : Generating random samples
- `runpf_complete.m` : MATPOWER codes combined together to avid downloading 
- `Sampling_Jaco.m`  : Cover to 'runpf' for obtaining power flow datasets
- `main_final.m`     : Main file to run


## Results
- `R{1, hr_idx}.f`  : Index of probabilisitically feasible injections of battery at the `bus_ESS`node
- `R{1, hr_idx}.Xb_f` Nx2 double : Probabilisitically feasible injections of battery at the `bus_ESS`node  (N is feasible points)
- `R{1, hr_idx}.s_limits` : Dx2 double : Minimum and Maximum injection limits at that hour at each node P and Q (D = Dimension loadbus x 2)
- `R{1, hr_idx}.ind_cap`: Indvidual PV injection (max)
- `R{1, hr_idx}.XPb` : The physical battery space of real power (grid for min and max possible P)
- `R{1, hr_idx}.XQb` : The physical battery space of reactive power (grid for min and max possible Q)

## Dependencies: 
GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION Toolbox version 4.2 for GNU Octave 3.2.x and Matlab 7.x and higher.
Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-06-11.
[Link](http://www.gaussianprocess.org/gpml/code/matlab/doc/)

* Newer GPML toolbox might require matching hyperparameter initialization. To avoid such issues, I have uploaded the GPML Toolbox version [here](https://github.com/parikshit-pareek/Closed-form-Power-Flow/tree/master/GPML). 

- Look [here](https://github.com/parikshit-pareek/Closed-form-Power-Flow) for more details  
- Other variable names are self-explaintory
- All the codes are tested on MATLAB R2020b
