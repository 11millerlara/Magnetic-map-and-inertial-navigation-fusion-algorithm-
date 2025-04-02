# AIPF Simulation Project Description

This project is based on the Adaptive Integrated Particle Filter (AIPF) algorithm, which fuses inertial navigation data and magnetic map data to achieve matching-based positioning in simulation. The project mainly consists of the following parts:

## 1. Main Functions
- **D2_igrf_pos_main.m**: Simulation verification program based on the IGRF model.  
- **D2_indoor_main.m**: Simulation verification program using indoor measurement data.  
- **D2_33test_main.m**: Simulation verification program based on UAV marine magnetic survey data.  

## 2. multi_algorithm Folder
Contains multiple particle filter algorithms for comparative experiments, including:
- Particle Filter (PF)  
- Extended Kalman Particle Filter (EKPF)  
- Adaptive Firefly Algorithm (AOFA)  
- Adaptive Integrated Particle Filter (AIPF)  

## 3. magnetic_map_data Folder
Stores magnetic map data interpolated by Gaussian Process Regression (GPR) for positioning and matching.

## 4. BPIT Folder
Contains the BPIT table built upon the standard normal distribution, used for table lookups in the AIPF algorithm.

## 5. cal_pos Folder
Provides functions for generating simulated test trajectories, offering data points for algorithm validation.

## 6. paper Folder
Holds papers and related materials published based on the experimental results of this project for reference and further research.

---

## Usage
1. Download the entire project and place it in your MATLAB working directory.  
2. In the root directory of the project, run the main function for the desired scenario (e.g., D2_igrf_pos_main.m) to start the simulation.  
3. Logs and results will be displayed in the MATLAB Command Window and in plot windows, calling the related functions or data from various folders as needed.

For users who need to modify or extend the algorithm, you can do so directly in the `multi_algorithm` folder and reference your updated scripts in the main functions.