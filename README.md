# EulerianFluidSimulation
Physically-based Eulerian (Grid-based) Fluid Simulation: Smoke Simulator  
Eulerian Smoke Simulation & Mathematical Foundation: <[FLUID SIMULATION SIGGRAPH 2007 Course Notes](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)>

### Installation
  - *libpng*: http://libpng.org/pub/png/libpng.html
  - *zlib*: http://zlib.net/
  
### Usage
  - Enter:
  ```
  cd EulerianFluidSimulation
  ```
  - Build:
    - Windows (*cygwin*):
      ```
      make -fMakefile.cygwin
      ```
    - Linux:
      ```
      make -fMakefile.linux
      ```
    - Mac OS:
      ```
      make -fMakefile.mac
      ```  
  - Run:
    ```
    ./FLUID_3D
    ```

### Result (.png)
  - Obstacle: ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0030.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0040.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0050.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0060.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0070.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0080.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/with_density_fullxy_0090.png)
  - No obstacle: ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0030.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0040.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0050.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0060.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0070.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0080.png)
  ![image](https://raw.githubusercontent.com/ironyoung/EulerianFluidSimulation/master/result/without_density_fullxy_0090.png)
