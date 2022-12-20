![GitHub branch checks state](https://img.shields.io/github/checks-status/chi-tech/chi-tech/development?label=Development)

# Chi-Tech #

A large-scale scientific simulation engine being developed at Texas A&M University 
as part of a research study.

Chat with us on <a href="https://discord.gg/AJHHeAA">Discord<img src="ChiDoc/HTMLimages/DiscordLogo.png" width="24" title="Chi-Tech on Discord" align="center"></a>

## Installation instructions

We have a set of easy instructions for users running Ubuntu 18.04+ (including WSL
with Ubuntu 18.04, as well as newer Ubuntu LTS)

Easy scripts - [Easy Linux instructions](ChiDoc/Install_ubuntu_easy.md)

If the automated installation of dependencies fails, or if you don't have a system
supporting it, then follow the instructions below: 

For Linux machines - [Linux installation instructions](ChiDoc/Install_linux.md)  
For MacOS machines - [MacOS installation instructions](ChiDoc/Install_macos.md)

## Documentation

The preferred method of accessing the documentation is to clone the repo and
build the doxygen-based documentation. However the documentation
is also periodically updated at the following link:

[Online Documentation](https://chi-tech.github.io)

# Showcase
## 1. Capable of sweeps on Polyhedral meshes
Sphere embedded within a box:
- Concave cells induce cyclic dependencies between cells.
- Non-cutting Volumetric partitioning induces cyclic dependencies
  between processors.
- Data structures allows for sweeping will all forms of cycles.

![yes](ChiDoc/HTMLimages/CoolPics/SOD_threshold.png)
3D Polyhedral mesh generated with STAR-CCM+.
![yes](ChiDoc/HTMLimages/CoolPics/SOD_slice.png)
Slice of the solution.
![yes](ChiDoc/HTMLimages/CoolPics/SOD_partitioning.png)
Arbitrary, non-cutting, KBA-style partitioning used.

## 2. C5G7 Criticality Benchmark with 768 processors
The famous reactor benchmark C5G7:
- 7 energy groups
- 200 directions
- 454,491 cells
- Ran on 768 processors
- Took only 18 minutes to complete
- Used 584 GB of memory
- `k_eff` within 100 pcm


![yes](ChiDoc/HTMLimages/CoolPics/C5G7_materials.png)
Closeup view of the mesh used. Colors represent materials.
![yes](ChiDoc/HTMLimages/CoolPics/C5G7_group0.png)
Energy group 0 solution.
![yes](ChiDoc/HTMLimages/CoolPics/C5G7_group6.png)
Energy group 6 solution.
![yes](ChiDoc/HTMLimages/CoolPics/C5G7_partition768.png)
ParMETIS partitioning of the mesh (768 processors).
![yes](ChiDoc/HTMLimages/CoolPics/C5G7_partition768b.png)
Closeup of the ParMETIS partitioning with the mesh visible.

## 2. Real world simulations
Center for Exascale Radiation Transport (CERT) simulated, and
compared to experiment, a graphite pile with a high energy neutron
source. This simulation used:
- ~172 energy groups. 
- over 3000 directions. 
- ~500k cells.
- Over 100k processors for some simulations.

![yes](ChiDoc/HTMLimages/CoolPics/CERTSim.png)