# Installation and Run Instructions

## Installation Instructions

We have a set of easy instructions for users running Ubuntu 18.04+ (including WSL
with Ubuntu 18.04, as well as newer Ubuntu LTS)

Easy scripts - [Easy Linux instructions](./Install_ubuntu_easy.md)

If the automated installation of dependencies fails, or if you don't have a system
supporting it, then follow the instructions below: 

For Linux machines - [Linux installation instructions](./Install_linux.md)  
For MacOS machines - [MacOS installation instructions](./Install_macos.md)


## Folder Layout
Now that the code has been installed, let us comment on the subdirectory structure:

- [ChiDoc](./ChiDoc) contains 
  - scripts to generate the code documentation (these scripts were actually used during the installation instructions, above),
  - the generated documentation;
- [ChiTech](./ChiTech) contains the source code that is agnostic of the physics, such as math, mesh, MPI, ...
- [ChiModules](./ChiModules) contains the source code for specialized physics modules,
- [ChiTest](./ChiTest) contains tests for the regression suite,
  - TODO: Jean to re-organize this folder.
- [ChiResources](./ChiResources) contains ???
  - TODO: clean-up needed. Decide where to place TestObjects
- [ThirdParty](./ThirdParty) contains external library source code, such as VTK
- [Tutorials](./Tutorials) contains some tutorials. 
  - TODO: this seems sparse. I want to add more user tutorial examples and created more dedicated pages for that.

TODO: At some point, rename folders to doc, modules, resources, src, tests, external

## Run Instructions

###  Serial run (interactive)
Simply type:
```bash 
path_to_exec/ChiTech  path_to_input_file/input_filename.lua
```
where 
- ```path_to_exec``` is the path to the directory where the executable was created and saved,
- ```ChiTech``` is the (default) filename for the executable,
- ```path_to_input_file``` is the path to the directory where the input file is located,
- ```input_filename.lua``` is the name of the LUA input file

It is also possible to run in batch mode on a cluster or a supercomputer. See below.

### Parallel run (interactive)
In interactive mode (i.e., on your own machine), simply type
```bash 
mpiexec -n N path_to_exec/ChiTech  path_to_input_file/input_filename.lua
```
where ```N``` is the requested number of parallel processes.

### Parallel run (batch)

In batch mode, you have to use the job submission of the cluster you intend 
on using.

For example, a slurm submission script may look like this:
```bash
slurm example
```
