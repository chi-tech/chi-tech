# User Tutorials

## Lua for users

ChiTech input filenames are written in Lua. 
- Some of the code objects are directly accessible in the Lua console
- Portion of the input can often be generated using Lua scripting

Hence, we recommend the users to familiarize themselves with Lua. 
See the [official lua page](https://www.lua.org/) for more details.

TODO: A table of all existing lua functions accessible from the lua input file would be good here.

## Meshing

ChiTech handles many types of meshes, including arbitrary polyhedral meshes. 
Complex meshes are often generated externally and loaded into ChiTech. 
Hence, it is important for users to be familiar with CAD and mesh generators.

We provide some [meshing tutorial here](https://chi-tech.github.io/d4/d25/_mesh_tutorials.html) 
but recommend users seek additional knowledge when needed.

## 
A meshing section, including using most of the features currently in chi-tech and a link to our other existing meshing tutorials.
i.      I am also considering porting the existing meshing tutorials how of the doxygen and have them in markdown. Why? I prefer the doxygen to be only the source code documentation
For the same reasons as above, moving out the 3 explained tutorials out of doxy and into markdown. Adding a lot more tutorials, including explanation for existing regression files and adding input files created by Ian and Patrick during their dissertation. Since this can easily add up to several dozens of examples, organizing them with connection and level of complexity as is done in deal.ii could be worth the investment.