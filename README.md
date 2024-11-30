# KJpost.jl

## Post-processing functions for KJ

[KJ](https://github.com/pvermees/KJ.jl) is a free and open data
reduction software package for LA-ICP-MS written in
[Julia](https://julialang.org/). KJpost is an extension that adds
post-processing functionality to KJ.  At the moment this is limited to
internal-isochron regression. Further functions will be added in the
future.

## Installation

Enter the following commands at the Julia console (a.k.a REPL):

```
import Pkg
Pkg.add(url="https://github.com/pvermees/KJ.jl.git")
Pkg.add(url="https://github.com/pvermees/KJpost.jl.git")
```

## Minimal working example

```
julia> using KJ, KJpost
julia> KJ(KJpost)
----------
 KJ 0.0.1
----------

r: Read data files[*]
m: Specify the method[*]
t: Tabulate the samples
s: Mark mineral standards[*]
g: Mark reference glasses[*]
v: View and adjust each sample
p: Process the data[*]
e: Export the results
l: Logs and templates
o: Options
u: Update
c: Clear
a: Extra
x: Exit
?: Help
a

i: Internal isochron
x: Exit
?: Help
x
julia> 
```