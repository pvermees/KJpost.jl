# PTpost.jl

## Post-processing functions for Plasmatrace

Plasmatrace is a free and open data reduction software package for
LA-ICP-MS written in [Julia](https://julialang.org/). PTpost is an
extension that adds post-processing functionality to Plasmatrace.  At
the moment this is limited to internal-isochron regression. Further
functions will be added in the future.

## Installation

Enter the following commands at the Julia console (a.k.a REPL):

```
import Pkg
Pkg.add(url="https://github.com/pvermees/Plasmatrace.jl.git")
Pkg.add(url="https://github.com/pvermees/PTpost.jl.git")
```

## Minimal working example

```
julia> using Plasmatrace, PTpost
julia> PT(PTpost)
-------------------
 Plasmatrace 0.7.4
-------------------

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