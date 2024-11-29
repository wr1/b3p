# Introduction

*b3p* is a preprocessor for wind turbine blade models. 

``` mermaid 
graph TD
    A[Blade input file] --> B[b3p build]
    B --> FEA[3D FEA]
    B --> 2D[Section meshes]
    B --> BEM 
```

 <!-- python package that translates a blade input file to allow running various solvers such as *ccx*, *anba* or *ccblade*. -->


