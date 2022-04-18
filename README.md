# MultiscaleGraphSignalTransforms.jl
| Doc | Build | Test |
|------|-------|------|
| [![](https://img.shields.io/badge/docs-passing-success)](https://ucd4ids.github.io/MultiscaleGraphSignalTransforms.jl/dev/) | [![CI](https://github.com/UCD4IDS/MultiscaleGraphSignalTransforms.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/UCD4IDS/MultiscaleGraphSignalTransforms.jl/actions) | [![codecov](https://codecov.io/gh/UCD4IDS/MultiscaleGraphSignalTransforms.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/UCD4IDS/MultiscaleGraphSignalTransforms.jl) |


![Haar-Walsh on R vs on Graph](./GHWT.png "Haar-Walsh on R vs on Graph")

## COPYRIGHT

Copyright 2015-2021 The Regents of the University of California

Implemented by Jeff Irion, Haotian Li, Naoki Saito, and Yiqun Shao


## SETUP

To install the MultiscaleGraphSignalTransforms.jl, run
```julia
julia> import Pkg; Pkg.add("MultiscaleGraphSignalTransforms")
julia> using MultiscaleGraphSignalTransforms
```

## GETTING STARTED

Currently, you can run a set of very small tests via ```] test MultiscaleGraphSignalTransforms```; see the actual file ```test/runtest.jl``` as well as [the documentation](https://ucd4ids.github.io/MultiscaleGraphSignalTransforms.jl/dev) for further information.

## REFERENCES

1. J. Irion and N. Saito, [Hierarchical graph Laplacian eigen transforms](https://www.math.ucdavis.edu/~saito/publications/hglets.html), *Japan SIAM Letters*, vol. 6, pp. 21-24, 2014.

2. J. Irion and N. Saito, [The generalized Haar-Walsh transform](https://www.math.ucdavis.edu/~saito/publications/ghwt.html), *Proc. 2014 IEEE Statistical Signal Processing Workshop*, pp. 488-491, 2014.

3. J. Irion and N. Saito, [Applied and computational harmonic analysis on
graphs and networks](https://www.math.ucdavis.edu/~saito/publications/spie15.html), *Wavelets and Sparsity XVI*, (M. Papadakis, V. K. Goyal, D. Van De Ville, eds.), *Proc. SPIE 9597*, Paper #95971F, Invited paper, 2015.

4. J. Irion, [Multiscale Transforms for Signals on Graphs: Methods and Applications](https://jefflirion.github.io/publications_and_presentations/irion_dissertation.pdf), Ph.D. dissertation, University of California, Davis, Dec. 2015.

5. J. Irion and N. Saito, [Learning sparsity and structure of matrices with multiscale graph basis dictionaries](https://www.math.ucdavis.edu/~saito/publications/matanal.html), *Proc. 2016 IEEE 26th International Workshop on Machine Learning for Signal Processing (MLSP)*, (A. Uncini, K. Diamantaras, F. A. N. Palmieri, and J. Larsen, eds.), 2016.

6. J. Irion and N. Saito, [Efficient approximation and denoising of graph signals using the multiscale basis dictionaries](https://www.math.ucdavis.edu/~saito/publications/eadgsumbd.html), *IEEE Transactions on Signal and Information Processing over Networks*, Vol. 3, no. 3, pp. 607-616, 2017.

7. N. Saito, [How can we naturally order and organize graph Laplacian eigenvectors?](https://www.math.ucdavis.edu/~saito/publications/lapeigport.html) *Proc. 2018 IEEE Workshop on Statistical Signal Processing*, pp. 483-487, 2018.

8. Y. Shao and N. Saito, [The extended Generalized Haar-Walsh Transform and applications](https://www.math.ucdavis.edu/~saito/publications/eghwt.html), *Wavelets and Sparsity XVIII*, (D. Van De Ville, M. Papadakis, and Y. M. Lu, eds.), *Proc. SPIE 11138*, Paper #111380C, 2019.

9. H. Li and N. Saito, [Metrics of graph Laplacian eigenvectors](https://www.math.ucdavis.edu/~saito/publications/metgraphlap.html), *Wavelets and Sparsity XVIII*, (D. Van De Ville, M. Papadakis, and Y. M. Lu, eds.), *Proc. SPIE 11138*, Paper #111381K, 2019.

10. Y. Shao, [The Extended Generalized Haar-Walsh Transform and Applications](https://www.math.ucdavis.edu/~tdenena/dissertations/202008_Shao_Yiqun_dissertation.pdf), Ph.D. dissertation, University of California, Davis, Sep. 2020.

11. C. Alexander, H. Li and N. Saito, [Natural graph wavelet packet dictionaries](https://www.math.ucdavis.edu/~saito/publications/ngwp.html), *J. Fourier Anal. Appl.*, vol. 27, Article \#41, 2021.

12. H. Li, Natural Graph Wavelet Dictionaries: Methods and Applications, Ph.D. dissertation, University of California, Davis, Jun. 2021.

13. N. Saito and Y. Shao, [eGHWT: The extended Generalized Haar-Walsh Transform](https://www.math.ucdavis.edu/~saito/publications/eghwt21.html), *J. Math. Imaging Vis.*, vol. 64, no. 3, pp. 261-283, 2022.
