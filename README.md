[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# MagicHelpR

This package was written mostly for personal convenience to help in performing QTL 
mapping in the Arabidopsis multiparent advanced generation intercross (MAGIC) 
lines (described in 
[Kover P .. Mott R 2009](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000551)).

This might however be useful for anyone using the Arabidopsis MAGIC lines. It
depends Richard Mott's `happy.hbrem` package, available at:
http://mtweb.cs.ucl.ac.uk/mus/www/magic/

The package implements the same linear model approach in the available 
scripts, but with some more flexibility for specifying models and outputs 
tabular data that is convenient for downstream analysis and plotting.

In the long term, a better way to carry out these tests is using the 
[R/qtl2 package](https://github.com/rqtl/qtl2), which is under development, but 
now supports analysis of the Arabidopsis MAGIC lines (see 
[here](https://groups.google.com/forum/#!topic/rqtl-disc/pwVi_Igr9zk)). I hope 
to soon implement some functions to convert the MAGIC genotypes to R/qtl2 format.

To install this package, launch R and issue the following commands:

```
#If not already installed, install devtools:
# install.packages("devtools")
devtools::install_github("tavareshugo/MagicHelpR")
```

There's also a basic [introduction vignette](https://tavareshugo.github.io/MagicHelpR//articles/introduction.html).
