[![Project Status: Abandoned – Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author(s) do not intend on continuing development.](https://www.repostatus.org/badges/latest/abandoned.svg)](https://www.repostatus.org/#abandoned)
[![Build Status](https://travis-ci.org/tavareshugo/MagicHelpR.svg?branch=master)](https://travis-ci.org/tavareshugo/MagicHelpR)

:warning: I would strongly advise against using this package! 
Instead, you can use the [R/qtl2](https://kbroman.org/qtl2/) for analysis of MAGIC data. 
For Arabidopsis specifically, I've written a data package ([`atMAGIC`](https://github.com/tavareshugo/atMAGIC)), 
which provides a `cross2` object to use with that package. Alternatively there's also data 
available from the [R/qtl2data](https://github.com/rqtl/qtl2data/tree/master/ArabMAGIC) repository.

# MagicHelpR

This package was written mostly for personal convenience to help in performing QTL 
mapping in the Arabidopsis multiparent advanced generation intercross (MAGIC) 
lines (described in 
[Kover P .. Mott R 2009](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000551)).

This might however be useful for anyone using the Arabidopsis MAGIC lines. It
depends on [the `happy.hbrem` package](http://mtweb.cs.ucl.ac.uk/mus/www/magic/), 
which can be installed from https://github.com/tavareshugo/happy.hbrem

The package implements the same linear model approach in the available 
scripts, but with some more flexibility for specifying models and outputs 
tabular data that is convenient for downstream analysis and plotting.

To install this package, launch R and issue the following commands:

```r
if(!require("remotes")) install.packages("remotes")
if(!require("happy.hbrem")) remotes::install_github("tavareshugo/happy.hbrem")

remotes::install_github("tavareshugo/MagicHelpR")

# If the above fails (it has done for me in Ubuntu), then do this instead
remotes::install_git("git://github.com/tavareshugo/MagicHelpR.git")
```

There's also a basic [introduction vignette](https://tavareshugo.github.io/MagicHelpR//articles/introduction.html).
