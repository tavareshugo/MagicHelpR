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

To install this package, launch R and issue the following commands:

```
#If not already installed, install devtools:
# install.packages("devtools")
devtools::install_github("tavareshugo/MagicHelpR")
```

There's also a basic [introduction vignette](https://tavareshugo.github.io/MagicHelpR/articles/introduction.html).
