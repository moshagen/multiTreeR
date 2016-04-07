multiTreeR
==========

multiTreeR is an R interface for multiTree. For more details, see <http://psycho3.uni-mannheim.de/multiTree>

Installation
------------

To install multiTreeR, run:

``` r
install.packages("devtools")
library(devtools)
install_github("moshagen/multiTreeR")
```

The Java Runtime Environment has to be installed on the target machine, get it here: <http://www.java.com>

Using multiTreeR
----------------

1.  Define the path to an existing MPT model file in .eqn format: `eqn <- 'source.eqn'`
2.  Either provide a data frame or define the path to an mdt: `mdt <- 'source.mdt'`
3.  Either define the path to an restrictions file or provide a list set parameter-restrictions: `restr <- list('a=g','D1=D2')`
4.  To run the analysis, call `res <- doMT(eqn, mdt, restrictions=restr)`
5.  Inspect results by calling `summary(res)` or by peeking in the results object `res$paramEst`

### References

-   Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. *Behavior Research Methods, 42,* 42–54. <http://doi.org/10.3758/BRM.42.1.42>
