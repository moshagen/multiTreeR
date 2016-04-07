# define eqn and mdt file
fEqn <- system.file("extdata", "source.eqn", package="multiTreeR")
fMdt <- system.file("extdata", "source.mdt", package="multiTreeR")

# fit model using restrictions a=g and D1=D2
res <- doMT(fEqn, fMdt, restrictions = list('a=g','D1=D2'))

# show results
summary(res)
