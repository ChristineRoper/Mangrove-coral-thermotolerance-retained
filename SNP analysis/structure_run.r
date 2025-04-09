library(dartR)

load("gl_for_structure.RData")

sr <- gl.run.structure(
  gl_for_structure,
  burnin = 10000,
  numreps = 10000,
  save2tmp = TRUE,

  # Documentation of options here: https://rdrr.io/cran/strataG/man/structure.html#heading-2
  k.range = 1:10,
  num.k.rep = 5, # number of replicates for each value in k.range.

  # The location of the external STRUCUTRE program to run
  exec = "./STRUCTURE/structure"
)

save(sr, file = "structure_run.RData")
