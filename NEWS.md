# HT12ProcessoR 0.0.18

* Added a `NEWS.md` file to track changes to the package.
* New parameter round4ANOVAcheck for control of function checkBatchEffects() -->
  calcAnovaSentrixRunSpecialbatchViaMatrixEQTL2() -->
  runMAtrixEQTLAnova(round4ANOVAcheck = 5) to allow slightly different results when e.g. using completely nested covariate structure


# HT12ProcessoR 0.0.19
 * fix runMAtrixEQTLAnova()

# HT12ProcessoR 0.0.20
 * fix writeFilesTosend() to allow when a subgroup is always sampleinfo$in_study = F

# HT12ProcessoR 0.0.21
 * making fix writeFilesTosend() save for data.table and data.frame
