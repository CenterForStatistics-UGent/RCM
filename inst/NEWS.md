
# 0.2.0

  - Importance parameters \(\psi\) are no longer calculated when
    non-parametric response functions are used.

# 0.3.0

  - Importance parameters \(\psi\) are enabled again when non-parametric
    response functions are used, but not used for plotting.
  - 2D sample plots for constrained ordination with non-parametric
    response functions have been disabled, as they are not
    interpretable. Variable plots are the only 2D plots still allowed
  - Explained deviance and inertia can be plotted on the axes rahter
    than the \(\psi\)’s using the “plotPsi” argument to the plot.RCM()
    function.
  - Possibility to provide lower dimensional fits has been disabled.
    *RCM* is fast enough to fit the whole model.

# 1.0.0

  - Release on BioConductor

# 1.0.1

  - Bug fix in buildCovMat() to avoid false warning
  - Check for alias structure in confounder and covariate matrices

# 1.2.0

  - Missing values in count matrix are now allowed. They simply do not
    contribute to the parameter estimation, but the rest of the row (or
    column) is still used.

# 1.2.1

  - Vertical reference line in residual plot
  - Bug fix for problematic variable names

# 1.2.2

  - Moving the online manual information to the vignette

# 1.2.3

  - Rename *a* and *b* to *rowExp* and *colExp* to avoid partial
    matching
  - Allow *rowExp* and *colExp* to be adapted for constrained
    correspondence analysis starting values as well

# 1.2.4

  - Adding a new inflVar variable to disambiguate in the influence
    plotting
  - More argument checking + tests for the plot.RCM function

# 1.5.2

  - Avoid returning nulls for residualPlot

# 1.5.4

  - A note in the vignette and in the help file of plot.RCM regarding
    limited number of combinations of constraining variables. Also a
    warning is now thrown

# 1.5.5

  - Bug fix for higher dimension residualPlot function, and tests for
    this function

# 1.5.6

  - Replace deprecated guides( =FALSE) by guides(=“none”)

# 1.5.7

  - Update vignette to number table of contents

# 1.11.0

  - Added FAQ section in vignette with first frequent question on number
    of samples not shown.
  - Fixed bugs for plots of data with missing values, and added tests.

# 1.11.2

  - Explicitly import stats::model.matrix, and only load necessary VGAM
    functions

# 1.11.3

  - For the unconstrained models: fit feature models one by one and
    Gram-Schmidt orthogonalize and center afterwards, rather than using
    Lagrange multipliers and huge Jacobian matrices. This will use less
    memory and speed up computations, but *may yield slightly different
    solutions*. Nothing changes for the constrained models.
    
# 1.11.4

 - Introduction of permanova testing for user-supplied groups using the _permanova_ function.
