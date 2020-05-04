
0.2.0
=====

-   Importance parameters *ψ* are no longer calculated when non-parametric response functions are used.

0.3.0
=====

-   Importance parameters *ψ* are enabled again when non-parametric response functions are used, but not used for plotting.
-   2D sample plots for constrained ordination with non-parametric response functions have been disabled, as they are not interpretable. Variable plots are the only 2D plots still allowed
-   Explained deviance and inertia can be plotted on the axes rahter than the *ψ*'s using the "plotPsi" argument to the plot.RCM() function.
-   Possibility to provide lower dimensional fits has been disabled. *RCM* is fast enough to fit the whole model.

1.0.0
=====

-   Release on BioConductor

1.0.1
=====

-   Bug fix in buildCovMat() to avoid false warning
-   Check for alias structure in confounder and covariate matrices

1.2.0
=====

-   Missing values in count matrix are now allowed. They simply do not contribute to the parameter estimation, but the rest of the row (or column) is still used.

1.2.1
=====

-   Vertical reference line in residual plot
-   Bug fix for problematic variable names

1.2.2
=====

-   Moving the online manual information to the vignette

1.2.3
=====

-   Rename *a* and *b* to *rowExp* and *colExp* to avoid partial matching
-   Allow *rowExp* and *colExp* to be adapted for constrained correspondence analysis starting values as well

1.2.4
=====

-   Adding a new inflVar variable to disambiguate in the influence plotting
-   More argument checking + tests for the plot.RCM function
