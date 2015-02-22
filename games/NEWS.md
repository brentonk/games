# games 1.1.2

* Fixed bugs in `egame12()` and `egame122()` that caused them to fail when the left-hand side was constructed with binary variables of class `"integer"` (thanks to Tyson Chatagnier for the bug report)

* Streamlined `Mode()` function according to the Stack Overflow answer <http://stackoverflow.com/a/8189441/143383>

* Reduced computation time of `egame122()` and `leblang2003` examples


# games 1.1-1

* Includes citation to published article in *Journal of Statistical Software*


# games 1.1-0

* `ultimatum()` has new argument `minOffer`

* Fixed bug where `predict.ultimatum()` computed E[offer] incorrectly; now using simulation due to absence of closed-form expression


# games 1.0-5

(thanks to Journal of Statistical Software reviewers for suggestions)

* All model classes gain new element `localID`, a logical indicator for whether the fitted model is locally identified (i.e., negative definite Hessian).  A warning is issued during fitting, and in `print` and `summary` methods, when a model is not locally identified.

* `profile.game()` gains argument `use.se`, allowing `dist` to be specified absolutely instead of as the number of standard errors (useful for models where local identification fails and estimated standard errors are ill-defined)


# games 1.0-4

* Fixed bug in the computation of the Vuong test statistic in `vuong()`


# games 1.0-3

(thanks to Julian Wucherpfennig for the suggestion)

* `predProbs()` has new argument `type` to select whether to generate probabilities for outcomes (the default) or actions

* The name of the corresponding argument in the `predict` methods has been changed from `probs` to `type` (for consistency with `predict.glm()` and others)


# games 1.0-2

* Fixed typo in `egame12()` documentation


# games 1.0-1

* Added arguments `caption` and `label` to `latexTable()`


# games 1.0-0

## Major additions

* New model `egame123()`, a three-player extensive form game (and associated simulated dataset `data_123`)

* New functions `clarke()` and `vuong()` for non-nested model testing

* New interactive function `makeFormulas()` for constructing model formulas

* Added dataset `leblang2003`, for replication of David Leblang's 2003 article "To Devalue or to Defend? The Political Economy of Exchange Rate Policy"

* All fitting functions now include argument `method` for selection of optimization algorithm from any of the options available in the `maxLik` package

* New dependencies: package `stringr`

## Minor additions and changes

* Added option `report` to `profile.game()` for printing status bar while computing profile log-likelihoods

* Fitted model objects' `convergence` component now contains information on the optimization routine used and the number of iterations to convergence

* Fitted `"ultimatum"` objects include both offer and acceptance (if applicable) in `y` component as matrix, rather than separately

* Log-likelihood component of fitted `"ultimatum"` objects where `outcome == "both"` includes attributes with individual log-likelihood contributions from offers and from acceptance

* Fitted models include component `xlevels` to store which levels of each factor variable were used in fitting (similar to the attribute of the same name in `"lm"` and `"glm"` objects)

* For consistency with previous software and published results, utilities fixed to 0 (e.g., u23 in `egame12()`) are now assumed to have a stochastic component in the private-error case

* For consistency, the variable name for a utility equation containing only a constant is now displayed in the form `"u1(1):(Intercept)"`, rather than being truncated to `"u1(1)"`

## Bug fixes

* `predict.game()` works even if `newdata` does not contain variables corresponding to the left-hand side of the original fit

* `predProbs()` uses the modal value of binary variables by default when constructing a regressor profile, rather than the first value that appears in the dataset

* `profile.game()` no longer refits starting at originally estimated coefficients (this would yield infinitesimal improvements in the log-likelihood from nearly identical coefficients, causing the non-convergence warning to be issued inappropriately)

* `plot.profile.game()` uses smarter defaults for y-axis limits

* In fitting functions, `formulas` specified as a list will no longer fail if the list has names

* `latexTable()` output now places the intercept in the first row


# games 0.6-0

## New features

* Introduced functions `profile.game()` and `plot.profile.game()` for diagnosing failure to converge to a global maximum

* All fitting functions (`egame12()`, `egame122()`, `ultimatum()`) now include an argument `profile`, for using profile output from a previous fit attempt to move closer to the global maximum

* New dataset `student_offers` for illustrating convergence failure in the ultimatum model

* Added indicator for whether gradient was used in fitting to the convergence element of `"game"` objects

* Fitted `"ultimatum"` objects now include the outcome of interest (offer or both), the acceptance vector (if any), and the offer tolerance

## Bug fixes

* Fixed error with starting value generation in ultimatum model when acceptance is not specified


# games 0.5-1

* Fixed error with starting value generation in ultimatum model when offer variable is of class `"integer"`

* Fixed numerical issues in ultimatum model


# games 0.5-0

* Initial release of package
