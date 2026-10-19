## Submission

This is the first submission of 'scramble', which implements the method of
Pfeiffer, Vana-Gür and Filzmoser (2025) <doi:10.1007/s11634-025-00656-3>.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local macOS 26.6 (aarch64), R 4.6.1, `R CMD check --as-cran`
  - with the 'torch' backend installed, running all examples
    (incl. `\donttest{}`), tests and the vignette
  - without the 'torch' backend, as on the CRAN check machines

## Use of 'torch'

The package builds on 'torch'. Its computational backend (libtorch) is
downloaded by the user after installation (`torch::install_torch()`) and is
not available on the CRAN check machines. Therefore:

* examples are wrapped in `if (torch::torch_is_installed())`,
* tests call `skip_if_not(torch::torch_is_installed())`,
* the vignette is evaluated only when the backend is available.

With the backend installed, all examples, tests and the vignette run without
errors or warnings.

The example that tunes the penalty parameter by Bayesian optimization takes
about 40 seconds and is therefore wrapped in `\donttest{}`.

## Reverse dependencies

None (new package).
