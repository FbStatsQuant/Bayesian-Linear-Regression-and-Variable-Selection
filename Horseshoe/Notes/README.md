# Horseshoe regression notes

- `horseshoe_linear_regression.tex`: editable, self-contained LaTeX source.
- `horseshoe_linear_regression.pdf`: compiled review draft.

The notes derive Gaussian horseshoe regression with a separate unpenalized
intercept, a noise-scaled coefficient prior, and a proper inverse-gamma noise
prior. They include all full conditionals, an ordered Gibbs sweep, numerical
linear algebra, prior calibration, posterior prediction, variable-selection
decisions, diagnostics, and alternative parameterizations. References are
embedded in the LaTeX source; no separate bibliography processor is needed.

## Rebuild in PowerShell

Run from this directory with a LaTeX distribution on `PATH`:

```powershell
New-Item -ItemType Directory -Force -Path '.build' | Out-Null
for ($pass = 1; $pass -le 3; $pass++) {
    pdflatex '-interaction=nonstopmode' '-halt-on-error' '-file-line-error' '-output-directory=.build' 'horseshoe_linear_regression.tex'
    if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed on pass $pass" }
}
Copy-Item -LiteralPath '.build/horseshoe_linear_regression.pdf' -Destination 'horseshoe_linear_regression.pdf'
```

Three passes resolve the table of contents, citations, page numbers, and cross
references. Inspect `.build/horseshoe_linear_regression.log` for warnings and
render the resulting PDF to check layout after revisions. Build intermediates
and local visual-review files belong in `.build/`, which is ignored by Git.

This is an explanatory document, with analytic examples only. No sampler
simulation, benchmark, or numerical validation has been run. Future
computational experiments follow the project's Owls execution and audit rules.
