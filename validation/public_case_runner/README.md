# Public Case Runner Regression

This workflow verifies the public case-runner interface rather than introducing a new
physical benchmark. The unified validation command:

1. runs both example cases from `/tmp`, outside the repository working
   directory;
2. checks the expected portable result files and structured metadata;
3. compares the generic quadratic Burton-Miller result at `ka=9.3558` with the
   fictitious-frequency regression node by node.

The comparison must have relative L2 difference no larger than `1e-12`. This
detects errors in TOML parsing, frequency derivation, boundary/incident-field
mapping, normalized mesh writing, namelist generation, or generic Fortran case
configuration. Analytical accuracy remains the responsibility of the dedicated
analytical validations.
