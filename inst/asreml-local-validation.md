# Local ASReml-R validation

PredictProR treats ASReml-R as an optional proprietary backend. The package
must install, load, and run its non-ASReml models when ASReml-R or an ASReml
license is unavailable.

## License configuration

Do not put a license-server address or license file in the package source. Set
the vendor-supported `vsni_LICENSE` environment variable outside the project:

```text
vsni_LICENSE=<port>@<host>
```

For the current shell, use `$env:vsni_LICENSE = "<port>@<host>"` in PowerShell
or `export vsni_LICENSE='<port>@<host>'` in a POSIX shell on macOS, Ubuntu, or
another Linux distribution. A user-level `.Renviron` entry is also portable
across R sessions. Protect that file as local configuration and do not commit
it.

The installed ASReml-R binary must match the operating system and the R
major/minor version. It may be installed in a user library or in the ignored
project library `.r-lib`.

## Fail-loud reference checks

Run from the package root:

```sh
Rscript tools/validate_asreml_reference.R
Rscript tools/validate_asreml_reference.R --full
```

The quick check validates license checkout, delta-method heritability standard
errors, CORGH reconstruction, and a live two-kernel FA1 MET fit. The full check
also exercises multi-trait and hybrid ASReml-R routes. Neither command prints
the configured server address.

These checks require network reachability to the configured license service.
They are intentionally opt-in so standard package tests remain portable and
license-independent.
