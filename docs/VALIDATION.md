# Validation status

The 1.1.1 candidate was tested locally on 25 October 2026. SHA-256 of
`validation/SOURCE_SHA256SUMS`: `d2ea7bc94626bc3b80082c659b4018bc04f7e4d993d8e62aea43603010fd9166`.
The manifest contains 54 code/configuration/example files, sorted by relative path, one
`sha256  path` line per file with LF endings. Documentation, citation/license files and validation
records are excluded to avoid a recursive hash. Run `sha256sum -c validation/SOURCE_SHA256SUMS`
from the repository root to check those exact bytes.

| Local environment | Configure/build | CTest | Install and installed checker |
|---|---|---|---|
| Windows 11, MSVC 19.50, GMP 6.3, Python 3.14 | pass | 10/10, 146.81 s | pass; release GMP DLL hashes match the release dependencies |
| Ubuntu 24.04 container, GCC 13.3, Python 3.12.3 | pass | 10/10, 94.37 s | pass |
| Ubuntu 24.04 container, Clang 18.1.3, Python 3.12.3 | pass | 10/10, 102.39 s | pass |
| Ubuntu 24.04 container, GCC 13.3, ASan + UBSan | pass | 10/10, 125.78 s; no sanitizer reports | pass |

Times are single runs under concurrent load, not performance guarantees. Containers used an existing
dependency image with networking disabled. Windows used the explicit release-library paths in README.

CTest covers all 31 examples, three validity checkers, instance-bound faithfulness, negative and
positive starting-isometry cases, hidden folds, refinements, 13 certificate mutation types,
200 fuzz instances, malformed input, figures (including 20 random cases and metamorphic checks),
and six stored-certificate checks. The new stored figures passed an independent geometry comparison
before their golden hashes were recorded. Both elementary built-in checker examples also pass.
Configuration fails if Python >= 3.11 is unavailable, instead of silently registering no tests.

**Pending release gate:** all five jobs of `.github/workflows/ci.yml` on the exact intended release
commit, including macOS AppleClang. Historical runs on other commits do not satisfy this gate.
No macOS run or hosted CI run on this new history has been claimed. The release remains unpublished.
