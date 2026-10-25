# Interface and display limitations

These limitations do not enlarge the mathematical input model in SPEC. Use a validity checker and
the evaluator together when asserting construction faithfulness.

| Area | Current behavior and usable workflow |
|---|---|
| Readers | Use standard JSON, ASCII digits in numeric strings, UTF-8 without a BOM, ASCII names, and no repeated closing polygon vertex. C++ accepts some nonstandard raw numbers and controls that Python rejects, and rejects non-ASCII Unicode escapes that Python accepts. Python accepts some Unicode numeric strings that C++ rejects. |
| Failure JSON | A C++ check failure mentioning a control character in a cell name can contain invalid JSON. Avoid controls in names; treat a nonzero exit as a failure even if its diagnostic cannot be parsed. |
| Resource status | Parse-time C++ allocation failures are reported as invalid input (exit 3), rather than out of memory (5). Other allocation failures and OS termination have distinct paths. |
| Order commands | Use plain nonnegative decimal iteration counts and seeds and the documented argument count. Numeric suffixes, negative iterations and extra arguments are not reliably rejected. Malformed instance input can return 4 instead of 3. Only construct validates point dimensions strictly. |
| Empty order | Omit `order` or supply a permutation. An explicit empty list defaults in construct but is rejected by the evaluator. |
| Binding | C++ check verifies the certificate's own geometry. The evaluator compares vertices and source images and ignores coverage, continuity and domain binding. Use `ref_check.py --instance` as well. Empty or incomplete cell lists cannot be certified by evaluator agreement alone. |
| Paths | Keep input and output paths distinct. Construct can overwrite its input. The figure tool writes through `OUT.svg.part`; neither that path nor OUT.svg should be an input or existing directory. An output failure may leave the part file. Use an explicit trusted `--exe` on Windows and paths representable in its current ANSI code page. |
| Windows install | The default vcpkg library search can choose debug GMP DLLs even for Release. Use the explicit release-library paths in README, or the release build-tree executable. The default installed DLLs may need a debug CRT. |
| Diagnostics | Outside-source errors do not name the source; too-small-L errors do not give a numerical bound. Inspect SPEC's exact inequalities or omit L to use the sufficient default. Missing-file messages in Python can omit the path. |
| Options | Unknown instance keys are ignored. Matrices use rows. Inspect the certificate's construction record to confirm the resolved choices. |
| Display | Coordinates are rounded to 0.001 px and range labels are approximate. Tiny cells and close labels can disappear or overlap; all cells remain in the SVG. Zoom and exact certificate inspection are necessary for small features. A checked figure verifies R1–R6, not its displayed construction record or any physical layer order. |
| Evaluation | There is no general point-evaluation CLI. Apply the listed affine map on a containing cell. A validity check guarantees agreement on shared boundaries. The faithfulness evaluator can return undecided. |

No formal verification or independent human mathematical review is claimed. This is exact rational
planar extension software; it does not optimize subdivisions or decide physical foldability.
