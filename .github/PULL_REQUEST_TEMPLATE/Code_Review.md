### Description

[Add a description of the introduced changes here.]

### Code review

* [ ] The code change is correct.
* [ ] The naming and place of new methods is clear and consistent **or** no new methods have been added.
* [ ] Ensure that new code and features use randomization where necessary to account for performance variability.
* [ ] The code is sufficiently documented.
* [ ] The coding style is OK, see `astylecfg`.

### Documentation and building

* [ ] The CHANGELOG is up to date (including API changes if present).
* [ ] The user documentation is up to date (doc/xternal.cpp, doc/inc/, installation instructions, ...).
* [ ] Both build systems and scripts/makedist.sh are up to date. Especially, newly added, renamed or removed source files have been added to, renamed in or removed from src/CMakeLists.txt.

### Testing

* [ ] ctest passes without errors.
* [ ] The performance impact has been checked **or** the changed code will not be executed by default.
* [ ] The new code is sufficiently covered by tests (perhaps, new coverage settings or new unit tests have been added).

### Does this pull request introduce an API change? :warning:

* [ ] No, **or** as far as possible, the code ensures backwards compatibility.
* [ ] No, **or** the `SOPLEX_APIVERSION` is updated.
