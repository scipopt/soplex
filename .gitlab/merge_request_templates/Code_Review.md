### Description

[Add a description of the introduced changes here.]

### Code review

* [ ] The code change is correct.
* [ ] The naming and place of new methods is clear and consistent **or** no new methods have been added.
* [ ] The code is sufficiently documented.
* [ ] The coding style is OK. (Use astyle.)

### Documentation and building

* [ ] The CHANGELOG is up to date (including API changes if present in this MR).
* [ ] The user documentation is up to date (doc/xternal.c, doc/inc/faq/, installation instructions, ...).
* [ ] Both build systems and makedist.sh are up to date. Especially, newly added, renamed or removed source files have been added to, renamed in or removed from src/CMakeLists.txt.

### Testing

* [ ] SCIP ctest has been checked (type some of `jenkins ctest`).
* [ ] SCIP debug runs have been checked (type some of `jenkins debug {short,minlp,mip,pb}`).
* [ ] The performance impact on SCIP has been checked (type some of `jenkins performance {mip,pb} (quick|continue|)`), **or** the changed code will not be executed by default.

### Does this merge introduce an API change? :warning:

* [ ] No, **or** as far as possible, the code ensures backwards compatibility.
* [ ] No, **or** the `SOPLEX_APIVERSION` will be updated. (Run `scripts/updateversion.sh -a` on branch **master** after the changes of this MR have arrived in the master branch.)
* [ ] The changes do not affect SCIP's LPI **or** the LPI has been adjusted.