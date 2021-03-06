### Description

[Add a description of the introduced changes here.]

### Code review

* [ ] Is the code change correct?
* [ ] Is the code sufficiently documented? Is the coding style OK (use astyle)?
* [ ] Is the naming and place of new methods and parameters clear and consistent?
* [ ] Do emphasis settings need to be adjusted?

### Documentation and building

* [ ] Are CHANGELOG entries added?
* [ ] Is necessary user documentation added (doc/xternal.c, doc/inc/faq/, installation instructions, ...)?
* [ ] Are new files added to makedist.sh and both build systems?  Updated dependencies via makedepend.sh?

### Testing

* [ ] Has ctest been checked: `jenkins ctest`?
* [ ] Has performance als MI(NL)P impact been checked on mi(nl)pdev-solvable (if affecting default)?
* [ ] Are coverage settings added to test new code?
* [ ] Have unit tests been added (if necessary)?

### Does this merge introduce an API change? :warning:

* [ ] Look for a satisfactory solution that ensures backwards compatibility.
* [ ] Document interface changes in the CHANGELOG.
* [ ] Increase SOPLEX_APIVERSION after the merge.
* [ ] Tag this MR with the label 'default parameter' and inform one of the developers responsible for SAP (default: Jakob) if a parameter was added/deleted/changed.
* [ ] if any LPI-changes ocurred, ensure to run the SCIP LPI unittests and the SCIP ctest before merging.