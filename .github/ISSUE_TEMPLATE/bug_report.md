---
name: Bug report
about: Help us track down bugs in SOPLEX

---

Welcome to SOPLEX and thank you for reporting a bug!

Please read the following before posting a new bug report:

- If you have a question or are unsure if the behavior you're experiencing is expected, also feel free to reach out on the [mailing list](http://listserv.zib.de/mailman/listinfo/scip/).

- If you are reasonably confident your issue is a bug in SOPLEX, this is the right place. Be sure to include as much relevant information as possible, including a minimal reproducible example. See https://help.github.com/articles/basic-writing-and-formatting-syntax/ for background on how to format text and code on GitHub issues.

- If the issue occurs when running an instance, be sure to include it, along with any non-default setting. See [the documentation on parameters](https://www.scipopt.org/doc/html/PARAMETERS.php)
If the issue occurs with a model built programmatically, include a minimal script producing the model.
You can use `SCIPwriteOrigProblem(scip, "my_model.lp", NULL, FALSE)` to write the model to a file.

- If the problem is SOPLEX producing an incorrect solution, be sure to include a reference solution as a sol file, produced either by SOPLEX or another solver.
You can also use `SCIPprintSol(scip, SCIPgetBestSol(scip), file, FALSE)` to print a solution to a file.

Thanks for contributing to SOPLEX!
