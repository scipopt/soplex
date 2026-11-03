---
name: Bug report
about: Help us track down bugs in SoPlex

---

Welcome to SoPlex and thank you for reporting a bug!

Please read the following before posting a new bug report:

- If you have a question or are unsure if the behavior you're experiencing is expected, also feel free to reach out on the [mailing list](https://listserv.zib.de/mailman/listinfo/soplex/).

- If you are reasonably confident your issue is a bug in SoPlex, this is the right place. Be sure to include as much relevant information as possible, including a minimal reproducible example. See https://help.github.com/articles/basic-writing-and-formatting-syntax/ for background on how to format text and code on GitHub issues.

- If the issue occurs when running an instance, be sure to include it, along with any non-default setting. See [the documentation on parameters](https://soplex.zib.de/doc/html/PARS.php).
If the issue occurs with a model built programmatically, include a minimal script producing the model.
You can use `writeStateReal("my_model", NULL, NULL, FALSE)` to write the model, the basis, and the settings to files.

- If the problem is SoPlex producing an incorrect solution, be sure to include a reference solution as a sol file, produced either by SoPlex or another solver.
You can also use the `-x=<primfile>` and `-y=<dualfile>` options to write the primal solution and the dual solution to files.

Thanks for contributing to SoPlex!
