# Testsuite

In general the execution of the testruns is the following:

- To start a testrun call `make test` for a local one, `make testcluster` for one on the cluster.
  (These are for running soplex tests but variants exist for qsoptex and perplex)
- The `Makefile` will call the corresponding `check.sh` or `check_cluster_<solver>.sh` script
  with the correct variables. These scripts will
  + configure the environment variables for local and cluster runs `configuration_set.sh`, `configuration_cluster.sh`
  + configure the test output files such as the .eval, the .tmp and the .set files.
  + cluster run scripts have a `waitcluster.sh` script in the loop that makes the scripts wait
    instead of overloading the cluster queue with jobs.
    + run `test.sh` on the instances and call `evaluation.py` for local tests or
    + `runcluster_<solver>.sh` to execute the tests with the solver, later `evalcheck_cluster.sh` has to be called to evaluate the logfiles which in turn calls `evaluation.py`.
  + `evaluation.py` parses the logfiles into the result table and checks the solutions recorded in them.

## local make targets

### solve a number of standard settings on a short testset

make check
  - `check.sh`
    + `configuration_set.sh`
    + `test.sh`
      - `evaluation.py`

### solve a testset locally

make test
  - `test.sh`
    + `configuration_set.sh`
    + `evaluation.py`

## cluster make targets

### solve a testset on the cluster with soplex

make testcluster
  - `check_cluster.sh`
    + `configuration_cluster.sh`
    + `configuration_set.sh`
    + `runcluster.sh`

### solve a testset on the cluster with perplex

make testclusterperplex
  - `check_cluster_perplex.sh`
    + `configuration_cluster.sh`
    + `configuration_set.sh`
    + `runcluster_perplex.sh`

### solve a testset on the cluster with qsoptex

make testclusterqsoptex
  - `check_cluster_qsoptex.sh`
    + `configuration_cluster.sh`
    + `configuration_set.sh`
    + `runcluster_qsoptex.sh`

## other Scripts

### Compare several testruns on their detailed per-instance fields.

- `compare.py`
  Compares several runs on any configurable per-instance fields (default `solvetime` and `iters`),
  writes performance profiles (`.prf`) and shifted geometric means, and reads the `.json` files that
  `evaluation.py` writes next to each `.res`, e.g. `compare.py results/check.*.json`.

# Files

## AWK files

- `make_solu_file.awk`
  Takes two files: <CPLEX log> <Perplex log> (order is critical!) and generates a .solu file.

## Bash Scripts

- `check.sh`
  Solves a number of standard settings with a short testset
  Called by `make check` from soplex root

- `test.sh`
  solve a given testset with given settings and time limit parameters:
    1: name of testset (has to be in check/testset)
    2: path to soplex executable
    3: name of settings (has to be in settings)
    4: time limit
    5: results directory

- `check_cluster*.sh`
  Call with "make testcluster", "make testclusterperplex" and "make testclusterqsoptex"
  The queue is passed via $QUEUE (possibly defined in a local makefile in soplex/make/local).
  For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
  with valid time measurements should be executed, this number should be chosen in such a way
  that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
  of each computer. Of course, the value depends on the specific computer/queue.
  With `EXCLUSIVE=auto` multiple instances are packed onto one exclusive node, each bound to its own
  physical core, and the number of instances per node is computed automatically. Combined with
  `DEBUGTOOL=perf` the solving time is frequency-normalized for more reproducible time measurements.

### Configuration

- `configuration_cluster.sh`
  Configures environment variables for cluster runs (queue, constraint, account, cpu and target frequency).
  It is to be invoked inside a `check_cluster*.sh` script.

- `configuration_set.sh`
  Configures environment variables that are needed for test runs both on the cluster and locally.
  It is to be invoked inside a `check(_cluster)*.sh` script.

### Running

- `runcluster*.sh`
  The script executing the solver on one instance and producing the logfiles.
  Can be executed either locally or on a cluster node.
  Is to be invoked inside a `check_cluster*.sh` script.

### Evaluation

- `evalcheck_cluster.sh`
  evaluate logfiles from testrun started by `make testcluster`. It injects the frequency-normalized
  time from any `.perf` files and calls `evaluation.py`. With `-U` the results are uploaded to
  rubberband.

- `evaluation.py`
  Parses a test run `.out` file into the `.res` result table, checks the recorded solutions, and writes a
  per-instance `.json` with all parsed fields next to the `.res`.
  Called by `test.sh` and `evalcheck_cluster.sh`.

### helper

- `waitcluster.sh`
  In order to not overload the cluster, no jobs are submitted if the queue is too full
  instead, this script waits until the queue load falls under a threshold and returns
  for the calling script to continue submitting jobs.

- `make_solu_file.sh`
  Shell script for generating .solu file from a .test file using CPLEX and perplex, calls
  + `make_solu_file.awk`

## folders and other files

- `CMakeLists.txt` is part of the cmake system, this file's main purpose is to add tests.

- In the `testset/` directory reside the *.test and *.solu files to be specified via `TEST=short`.
  A `.test` file lists problem files, one file per line, absolute and relative paths to its location.
  A `.solu` file with the same basename as the `.test` file contains information about feasibility and best known objective value.
  It is optional for a testrun.

- In the `instances` direcotry are some example instances that can be solved with SoPlex.
