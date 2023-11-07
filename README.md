## Instructions (to group members)

Always check at least that there are no compilation errors. Simply remove the `LTE-Sim` executable and do `make -j8` and see if the compilation succeeds.

### Git operations

For small changes that do not need review, directly modify the `master` branch. However, always pull before pushing. Do `git fetch origin` and `git merge origin/master` and solve any potential conflicts, then do `git push -u origin`.

For large changes that need review, do it on a different branch.

```bash
git fetch origin
git merge origin/master  # Update to the latest master branch
git checkout -b [feature-branch]  # Create and checkout to a new branch
```

Replace `[feature-branch]` with a customized branch name that is instructive. Make your changes in that branch. Fetch and merge the master branch frequently to keep your branch up-to-date and avoid too much conflicts. When you are done with your modifications, do `git push -u origin [feature-branch]`. Then go to GitHub and make a pull request.

### Style guide

Regardless of what format the simluator previously follows, we will stick to [Google C++ format](https://google.github.io/styleguide/cppguide.html) for any further changes made. One should install `clang-format` and format the code prior to any code commit.

```bash
sudo apt install clang-format  # Install clang formatting tool
git clang-format --style=file  # Do this prior to any commit
```

## Week of 10/30

Added customizable objective (metric) to the RadioSaber scheduler. This can be accessed via `SingleCellWithI` and setting the inter-slice scheduler argument. 9 stands for the original RadioSaber (with spectral efficiency as objective), 91 stands for RadioSaber with proportional fairness, and 92 stands for RadioSaber with M-LWDF (modified largest weighted delay first).

Added a trial experiment script, adapted from RadioSaber but with different interslice metrics.

## Week of 11/06

Prepared configurations for the comparison test. Refactored and implemented the plots for throughput, RB/s, flow completion time, and queuing delay. To run the scripts:

```bash
cd ran-sched-experiments/contradict-objective
./run.sh                                   # Create all logs (test included)
./plot.py --test                           # Plot test log
./plot.py --intra-alg all --flow-type all  # Plot all experiment logs
```
