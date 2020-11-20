# Branch-and-bound tree visualization

In scip, set visual/vbcfilename=vcbfilename and visual/dispsols=TRUE. After you have solved a problem, you can do

`ruby vbc2dot.rb vbcfilename [options]`

This will generate you a visualization of the branch-and-bound-tree of the problem.
