BEGIN lp2.spc   for warm starting Quad MINOS
   Rows                     80000
   Columns                 100000
   Elements               1000000

   Old basis file              11
   New basis file              12
   Save Frequency          100000

   Print level                  1
   Print frequency            100

   Scale option                 2
   Iteration limit        1000000
   Expand frequency        100000

   Feasibility tol          1e-15
   Optimality  tol          1e-15

   LU Factor tol             10.0
   LU Update tol             10.0
*  LU               Rook Pivoting
*  LU Singularity tol       1e-15

   Solution                    No
End
BEGIN lp3.spc   for 2nd restart
   Rows                     80000
   Columns                 100000
   Elements               1000000
 
*  Old basis file              12  * Warm start used instead
   New basis file              13
   Save Frequency          100000

   Print level                  1
   Print frequency            100

   Scale option                 0
   Iteration limit        1000000
   Expand frequency        100000

   Feasibility tol          1e-15
   Optimality  tol          1e-15

   LU Factor tol              5.0
   LU Update tol              5.0
*  LU               Rook Pivoting
*  LU Singularity tol       1e-15

   Report file                 81  * Export solution to file fort.81
   Solution                    No
   Solution file                9  * Print solution in E format
End
