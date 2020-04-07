BEGIN FBA.spc   for cold start on tough FBA problems
   Rows                     80000
   Columns                 100000
   Elements               1000000
 
*  Old basis file              11
   New basis file              11
   Save Frequency            1000 

   Print level                  1
   Print frequency            100

   Scale option                 2
   Iteration limit        2000000
   Expand frequency        100000

   Feasibility tol           1e-6
   Optimality  tol           1e-6

   LU               Rook Pivoting
   LU Factor tol              4.0
   LU Update tol              4.0

*  Report file                 81  * Export solution to file fort.81
   Solution                    No
   Solution file                9  * Print solution in E format
End
