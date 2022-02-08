************************************************************************
file with basedata            : mm2_.bas
initial value random generator: 1424959589
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  7
horizon                       :  40
RESOURCES
  - renewable                 :  1   R
  - nonrenewable              :  0   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     5      na      na       na       na
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          1           2
   2        2          4           3   4   5   6
   3        2          1           7
   4        1          1           7
   5        1          7           7
   6        1          1           7
   7        1          0
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1
------------------------------------------------------------------------
  1      1     0       0
  2      1     5       4
         2     4       4
  3      1    11       1
         2    12       2
  4      1     5       2
  5      1     4       1
  6      1     1       3
  7      1     0       0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1
    4
************************************************************************