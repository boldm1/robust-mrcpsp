************************************************************************
file with basedata            : md378_.bas
initial value random generator: 357464772
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  22
horizon                       :  163
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     20      0       29        4       29
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6   7
   3        3          3           5   6   8
   4        3          2          10  14
   5        3          1          19
   6        3          3          10  11  20
   7        3          3           8  13  20
   8        3          2           9  12
   9        3          3          16  17  19
  10        3          2          15  17
  11        3          3          13  15  17
  12        3          2          14  21
  13        3          2          14  21
  14        3          1          19
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        3          1          21
  19        3          1          22
  20        3          1          22
  21        3          1          22
  22        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       0    3    2    3
         2     9       0    3    2    2
         3    10       4    0    1    1
  3      1     1       9    0    6    3
         2     7       9    0    4    2
         3    10       0    3    3    2
  4      1     4       8    0    7    2
         2     6       0    7    5    2
         3     9       3    0    4    2
  5      1     1       8    0    6   10
         2     2       0    6    6   10
         3     7       0    5    4   10
  6      1     4       0   10    3    3
         2     5       3    0    2    3
         3     6       0    7    1    3
  7      1     3       8    0    9    6
         2     8       7    0    8    6
         3    10       0    9    8    3
  8      1     5       0    8    7    9
         2     6       0    7    5    9
         3     6       0    8    2    7
  9      1     5       0    7    8    7
         2     9       6    0    8    4
         3     9       0    5    7    2
 10      1     1       0    7    9    7
         2     7       0    5    8    6
         3    10       7    0    8    3
 11      1     6       3    0    9    5
         2     8       0    8    5    2
         3     8       1    0    5    2
 12      1     1       0    9    2    8
         2     1       4    0    2    6
         3     7       0   10    2    6
 13      1     5       0    1    8    3
         2     8       0    1    4    3
         3    10       3    0    1    2
 14      1     4       0   10    6    5
         2     5       0    9    5    3
         3     9       0    9    5    1
 15      1     4       0    3    2   10
         2     6       0    2    1   10
         3     6       0    1    2    9
 16      1     1      10    0    8    6
         2     6       0    4    7    4
         3     7       9    0    7    3
 17      1     4       0    8    6    9
         2     6       0    8    2    8
         3     6       0    7    4    8
 18      1     4       0    5    6    9
         2     6       0    4    5    8
         3     6       0    5    5    6
 19      1     7       0    7    6    6
         2     9       0    5    4    5
         3    10       0    4    4    5
 20      1     5       0    5    4    8
         2     5       0    6    3    8
         3     7       4    0    2    7
 21      1     3       5    0   10    5
         2     4       4    0    7    5
         3    10       3    0    7    4
 22      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   14   20  124  124
************************************************************************
