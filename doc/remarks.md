# Remarks about the project

## 2025/05/13

Until now there was a variable named sym (boolean) to hodl the information
about the sparse matrix symmetry. It seems that a better approach is not
to deal with symmetry this way. In some situations the library user wants
to deal with all the matrices terms, even if the matrix is symmetric; in
other, he/she rathers dealing only with the upper/lower terms. So, it seems
that the best approach is to have a variable to hold one of the three options:
'full', 'upper', 'lower'. Maybe a good name for this variable is 'storage'.
This, together with mtype variable may be a good way to deal with the 
sparse\_matrix storage logic. 

## 2025/05/19

Probably the best algorithm to be used for sorting will be mergesort
https://stackoverflow.com/questions/53269004/heap-sort-vs-merge-sort-in-speed
(O(n.logn), but probably faster due to caching)

## 2025/05/20

After reading more about the sorting methods, probably the best choice
is the quicksort algorithm. A good implementation is available in Numerical
recipes in Fortran 90 (it is necessary to adapt some things, but this
implementation also deal with the issue of the worst case efficiency
(according to the volume I versio, Numerical Recipes in Fortran 77).

Considering n=600 (probably one of the largests cases for a sparse line
in normal cases.)
Merge sort: time efficiency -> O(n log_2(n)) -> aprox. 5537 operations
Quick sort: time efficiency -> O(n log_10(n)) -> aprox. 1667 operations

## 2025/05/21

The sorting subroutine will be coded to return an array of integers
representing the position the index should be moved such that they would be
in ascending order. For doing this, an array with integers, ranging from
1 to n (n being the number of indexes in the array) is necessary.

Considering this array from 1 to n will be necessary several times (each
time the sorting is necessary), we tested how much time would be necessary
for creating this '1 to n' array (comparing the time necessary to copy it
from a global variable previously created, for instace).

A test for evaluating this time difference was coded in the ptest.f90 file.
The results indicate that it does not make a significant difference between
creating it from scratch each time or copying it from the global array.

The results we found:

$ make run-ptest
gfortran -O2 src/test/ptests.f90 -o bin/ptests
For an index array of size  100
 - Time to create index array (1000 times):    0.00001000 seconds
 - Time to copy index array (1000 times):      0.00002300 seconds
For an index array of size  200
 - Time to create index array (1000 times):    0.00002700 seconds
 - Time to copy index array (1000 times):      0.00000900 seconds
For an index array of size  300
 - Time to create index array (1000 times):    0.00005300 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds
For an index array of size  400
 - Time to create index array (1000 times):    0.00005200 seconds
 - Time to copy index array (1000 times):      0.00002000 seconds
For an index array of size  500
 - Time to create index array (1000 times):    0.00003600 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds
For an index array of size  600
 - Time to create index array (1000 times):    0.00002000 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds
For an index array of size  700
 - Time to create index array (1000 times):    0.00002000 seconds
 - Time to copy index array (1000 times):      0.00002000 seconds
For an index array of size  800
 - Time to create index array (1000 times):    0.00001900 seconds
 - Time to copy index array (1000 times):      0.00002000 seconds
For an index array of size  900
 - Time to create index array (1000 times):    0.00001900 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds
For an index array of size 1000
 - Time to create index array (1000 times):    0.00002000 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds

Even though for small matrices there is some difference, it is not
signficant (at least for this computer in which we performed the
tests). (Notice that for no optimization, i.e., using the -O0
option in gfortran, there is a greater difference. However, if one
uses no optimization, probably performance is not something in
such a situation).


