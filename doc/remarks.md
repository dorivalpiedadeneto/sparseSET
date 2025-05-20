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
in normal cases.
Merge sort: time efficiency -> O(n log_2(n)) -> aprox. 5537 operations
Quick sort: time efficiency -> O(n log_10(n)) -> aprox. 1667 operations

