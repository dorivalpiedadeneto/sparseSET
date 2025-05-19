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

