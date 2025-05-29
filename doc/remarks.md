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

```
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

For an index array of size 1000*
 - Time to create index array (1000 times):    0.00002000 seconds
 - Time to copy index array (1000 times):      0.00001900 seconds
```

Even though for small matrices there is some difference, it is not
signficant (at least for this computer in which we performed the
tests). (Notice that for no optimization, i.e., using the -O0
option in gfortran, there is a greater difference. However, if one
uses no optimization, probably performance is not something relevant
in such a situation).


## About the sorting algorithm

Due to an error in the adapted sorting routine, we tested another
implementation (a recursive one). Now that the bug in the Numerical
Recipes based function is correctly working, it is possible to
compare their performance.

The results for a 600 length line.

```
 Recursive version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort         600  terms:   1.8500000000000005E-004  (s)

 Numerical Recipes version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort         600  terms:   7.9999999999999776E-006  (s)

```

The value 600 is the length of a sparse\_line before assembling it.
(we believe this is maximum length value for a common case of meshes
of high order solid element modelling; generally 100 ou 200 would be
typical values)

In this case, considering the time to sort one million lines, for
the recursive subroutine would be 185 seconds; the time to sort using
the Numerical Recipes based function would be 8 seconds.

Testing other values (only to undestand how perfomance changes)

For length 100:

```
Recursive version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort         100  terms:   4.9999999999999697E-005  (s)

 Numerical Recipes version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort         100  terms:   7.0000000000000617E-006  (s)
 ```

For length 1000:

```
Recursive version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort        1000  terms:   5.9999999999999984E-004  (s)

 Numerical Recipes version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort        1000  terms:   1.0999999999999725E-005  (s)

```

For length 10000:

```
Recursive version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort       10000  terms:   3.4148999999999999E-002  (s)

 Numerical Recipes version
 Testing quicksort: worked? (T/F?) -> T
 Time to sort       10000  terms:   1.5300000000000036E-004  (s)
```

**In summary:**

**Trs**: Time using recursive solution

**Tnr**: Time using Numerical Recipes solution

|  length    |    Trs (s)   |    Tnr (s)   |   Trs/Tnr | 
|------------|--------------|--------------|-----------|
|    100     |    5.00E-5   |    8.00E-6   |   6.25    |
|    600     |    1.85E-4   |    7.00E-6   |  26.42    |
|   1000     |    6.00E-4   |    1.10E-5   |  54.54    |
|  10000     |    3.41E-2   |    1.53E-4   | 222.88    |

## 2025/05/29

### Some performance test to evaluate the current implementation.

The sorting function returns a integer vector with the order of the
sorted values. Them, the sparse lines indexes and values are copied,
sorted using this vector, summed and copied back to the sparse line.

All this copies may be something bad for the performance.

To evaluate how much the whole process (sorting, copying and summing)
can be bad for the performance, let's compare the times for just
sorting and for the process:

_(using the subroutines  from ptest.f90) _

```
gfortran -O2 src/test/ptests.f90 -o bin/ptests
  25k lines, terms/lines: 100
  -> Time to assemble (s):  0.22111700000000001
  100k lines, terms/lines: 200
  -> Time to assemble (s):   2.7778140000000002
  100k lines, terms/lines: 600
  -> Time to assemble (s):   20.100937000000002
  25k lines, terms/lines: 100 (only sorting)
  -> Time to sort (s):   5.6910000000002015E-002
  100k lines, terms/lines: 200 (only sorting)
  -> Time to sort (s):  0.49897200000000197
  100k lines, terms/lines: 600 (only sorting)
  -> Time to sort (s):   1.8209410000000013
```

In fact, the adopted approach presents a bad performance.


A new method to order both indexes and values at the same time was 
created in ptest. It was implemented as a subroutine that only
swaps terms (instead of creating a new integer vector and returning
it ordered as a function). The time difference is huge.
First it seemed that the new method is wrong (due to such a great
time difference). The first tests indicate that it is right.

Next step: implement it in sparseSET, test it, and if it is indeed
right... well, thats great news.

Results obtained for the new method:


```
  Now testing version 2 of assemble sparse line!!!
  25k lines, terms/lines: 100
  -> Time to assemble (s):   2.7800000000000047E-004
  100k lines, terms/lines: 200
  -> Time to assemble (s):   9.6000000000001640E-004
  100k lines, terms/lines: 600
  -> Time to assemble (s):   2.0219999999999683E-003
```

(or something was terrible wrong before, or now! Let's test it better!)
