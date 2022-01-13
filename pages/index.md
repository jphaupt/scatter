title: TODOs
---

@todo
In rough order:

- do a total cross section calculation, since that's relatively easy to do
(and you've done it before)
remember:
    - clean code
    - TDD
keep in mind (maybe do a bit later)
 - OOP if possible
 - function pointers? e.g. in case we want something other than Lenard-Jones
For this part, I will largely refer to *Computational Physics* by Thijssen.

- regression testing
https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html
https://stackoverflow.com/questions/19435375/simple-integration-testing-using-bash-with-cmake-and-ctest
probably use `add_test(... COMMAND ..)` and use a Python command to compare two (hopefully identical up to precision) text documents. 

- parallelise the cross section calculation (MPI probably), test with PFUnit

- extend to partial cross section

- extend to loss rates
@endtodo

## Extra/Low-Priority
@todo
some extra things to do that aren't really important:

- reach: embed required libraries (pFUnit, JSON-Fortran, stdlib?)

- Optional/default CMake options, `OPTION(FOO "Foo Option" OFF)`

- publish resulting code to github pages
@endtodo
