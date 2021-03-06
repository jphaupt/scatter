title: TODOs
---

@todo
In rough order:

- [x] do a total cross section calculation, since that's relatively easy to do
(and you've done it before)
remember:
    - clean code
    - remember TDD (use this as a learning opportunity)
 - OOP if possible (and done)
 - Not using function pointers, but instead using class objects, which I think is a lot cleaner anyway. If the potential changes, we can simply write a new class.

- [ ] compare against the paper: https://aip.scitation.org/doi/pdf/10.1063/1.438414
    - I tried comparing against Thijssen's code and thought I did something wrong, but when I looked at the paper I realised that my cross-section was a lot better/closer to experiment than the solution given by the book... (probably because I make fewer approximations)
    - Just compare against FIG 3d (convert sigma to rho^2, and scale x-axis to be a function of E=1/2 mg^2)
    - easiest way to get data would probably just to OCR the data from the paper

- regression testing
https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html
https://stackoverflow.com/questions/19435375/simple-integration-testing-using-bash-with-cmake-and-ctest
probably use `add_test(... COMMAND ..)` and use a Python command to compare two (hopefully identical up to precision) text documents.
https://cmake.org/cmake/help/latest/command/exec_program.html
https://stackoverflow.com/questions/25437632/cmake-run-built-executable-before-building-another

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
