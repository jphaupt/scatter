# scatter

Quantum scattering, and (eventually) loss rate calculations for elastic collisions. 

## Instructions

Please download and install: 
- PFUnit
  - For now this is mandatory but I plan to make it optional in case you do not wish to enable testing.
- JSON-Fortran
  - For now you must build this separately. Please make sure you select the right precision. At some point, I want to embed this into my project, but as a hacky temporary workaround, I'm going to just ask you to do it yourself. :)
  - ```
  mkdir build && cd build
  cmake -DSKIP_DOC_GEN=TRUE .. # optionally also e.g. -DJSON_REAL_KIND=REAL32
  make
  make install
  ```
  optionally also add to your `PATH`. Maybe also use `-DCMAKE_INSTALL_PREFIX=...`

You can also get automatically-generated documentation using FORD:
```
ford scatter_doc.md
```