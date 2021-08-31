# bundled-HFK
A `C++` library for calculating knot Floer homology via Ozsváth and Szabó's
bordered approach.

This code was written during my (Isaac Ren) stay at Inria Sophia
Antipolis-Méditerranée, doing an internship under the suprvision of Clément
Maria, from March to August 2021.

### Note on the license
This program is released under
[GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html), or (at your option)
any later version as published by the Free Software Foundation.

This program uses bits of Ozsváth and Szabó's program
[`ComputeHFKv2`](https://web.math.princeton.edu/~szabo/HFKcalc.html). The
original authors have a plan to release an updated version of the code on
Github. A slightly modified version of `ComputeHFKv2` also appears in a
`Python` module for [`SnapPy`](https://snappy.math.uic.edu/) under the name
[`knot_floer_homology`](https://github.com/3-manifolds/knot_floer_homology),
and is released under GPL v2.0.

### Running examples
The `bundled-HFK` library is made of only header files. In order to see the
code in action, we have written up compilable `.cpp` files in `/example`. From
any example folder, execute
```
sh compile.sh
```
The compiled file is usually named `bundled-hfk-example`. Executing the file
may require extra arguments; this is specified in a README file in the example
folder. Execute
```
./bundled-hfk-example [arguments]
```
to run the example.