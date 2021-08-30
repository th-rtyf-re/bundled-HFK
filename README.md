# bundled-HFK
A `C++` library for calculating knot Floer homology via Ozsváth and Szabó's
bordered approach.

This code was written during my (Isaac Ren) stay at Inria Sophia
Antipolis-Méditerranée, doing an internship under the suprvision of Clément
Maria, from March to August 2021.

### Note on the license
This code uses bits of Ozsváth and Szabó's program
[`ComputeHFKv2`](https://web.math.princeton.edu/~szabo/HFKcalc.html). They did
not specify a license for distribution, so I chose to distribute this code
under GPL v3.0. I have contacted the original authors, and when they reply, I
will update the license accordingly.

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