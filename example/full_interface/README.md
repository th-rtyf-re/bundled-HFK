# Bundled HFK full interface
### Running the example

Before compiling,
[build Regina from source](https://regina-normal.github.io/source.html) and
download [`ComputeHFKv2`](https://web.math.princeton.edu/~szabo/HFKcalc.html)
to the `utilities` folder. To compile, run
```
sh compile.sh
```
**Note:** The file `compile.sh` looks like:
```
g++ -O3 \
  Full_interface.cpp \
  `regina-engine-config --cflags --libs` \
  -I ../../src \
  -I ../../utilities/ComputeHFKv2 \
  -o bundled-hfk-full-interface
```
The third line comes from
[Regina's website](https://regina-normal.github.io/buildtips.html#linking),
and may depend on where you install Regina.

The compiled file `bundled-hfk-full-interface` can run with several options.
Here is the options list:
```
-me, --morse-events <csv file>
-pd, --planar-diagram <txt file>
-rs, --regina-signature <csv file> [<knot number/start>] [<knot end>]
```

### Morse events

The option `-me` takes one argument, a `CSV` file representing a knot diagram.
The file should be of the form
```
event0,position0
event1,position1
...
```
where each line corresponds to a knot slice, and has two integers. The first
integer is the position of the event, and the second integer is the type of
the event. The possible events are:
- 0: positive crossing
- 1: negative crossing
- 2: local maximum
- 3: local minimum
- 4: global minimum

It is required that all local minima appear in position `0`.

For example,
```
2,0
2,0
0,1
0,1
0,1
3,0
4,0
```
represents a right-handed trefoil knot. It can be found at
`../../data/csv/right_handed_trefoil.csv`. Thus, executing
```
./bundled-hfk-example -me ../../data/csv/right_handed_trefoil.csv
```
should produce the knot Floer homology `t^{-1}q^{-2} + q^{-1} + t`.

### Planar diagram

The option `-pd` takes one argument, a text file containing a planar diagram
as described in Osváth and Szabó's program
[`ComputeHFKv2`](https://web.math.princeton.edu/~szabo/HFKcalc.html). See there
for a description of the syntax.

The download link on that page also provides examples of planar diagram text
files, which you can use.

### Regina signature

The option `-rs` takes a `CSV` file of Regina knot signatures, coming from
Benjamin Burton's [knot tables](https://regina-normal.github.io/data.html#knots).
You can download the knot tables and run the program on the files.

Without extra arguments, the program will compute the knot Floer homology of
each knot in the file. You can specify a start point, as well as an end point.
The start is inclusive and the end is exclusive.

### Visualize the example
By default, the macro `BUNDLE_HFK_DRAW_` is defined in `csv_to_hfk.cpp`. If
this is the case, after executing the code, two files will be created:
`differential_suffix_forest.tex` and `knot_diagrams.tex`. Compile the `LaTeX`
file `drawings.tex` (provided the considered knots are not too big) to
visualize all the considered knots, as well as the differential suffix forest
at each step of the calculation for the last knot. Note that this may reach
`LaTeX`'s limits even for relatively small knots, e.g. with 10 crossings.