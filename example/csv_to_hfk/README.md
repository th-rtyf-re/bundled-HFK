# CSV to HFK
### Running the example
Compile by executing
```
sh compile.sh
```
The compiled file `bundled-hfk-example` takes one argument, a `CSV` file
representing a knot diagram. The file should be of the form
```
position0,event0
position1,event1
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

It is not required to add the global minimum; a local minimum will be converted
to global if needed. However, it is required that all local minima appear in
position `0`.

As an example,
```
0,2
0,2
1,0
1,0
1,0
0,3
0,3
```
represents a right-handed trefoil knot. It can be found at
`../../data/csv/right_handed_trefoil.csv`. Thus, executing
```
./bundled-hfk-example ../../data/csv/right_handed_trefoil.csv
```
should produce the knot Floer homology `t^{-1}q^{-2} + q^{-1} + t`.

### Visualize the example
By default, the macro `DRAW` is defined in `csv_to_hfk.cpp`. If this is the
case, after executing the code, two files will be created:
`differential_suffix_tree.tex` and `knot_diagram.tex`. Compile the `LaTeX` file
`drawings.tex` (provided the studied knot is not too big) to visualize the
knot, as well as the differential suffix forest at each step of the
calculation.