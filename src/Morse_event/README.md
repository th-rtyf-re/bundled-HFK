# Morse events
A Morse event is the given of a knot slice and its associated DA-bimodule.

### How to add a Morse event
It is possible to add any type of knot slice to our knot diagrams, as long as
the DA-bimodule structure is implementable. For example, we can add a Morse
event for marked minima in order to calculate link Floer homology.

In order to add a Morse event, follow these steps:
1. Add an item to the Event enum class in `Morse_event.h`
2. Write a `.h` file with a `Knot_slice`-derived class and a
`DA_bimodule`-derived class, respectively implementing the topological
structure and the algebraic structure of the Morse event.
3. In `Event_registers.h`, add the `Knot_slice`-derived class to
`Knot_slice_register` and the `DA_bimodule`-derived class to
`DA_bimodule_register`.