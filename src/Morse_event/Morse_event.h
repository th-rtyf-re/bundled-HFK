#ifndef MORSE_EVENT_H_
#define MORSE_EVENT_H_

/* Morse events.
 *
 * This enum class lets us declare the various knot slices we use in our knot
 * diagrams.
 */
enum class Morse_event {
  positive_crossing,
  negative_crossing,
  local_maximum,
  local_minimum,
  global_minimum
};

#endif  // MORSE_EVENT_H_