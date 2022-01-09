g++ -O3 \
  Full_interface.cpp \
  `regina-engine-config --cflags --libs` \
  -I ../../src \
  -I ../../utilities/ComputeHFKv2 \
  -o bundled-hfk-full-interface