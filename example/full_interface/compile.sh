g++ -O3 \
  Full_interface.cpp \
  `/usr/local/bin/regina-engine-config --cflags --libs` \
  -I ../../src \
  -I ../../../knothomology \
  -o bundled-hfk-full-interface