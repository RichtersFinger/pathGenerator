
# Simple Path Generator
The algorithm in `makerandompath.cpp` uses a composition of Fourier modes to generate a pseudo-random smooth path between the two points (0,0) and (1,1) in the plane.

Compile this using the command
   ```console
   $ g++ makerandompath.cpp -o run
   ```
and plot the result with the provided gnuplot script by entering
   ```console
   $ gnuplot plot.plt
   ```
after execution.

Example gallery:
![](/demo/examples.png?raw=true)
