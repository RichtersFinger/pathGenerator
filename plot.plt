set term pngcairo size 960,960

set size ratio -1
set xr [*<-0.05:1.05<*]
set yr [*<-0.05:1.05<*]

unset tics; unset label; unset border; unset grid; unset key

set obj 1 circ at 0,0 size 0.025 fs solid fc rgb '#000000' front
set obj 2 circ at 1,1 size 0.025 fs solid fc rgb '#000000' front

set output 'path.png'

pl 'output.dat' w li lw 2 lc rgb '#0060ad'

unset output