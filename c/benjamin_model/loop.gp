set xrange [1:100]
set yrange [1:100]
set cbrange [0:3]
#set palette defined (0 "white", 2 "blue")
p 'field.dat' u 1:2:3 w image, 'cm.dat' every ::2 u 1:2 w points lc "black" lt 7
while (1) {
  pause 0.3
  replot
}