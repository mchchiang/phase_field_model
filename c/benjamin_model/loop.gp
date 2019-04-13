set xrange [1:100]
set yrange [1:100]
set cbrange [0:3]
#set palette defined (0 "white", 2 "blue")
p 'output.dat' matrix w image, 'cm.dat' u 3:2 w points lc "white" lt 7
while (1) {
  stat = system("wc output.dat | awk '{print $1}'")
  if (stat > 0 ) {
    replot
  }
  pause 0.3
}