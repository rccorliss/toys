
#note:  this is saved to my repo for posterity in its current form, but can also be generated by ParseOGP.C if you set the 'build_gnuplot' flag to true
set terminal pdf size 60,30 linewidth 1 font 'Verdana,100'
set output 'all_spotchecks_multi_large.pdf'

#approximate values:
dx=-56.31;dy=-222.45;dphi=-1.2626


#petal 44
unset for [i=0:250] object i
set object 91 circle at   220.201,109.338 size first 5 fc rgb "orange" lw 15
set object 123 circle at   326.657,106.244 size first 12 fc rgb "red" lw 15
set object 123 circle at   326.657,106.244 size first 12 fc rgb "orange" lw 15
plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_44_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 44_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_44_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 44_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_44_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 44_3'


unset for [i=0:250] object i
unset for [i=0:250] object i

#petal 46
set object 116 circle at   316.560,43.681 size first 12 fc rgb "red" lw 15
set object 116 circle at   316.560,43.681 size first 12 fc rgb "orange" lw 15
set object 138 circle at   343.060,156.350 size first 12 fc rgb "orange" lw 15
set object 209 circle at   519.960,65.794 size first 12 fc rgb "red" lw 15

plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_46_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 46_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_46_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 46_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_46_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 46_3',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_46_4.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 46_4'

unset for [i=0:250] object i


plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_48_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 48_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_48_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 48_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_48_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 48_3'

unset for [i=0:250] object i


plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_49_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 49_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_49_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 49_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_49_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 49_3',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_49_4.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 49_4',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_49_5.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 49_5'

unset for [i=0:250] object i

#petal 74
set object 88 circle at   233.747,47.413 size first 5 fc rgb "red" lw 15
set object 88 circle at   233.747,47.413 size first 5 fc rgb "orange" lw 15
set object 94 circle at   246.183,83.625 size first 5 fc rgb "orange" lw 15
set object 140 circle at   368.287,154.271 size first 12 fc rgb "orange" lw 15


plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_74_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 74_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_74_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 74_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_74_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 74_3'
unset for [i=0:250] object i


plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_76_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 76_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_76_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 76_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_76_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 76_3',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_76_4.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 76_4',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_76_5.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 76_5'



plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_78_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 78_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_78_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 78_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_78_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 78_3'



plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_82_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 82_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_82_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 82_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_82_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 82_3'


plot  "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_87_1.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l title 'spotcheck 87_1',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_87_2.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 87_2',\
 "<grep 'Circle     7' '/Volumes/Enfain/Output\ Files/spotcheck_87_3.DAT' -A 350000" u (($1+dx)*cos(dphi)-($2+dy)*sin(dphi)):(($2+dy)*cos(dphi)+($1+dx)*sin(dphi)) w l  title 'spotcheck 87_3'
