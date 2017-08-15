#!/bin/csh 
foreach aa (*.cgm )
#foreach aa (gmeta)
#foreach aa (*PPI*)
#ctrans -d xwd -resolution 1500x1500 $aa > $aa'.xwd'
/opt/ncl-6.1.2-gnu/bin/ctrans -d xwd -resolution 2048x2048 $aa > $aa'.xwd'
#ctrans -d xwd -resolution 1024x1024 $aa > $aa'.xwd'
#ctrans -d xwd $aa > $aa'.xwd'
convert $aa'.xwd' $aa'.gif'

#rm -f $aa'.xwd'
rm -f $aa'.xwd'  $aa
end
