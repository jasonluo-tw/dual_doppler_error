#!/bin/csh

#ctrans -d ps.color -f font12 $1 > $1.ps
#psplit $1'.ps' $1'_'

#foreach gmt ($1'_'*.ps )
#convert $gmt $gmt'.gif'
#end


foreach cgm(*.cgm)
/opt/ncl-6.1.2-gnu/bin/ctrans -d ps.color -font $cgm > $cgm'.ps'
psplit $cgm'.ps' $cgm'_'
end

foreach gmt (*_*.ps )

convert $gmt $gmt'.gif'
end

rm -f *.ps
rm -f *.cgm
rm -f *.cgm.ps.*
