#!/bin/csh

ncargf77 calc_dual_domain_vari_2.f
set nn=1
while($nn <= 8)
a.out << inp
$nn
1000
inp
mv gmeta 'domain'$nn'.cgm'
@ nn += 1
end
t
#Xgmeta2gif.csh
