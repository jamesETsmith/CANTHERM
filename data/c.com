%nprocshared=4
%chk=c.chk
#p opt=(calcfc,tight,noeigentest) freq=savenormalmodes ub3lyp/cc-pvdz
nosymm scf=qc geom=connectivity iop(7/33=1)

Title Card Required

0 3
 C                 -1.53543306    1.69291336    0.00000000

 1

