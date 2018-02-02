%nprocshared=4
%chk=s.chk
#p opt=(calcfc,tight,noeigentest) freq=savenormalmodes ub3lyp/cc-pvdz
nosymm scf=qc geom=connectivity iop(7/33=1)

Title Card Required

0 1
 S                 -2.95275588   -0.06561680    0.00000000

 1

