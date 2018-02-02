%nprocshared=4
%chk=o.chk
#p opt=(calcfc,tight,noeigentest) freq=savenormalmodes ub3lyp/cc-pvdz
nosymm scf=qc geom=connectivity iop(7/33=1)

Title Card Required

0 3
 O                 -2.84776901    1.08923883    0.00000000

 1

