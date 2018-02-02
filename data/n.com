%nprocshared=4
%chk=n.chk
#p opt=(calcfc,tight,noeigentest) freq=savenormalmodes ub3lyp/cc-pvdz
nosymm scf=qc geom=connectivity iop(7/33=1)

Title Card Required

0 2
 N                  1.03674540   -1.08923883    0.00000000

 1

