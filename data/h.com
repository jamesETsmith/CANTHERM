%nprocshared=4
%chk=h.chk
#p opt=(calcfc,tight,noeigentest) freq=savenormalmodes ub3lyp/cc-pvdz
nosymm scf=qc geom=connectivity iop(7/33=1)

Title Card Required

0 2
 H                  4.81627293    0.77427820    0.00000000

 1

