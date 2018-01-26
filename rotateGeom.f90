subroutine rotateGeom(numAtoms, Geom,numRotors, rotors, dihedrals, Mass, Geomout)

  !rotors contains other pivot atom, this pivot atom and num atoms in rotor and list of atoms

  implicit none
  double precision Geom(100,3), diheds(20), Mass(100), Geomout(100,3), pivotcoord(3)
  integer rotors(100,20), numRotors, i, numAtoms, j, k

 !copy the geometry
  do i=1,numAtoms
     do j=1,3
        Geomout(i,j) = Geom(i,j)
     end do
  end do

 !move the 
  do i=1,numRotors

     do j=1,3
        pivotcoord(j) = Geom(rotors(i,2),j)
     end do
     
     do j=1,numAtoms
        do k=1,3
           Geomout(j,k) = Geomout(j,k) - pivotcoord(k)
        end do
     end do
