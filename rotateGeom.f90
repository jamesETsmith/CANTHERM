SUBROUTINE rotateGeom(numAtoms, Geom,numRotors, rotors, dihedrals, Mass, Geomout)

  !rotors contains other pivot atom, this pivot atom and num atoms in rotor and list of atoms

  IMPLICIT NONE
  DOUBLE PRECISION Geom(100,3), diheds(20), Mass(100), Geomout(100,3), pivotcoord(3)
  INTEGER rotors(100,20), numRotors, i, numAtoms, j, k

  !copy the geometry
  DO i=1,numAtoms
     DO j=1,3
        Geomout(i,j) = Geom(i,j)
     END DO
  END DO

  !move the
  DO i=1,numRotors

     DO j=1,3
        pivotcoord(j) = Geom(rotors(i,2),j)
     END DO

     DO j=1,numAtoms
        DO k=1,3
           Geomout(j,k) = Geomout(j,k) - pivotcoord(k)
        END DO
     END DO
