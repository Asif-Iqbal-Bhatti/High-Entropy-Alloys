program BCC_ScrewDislocation
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!!! USAGE ::: gfortran -o test dislocation.f90; ./test
!!! AUTHOR::: Asif Iqbal
!!! DATED ::: 28/04/2020
!!! GITHUB::: @asif_em
!!! USE AT YOUR OWN RISK. NOT EVEN IMPLIED WARRANTY WHATSOEVER
!!! CAREFULLY CHECK THE GEOMETRY BEFORE SUBMITTING TO DFT CALCULATION.
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

IMPLICIT NONE
!------------------------- 
integer                :: atom_cell, ncelltot, numatomstot, latt_par
integer                :: i, j, k, l, image_cell, iun, m, p, q
real(kind=8)           :: R_2, R_3, R_6, delta, r, d
real(kind=8)           :: burgers, bcc_lat
real(kind=8)           :: pivalue, factor
real(kind=8)           :: xx, yy, zz
real(kind=8)           :: ux, uy, theta1,theta2,theta3,theta4
real(kind=8)           :: d1, d2, d3, d4
integer,dimension(4)   :: N
real(kind=8),dimension(3)                 :: L0, L1, diag
real(kind=8),dimension(3)                 :: C1, C2, C3, C4
real(kind=8),dimension(3)                 :: unit_lat, supercell
real(kind=8),dimension(3,2)               :: dimbox
real(kind=8),dimension(3,3)               :: dimbox_pos
real(kind=8),dimension(:,:), allocatable  :: Cart_cord
real(kind=8),dimension(:,:,:),allocatable :: Tot_atom, Tot_perfect

!------------------------- PARAMETERS for (111) unit cell
R_6 = sqrt(6.0d0) ! [112]
R_2 = sqrt(2.0d0) ! [110]
R_3 = sqrt(3.0d0) ! [111]

!--- PI value encoded in radians
pivalue = 2.0d0*asin(1.0d0) 
!--- Enter the Ta or Nb lattice parameter obtained from DFT
bcc_lat = 3.31953d0
!------------------------- 
call random_number(r)

! -- Definition of unit cell
atom_cell = 6
allocate(Cart_cord(atom_cell,3))
unit_lat(1) = bcc_lat*R_6       ! a*[1-1-2]=X
unit_lat(2) = bcc_lat*R_2       ! a*[110]=Y
unit_lat(3) = bcc_lat*R_3/2.0d0 ! a*[1-11]=Z --> a*<111>/2

! -- Atomic positions within a unitcell (X,Y,Z) in Fractional coordinates
Cart_cord(1,1)=unit_lat(1)*0.0d0        ; Cart_cord(1,2)=unit_lat(2)*0.0d0  ; Cart_cord(1,3)=unit_lat(3)*0.0d0;
Cart_cord(2,1)=unit_lat(1)*(1.0d0/2.0d0); Cart_cord(2,2)=unit_lat(2)*(0.5d0); Cart_cord(2,3)=unit_lat(3)*(0.0d0)
Cart_cord(3,1)=unit_lat(1)*(1.0d0/3.0d0); Cart_cord(3,2)=unit_lat(2)*(0.0d0); Cart_cord(3,3)=unit_lat(3)*(2.0d0/3.0d0)
Cart_cord(4,1)=unit_lat(1)*(5.0d0/6.0d0); Cart_cord(4,2)=unit_lat(2)*(0.5d0); Cart_cord(4,3)=unit_lat(3)*(2.0d0/3.0d0)
Cart_cord(5,1)=unit_lat(1)*(1.0d0/6.0d0); Cart_cord(5,2)=unit_lat(2)*(0.5d0); Cart_cord(5,3)=unit_lat(3)*(1.0d0/3.0d0)
Cart_cord(6,1)=unit_lat(1)*(2.0d0/3.0d0); Cart_cord(6,2)=unit_lat(2)*(0.0d0); Cart_cord(6,3)=unit_lat(3)*(1.0d0/3.0d0)

!------------------------- Size of the unit cell !-------------------------
!!!              Change supercell vector according to your need         !!!
N(1)=7
N(2)=11  ! odd
N(3)=1

supercell(1)=N(1)*unit_lat(1) 
supercell(2)=N(2)*unit_lat(2)
supercell(3)=N(3)*unit_lat(3)

print'("X=a*[112]; ", "Y=a*[110]; ", "Z=a*[111]/2; ", 3F12.6)',unit_lat(1), unit_lat(2),unit_lat(3) 
print'("Supercell scaling N(X,Y,Z)", 3I9)', N(1), N(2), N(3)
print'("Scaling unit cells along NX,", F15.9)', supercell(1)
print'("Scaling unit cells along NY,", F15.9)', supercell(2)
print'("Scaling unit cells along NZ,", F15.9)', supercell(3)

ncelltot = N(1)*N(2)*N(3)
allocate(Tot_atom(ncelltot,atom_cell,3))
allocate(Tot_perfect(ncelltot,atom_cell,3))

!------------------------- Generate atomic positions for a supercell
image_cell = 0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   image_cell = image_cell + 1
   do l=1,atom_cell
     Tot_atom(image_cell,l,1) = Cart_cord(l,1) + (i-1)*unit_lat(1) ! Xi
     Tot_atom(image_cell,l,2) = Cart_cord(l,2) + (j-1)*unit_lat(2) ! Yi
     Tot_atom(image_cell,l,3) = Cart_cord(l,3) + (k-1)*unit_lat(3) ! Zi
   enddo
  enddo
 enddo
enddo
Tot_perfect = Tot_atom
print'("Total Number of atoms in the supercell", I8)', image_cell*6

print'(" ")'   
write(*,'(a)') 'Generating Screw Dislocations (dislo) at Position of dislocation line >>>'
!#***************Position of dislocation line at C1(X, Y, Z) for +b*********
C1(:) = 0.0d0; delta = 0.001D-01 ! to avoid on top of atom
C1(1) = 0.251d0*supercell(1) - delta 
C1(2) = 0.251d0*supercell(2) 			
C1(3) = supercell(3) 					
!#***************Position of dislocation line at C2(X, Y, Z) for -b*********
C2(:) = 0.0d0;
C2(1) = 0.751d0*supercell(1) + delta
C2(2) = 0.251d0*supercell(2) 				
C2(3) = supercell(3) 				
!#***************Position of dislocation line at C3(X, Y, Z) for -b*********
C3(:) = 0.0d0; 
C3(1) = 0.251d0*supercell(1) - delta 
C3(2) = 0.751d0*supercell(2) 				
C3(3) = supercell(3) 				
!#***************Position of dislocation line at C4(X, Y, Z) for +b*********
C4(:) = 0.0d0;
C4(1) = 0.751d0*supercell(1) + delta 
C4(2) = 0.751d0*supercell(2) 				
C4(3) = supercell(3) 			

print'("C1(X,Y,Z) for +b",3F14.6)', C1(1),C1(2),C1(3)
print'("C2(X,Y,Z) for -b",3F14.6)', C2(1),C2(2),C2(3)
print'("C3(X,Y,Z) for +b",3F14.6)', C3(1),C3(2),C3(3)
print'("C4(X,Y,Z) for -b",3F14.6)', C4(1),C4(2),C4(3)
!----------------------------------------------------------------------------

d1 = sqrt( (C1(1)-C2(1))**2 + (C1(2)-C2(2))**2 )
d2 = sqrt( (C1(1)-C3(1))**2 + (C1(2)-C3(2))**2 )
d3 = sqrt( (C1(1)-C4(1))**2 + (C1(2)-C4(2))**2 )
print'("d1=|C1-C2|, d2=|C1-C3|, d3=|C1-C4| >>> ",3F12.6 )', d1, d2, d3
print'("-Dist along NX-d1,(NX-d1)/2 >>> ",2F14.6 )', supercell(1)-d1, (supercell(1)-d1)/2
print'("-Dist along NY-d2,(NY-d2)/2 >>> ",2F14.6 )', supercell(2)-d2, (supercell(2)-d2)/2
print'(" ",2F14.6 )', C1(1) - sqrt( 28.459**2 - (C1(2) - 38.7815)**2 )

print'("-------------------------------------------------------------------------")'
print'("Constructing Edge lattice from unit cell >>> ")'
print'("Ju Li et al paper")'
print'("-> h1 = 7e1 >>> ",F20.6)'               , 7*unit_lat(1)
print'("-> h2 = 3.5e1+5.5e2+0.5e3 >>> ",3F14.6)', 3.5*unit_lat(2), 5.5*unit_lat(2), 0.5*unit_lat(2)
print'("-> h3 = 2e3 >>> ",F20.6)'               , 2*unit_lat(3) ! For HEAs twice the thickness
p=15; q=9
print'("L. Ventelon et al paper")'
print'("-> C1 = pe1-(1/3q)e3 >>> ",2F20.6)' , p*unit_lat(1), -(1/(3*q))*unit_lat(3)
print'("-> C2 = 0.5pe1+qe2+(0.5-1/6q)e3 >>> ",3F14.6)', 0.5*p*unit_lat(1), q*unit_lat(2), (0.5-1/(6*q))*unit_lat(3)
print'("-> C3 = e3 >>> ",F20.6)'               , unit_lat(3) ! For HEAs twice the thickness
print'("-------------------------------------------------------------------------")'

!----------------------------------------------------------------------------
!A periodic array is quadrupole, if the vector d linking the two disloca-
!tions of opposite signs is equal to 1/2 (u1 +u2), where u1 and u2 are the periodicity
!vectors of the simulation cell. This ensures that every dislocation is a sym-
!metry center of the array: fixing, as a convention, the origin at a dislocation center,
!if a dislocation b is located at the position r, there will also be a dislocation b in −r.

d = (supercell(1) + supercell(2))/2.0
print'("Vector d linking the two +b & -b  >>> ",2F14.6)', d, d-d1

L0=0.0; L1=0.0
L0(1) = supercell(1)
L1(2) = supercell(2)
diag = 0.5*(dot_product(L0,L0)+dot_product(L0,L1))/(dot_product(L0,L0)+dot_product(L1,L1)+4*dot_product(L0,L1))
	
! --------------------Screw dislocation displacement------------------------
burgers = bcc_lat*R_3/2.0 ! b=a<111>/2
image_cell = 0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   image_cell = image_cell + 1
   do l=1,atom_cell
	 
    !--- (b/2pi)*tan**(-1)(y/x) ""DEF:: numpy.arctan2(x1, x2) x1 = y; x2 = x ""
    xx = Tot_atom(image_cell,l,1)
    yy = Tot_atom(image_cell,l,2)
		
		! angle from X -> Y This is right
    theta1 = datan2( (yy-C1(2) ), (xx-C1(1)) )
    theta2 = datan2( (yy-C2(2) ), (xx-C2(1)) )
    theta3 = datan2( (yy-C3(2) ), (xx-C3(1)) )
    theta4 = datan2( (yy-C4(2) ), (xx-C4(1)) )
	
    !theta1 = datan2( (xx-C1(1)), (yy-C1(2) ) ) ! angle from Y -> X
    !theta2 = datan2( (xx-C2(1)), (yy-C2(2) ) )
    !theta3 = datan2( (xx-C3(1)), (yy-C3(2) ) )
    !theta4 = datan2( (xx-C4(1)), (yy-C4(2) ) )
		
    !--- Screw dislocation displacement in Z direction
    Tot_atom(image_cell,l,3) = Tot_atom(image_cell,l,3) + (burgers/(2.0d0*pivalue))*(theta1)
    Tot_atom(image_cell,l,3) = Tot_atom(image_cell,l,3) - (burgers/(2.0d0*pivalue))*(theta2)
    Tot_atom(image_cell,l,3) = Tot_atom(image_cell,l,3) - (burgers/(2.0d0*pivalue))*(theta3)
    Tot_atom(image_cell,l,3) = Tot_atom(image_cell,l,3) + (burgers/(2.0d0*pivalue))*(theta4)

   enddo
  enddo
 enddo
enddo 

!----------------------------------------------------------------------------------- 
!#********************************* Writing to a file ******************************
!-----------------------------------------------------------------------------------
dimbox(:,:)=0.0d0; dimbox_pos(:,:)=0.0d0
dimbox(1,2) = supercell(1); dimbox_pos(1,1) = supercell(1)
dimbox(2,2) = supercell(2); dimbox_pos(2,2) = supercell(2)
dimbox(3,2) = supercell(3); dimbox_pos(3,3) = supercell(3)

!#*************************************************************************
open(1,file='TaScrew_unrelaxed.lmp',status='REPLACE')
write(1,*) 'Position data for Fe File'
write(1,*) 
write(1,*) 6*ncelltot , ' atoms'
write(1,*) ' 1 atom types'
write(1,*) dimbox(1,1:2),' xlo xhi'
write(1,*) dimbox(2,1:2),' ylo yhi'
write(1,*) dimbox(3,1:2),' zlo zhi'
write(1,*)
write(1,*) 'Masses'
write(1,*)
write(1,*) '1 180.94788  # Ta'
write(1,*)
write(1,*) 'Atoms'
write(1,*)
m=1
iun=1
image_cell=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   image_cell = image_cell+1
   do l=1,6
    write(1,10) m, iun, Tot_atom(image_cell,l,1:3)
    m  = m+1
   enddo
  enddo
 enddo
enddo
close(1)

!#********************* POSCARScrew_unrelaxed *****************************
print'(" >>> Files have been generated ... <<< ")'
open(2,file='POSCARScrew_unrelaxed',status='REPLACE')
write(2,'(a)') 'Ta screw'
write(2,'(F9.6)') bcc_lat ! bcc_lat = 3.31953d0
write(2,*) dimbox_pos(1,1:3)/bcc_lat
write(2,*) dimbox_pos(2,1:3)/bcc_lat
write(2,*) dimbox_pos(3,1:3)/bcc_lat
write(2,'(a)') 'Ta'
write(2,'(I7)') 6*ncelltot 
write(2,'(a)')'Cartesian'

image_cell=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   image_cell = image_cell+1
   do l=1,6
    write(2,20) Tot_atom(image_cell,l,1:3)/bcc_lat
   enddo
  enddo
 enddo
enddo
close(2)

!#********************* POSCARScrew_perfect *****************************
open(3,file='POSCAR_perfect',status='REPLACE')
write(3,'(a)') 'Ta '
write(3,'(F9.6)') bcc_lat ! bcc_lat = 3.31953d0
write(3,*) dimbox_pos(1,1:3)/bcc_lat
write(3,*) dimbox_pos(2,1:3)/bcc_lat
write(3,*) dimbox_pos(3,1:3)/bcc_lat
write(3,'(a)') 'Ta'
write(3,'(I7)') 6*ncelltot 
write(3,'(a)')'Cartesian'

image_cell=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   image_cell = image_cell+1
   do l=1,6
    write(3,20) Tot_perfect(image_cell,l,1:3)/bcc_lat
   enddo
  enddo
 enddo
enddo
close(2)

!#*************************************************************************
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

deallocate(Cart_cord)
deallocate(Tot_atom)


10 format(i8, i8, 3(1x, e12.5))
20 format(3(1x, e15.9))


end program


	
	
