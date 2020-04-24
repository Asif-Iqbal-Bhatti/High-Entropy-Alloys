program BCC_Dislocation

!!! USAGE ::: gfortran -o test dislocation.f90

implicit none
!------------------------- 
integer                :: atomspercel, ncelltot, numatomstot, latt_par
integer                :: i, j, k, l, NumCellP, iun, m, numhide
integer,dimension(4)   :: N
real(kind=8)           :: R_2, R_3, R_6, delta
real(kind=8)           :: burgers, bcc_lat
real(kind=8)           :: pivalue, pisurtrois, factor, pioversix
real(kind=8)           :: radius, XDislo, YDislo, xx, yy, zz
real(kind=8)           :: ux, uy, theta1, theta2, distance
real(kind=8),dimension(3)                 :: C1, C2
real(kind=8),dimension(3)                 :: unit_lat
real(kind=8),dimension(3)                 :: UnitCellDim
real(kind=8),dimension(3,2)               :: dimbox
real(kind=8),dimension(3,3)               :: dimbox_pos
real(kind=8),dimension(:,:), allocatable  :: cell_coord
real(kind=8),dimension(:,:,:),allocatable :: AtomsTot

!------------------------- PARAMETERS for (111) unit cell
R_3 = sqrt(3.0d0) ! [111]
R_6 = sqrt(6.0d0) ! [112]
R_2 = sqrt(2.0d0) ! [110]
pivalue = 2.0d0*asin(1.0d0) !* PI value encoded in radians
pisurtrois = pivalue / 3.0d0
pioversix = pivalue / 6.0d0

! -- Enter the Ta or Nb lattice parameter obtained from DFT
bcc_lat = 3.31953d0
!------------------------- 

! -- Definition of unit cell
atomspercel = 6
allocate(cell_coord(atomspercel,3))
unit_lat(1) = bcc_lat*R_6       ! a*[112]=X
unit_lat(2) = bcc_lat*R_2       ! a*[110]=Y
unit_lat(3) = bcc_lat*R_3/2.0d0 ! a*[111]=Z --> a*<111>/2

! -- Atomic positions within a unitcell (X,Y,Z)
cell_coord(1,1)=unit_lat(1)*0.0d0        ; cell_coord(1,2)=unit_lat(2)*0.0d0  ; cell_coord(1,3)=unit_lat(3)*0.0d0;
cell_coord(2,1)=unit_lat(1)*(1.0d0/2.0d0); cell_coord(2,2)=unit_lat(2)*(0.5d0); cell_coord(2,3)=unit_lat(3)*(0.0d0)
cell_coord(3,1)=unit_lat(1)*(1.0d0/3.0d0); cell_coord(3,2)=unit_lat(2)*(0.0d0); cell_coord(3,3)=unit_lat(3)*(2.0d0/3.0d0)
cell_coord(4,1)=unit_lat(1)*(5.0d0/6.0d0); cell_coord(4,2)=unit_lat(2)*(0.5d0); cell_coord(4,3)=unit_lat(3)*(2.0d0/3.0d0)
cell_coord(5,1)=unit_lat(1)*(1.0d0/6.0d0); cell_coord(5,2)=unit_lat(2)*(0.5d0); cell_coord(5,3)=unit_lat(3)*(1.0d0/3.0d0)
cell_coord(6,1)=unit_lat(1)*(2.0d0/3.0d0); cell_coord(6,2)=unit_lat(2)*(0.0d0); cell_coord(6,3)=unit_lat(3)*(1.0d0/3.0d0)

!------------------------- Size of the unit cell !-------------------------

UnitCellDim(1)=unit_lat(1) 
UnitCellDim(2)=unit_lat(2)
UnitCellDim(3)=unit_lat(3)
N(1)=5 
N(2)=8 
N(3)=1
write(*,*)'X=[112]',unit_lat(1),'Y=[110]', unit_lat(2)
write(*,*)'Supercell scaling (X,Y,Z)', N(1), N(2), N(3)
write(*,*)'Scaling unit cells along X,', UnitCellDim(1)*N(1)
write(*,*)'Scaling unit cells along Y,', UnitCellDim(2)*N(2)
write(*,*)'Scaling unit cells along Z,', UnitCellDim(3)*N(3)
ncelltot = N(1)*N(2)*N(3)

allocate(AtomsTot(ncelltot,atomspercel,3))

!------------------------- Generate the atomic positions for a supercell
NumCellP = 0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   NumCellP = NumCellP + 1
   do l=1,atomspercel
     AtomsTot(NumCellP,l,1) = cell_coord(l,1) + (i-1)*UnitCellDim(1)
     AtomsTot(NumCellP,l,2) = cell_coord(l,2) + (j-1)*UnitCellDim(2)
     AtomsTot(NumCellP,l,3) = cell_coord(l,3) + (k-1)*UnitCellDim(3)
   enddo
  enddo
 enddo
enddo

!#***************Position of dislocation line at C1(X, Y, Z) for +b*********
write(*,'(a)') ''   
write(*,'(a)') 'Generating Screw Dislocations >>>'   

C1(:) = 0.0d0; delta = 0.001D-01 ! to avoid on top of atom
C1(1) = 0.30d0*N(1)*UnitCellDim(1) - delta !** X
C1(2) = 0.26d0*N(2)*UnitCellDim(2) 					!** Y
C1(3) = N(3)*UnitCellDim(3) 								!** Z

write(*,'(a)') 'Position of dislocation line at C1(X, Y, Z) for +b >>>'
write(*,*) C1(1),C1(2),C1(3)
!#**************************************************************************
!#***************Position of dislocation line at C2(X, Y, Z) for -b*********
C2(:) = 0.0d0; delta = 0.001D-01 ! to avoid on top of atom
C2(1) = 0.59d0*N(1)*UnitCellDim(1) - delta !** X
C2(2) = 0.72d0*N(2)*UnitCellDim(2) 					!** Y
C2(3) = N(3)*UnitCellDim(3) 								!** Z

write(*,'(a)') 'Position of dislocation line at C2(X, Y, Z) for -b >>>'
write(*,*) C2(1),C2(2),C2(3)
!#*************************************************************************

distance = sqrt( (C1(1)-C2(1))**2 + (C1(2)-C2(2))**2 )
write(*,*) 'Distance between Dipoles >>> ', distance

burgers = bcc_lat*R_3/2.0 ! b=a/2<111>

NumCellP = 0
numhide=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   NumCellP = NumCellP + 1
   do l=1,6
	 
    !* tan**(-1)(x/y)
    xx = AtomsTot(NumCellP,l,1)
    yy = AtomsTot(NumCellP,l,2)
		
    theta1 = datan2((xx-C1(1)),(yy-C1(2) ) )
    theta2 = datan2((xx-C2(1)),(yy-C2(2) ) )
		
    !*** screw dislocation displacement in Z direction
    AtomsTot(NumCellP,l,3) = AtomsTot(NumCellP,l,3) + burgers/2.0d0/pivalue*(theta1)
    AtomsTot(NumCellP,l,3) = AtomsTot(NumCellP,l,3) - burgers/2.0d0/pivalue*(theta2)

   enddo
  enddo
 enddo
enddo   

!#****************************** Writing to a file ************************
dimbox(:,:)=0.0d0; dimbox_pos(:,:)=0.0d0
dimbox(1,2) = N(1)*UnitCellDim(1); dimbox_pos(1,1) = N(1)*UnitCellDim(1)
dimbox(2,2) = N(2)*UnitCellDim(2); dimbox_pos(2,2) = N(2)*UnitCellDim(2)
dimbox(3,2) = N(3)*UnitCellDim(3); dimbox_pos(3,3) = N(3)*UnitCellDim(3)

!#*************************************************************************
open(1,file='Ta_lammps.lmp',status='REPLACE')
write(1,*) 'Position data for Fe File'
write(1,*) 
write(1,*) 6*ncelltot - numhide, ' atoms'
write(1,*) ' 1 atom types'
write(1,*) dimbox(1,1:2),' xlo xhi'
write(1,*) dimbox(2,1:2),' ylo yhi'
write(1,*) dimbox(3,1:2),' zlo zhi'
write(1,*)
write(1,*) 'Atoms'
write(1,*)
m=1
iun=1
NumCellP=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   NumCellP = NumCellP+1
   do l=1,6
    write(1,10) m, iun, AtomsTot(NumCellP,l,1:3)
    m  = m+1
   enddo
  enddo
 enddo
enddo
close(1)
!#*************************************************************************
open(2,file='POSCAR_test',status='REPLACE')
write(2,'(a)') 'Ta screw'
write(2,'(a)') '1.0'
write(2,*) dimbox_pos(1,1:3)
write(2,*) dimbox_pos(2,1:3)
write(2,*) dimbox_pos(3,1:3)
write(2,'(a)') 'Ta'
write(2,*) 6*ncelltot - numhide
write(2,'(a)')'Cartesian'

NumCellP=0
do i=1,N(1)
 do j=1,N(2)
  do k=1,N(3)
   NumCellP = NumCellP+1
   do l=1,6
    write(2,20) AtomsTot(NumCellP,l,1:3)
   enddo
  enddo
 enddo
enddo
close(2)
!#*************************************************************************

deallocate(cell_coord)
deallocate(AtomsTot)


10 format(i8, i8, 3(1x, e12.5))
20 format(3(1x, e15.9))

end program
