! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OLIVINEABS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! XSPEC local model for olivine absorption edge structure
! Contains ISMdust cross-sections for absorption with Fe-K edge for
! olivine, obtained from Rogantini et al. 2018
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine olivineabs(ear, ne, param, ifl, photar)
!
! The main routine to call all subroutines
!
implicit none
integer,parameter :: num_param = 2
integer,parameter :: nemod=13047 !Number of elements for each cross section.
integer :: ne, ifl, a
double precision :: moliv, rshift, emod(nemod), coemod(nemod)
double precision :: bxs(nemod), bener(nemod)
double precision :: zfac
real :: ear(0:ne), param(num_param), photar(ne)
logical :: startup=.true.
character (len=40) version
version='0.1'
 if(startup)then
  print *, ' '
  print *, 'OlivineAbs: High resolution model for silicate with Olivine Fe K',version
  print *, 'Continuum xsect from ISMdust silicate from Corrales+ 2016 (MNRAS, 458, 1345)'
  print *, 'Optical constants come from Draine 2003 (ApJ, 598, 1026) except'
  print *, 'Fe K absorption, which comes from Rogantini+ 2018 (A&A, 609, A22)'
  print *, 'WARNING: If used in conjunction with neutral metal absorption models'
  print *, '(e.g. TBabs, TBnew), be sure to change abundances'
  print *, ' '
  call read_cross_sections_olivine(nemod,bxs,ifl,bener)
  startup=.false.
 endif
! Model parameters
moliv = param(1)
rshift = param(2)
zfac = 1.d0/(1.d0+dble(rshift))

call absorption_olivine(moliv, zfac, emod, nemod, coemod,bxs,ifl,bener)
!
call map_to_grid_olivine(dble(ear),ne,emod,nemod,photar,coemod,ifl)
return
end subroutine olivineabs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_cross_sections_olivine(bnene,xabs,ifl,ener)
!
! This routine reads cross sections and puts them on a given grid
!
! It uses the X-Spec local variable/dictionary value ISMDUSTROOT
! to locate the data file. If this is not set then it uses
! the setting of the local_dir parameter below (which should be
! changed to match the location of the data file). The path -
! however it is given - should point to the location that contains
! the dust extinction templates, i.e. /path/to/ismdust/edge_files/
!
implicit none
integer,parameter :: out_unit=20
integer :: bnene, ifl, i, j, status
integer :: nemax
double precision :: ener(bnene), xabs(bnene)
character (*), parameter :: fileloc = 'olivine_abs.fits'
character (*), parameter :: olivinereadchat = 'olivineabs: reading from '
character (len=255 + 29) :: filename2 ! len(fileloc)
character (len=240) :: local_dir = './edge_files/'
character (len=255) :: ismdust_root = ''
character (len=len(olivinereadchat)+len(filename2)) :: chatmsg = ''
integer inunit,readwrite,blocksize
integer :: hdutype,colnum
integer :: felem=1, nulld=0
logical :: anynull
character (len=255) :: fgmstr
external :: fgmstr
nemax=13047 !Number of elements for each grain type cross section.
! Where do we look for the data?
ismdust_root = trim(fgmstr('ISMDUSTROOT'))
if (ismdust_root .EQ. '') then
ismdust_root = local_dir
endif
! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ismdust_root) // fileloc
chatmsg=olivinereadchat // filename2
call xwrite(chatmsg, 20)
! Get an unused Logical Unit Number to use to open the FITS file.
call ftgiou(inunit,status)
! Open the FITS file
call ftopen(inunit,filename2,readwrite,blocksize,status)
! Move to the extension 2 (the binary table)
call ftmahd(inunit,2,hdutype,status)

!Read in one energy grid (column 1)
colnum=1
do j=1,nemax
  call ftgcvd(inunit,colnum,j,felem,1,nulld,ener(j),anynull,status)
enddo

!Read in the absorption cross section information (column 3)
colnum=3
do j=1,nemax
  call ftgcvd(inunit,colnum,j,felem,1,nulld,xabs(j),anynull,status)
  !print *, xabs(j), ener(j)
enddo

! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the X-Spec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used. An alternative would be to use xwrite()
! with a low chatter level.
!
! This message could be displayed only once, but it is probaly worth
! repeating each time it is used.
if (status .ne. 0) then
write (*,*) 'ERROR: unable to read cross sections from ', filename2
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_cross_sections_olivine
! ======================================= !
subroutine absorption_olivine(moliv, zfac, e1, bnene, coeff, bxs2,ifl,bener)
!
! This is routine that calculates the optical depth given the column densities
! Finally returns the absorption coefficient exp(-tau)
!
implicit none
integer,parameter :: out_unit=20
integer :: bnene, ifl
integer :: i
double precision :: moliv
double precision :: bener(bnene), bxs2(bnene), e1(bnene)
double precision :: tau, coeff(bnene)
double precision :: zfac
real hphoto, gphoto
external hphoto, gphoto

! Calculates the optical depth and the extinction coefficient exp(-tau)
e1(1)=(bener(1)*zfac)/1.d3
do i=2,bnene
  e1(i)=(bener(i)*zfac)/1.d3
  tau=moliv * bxs2(i)
  coeff(i)=dexp(-tau)
enddo

end subroutine absorption_olivine
! ======================================= !
subroutine map_to_grid_olivine(new_en,nne,old_en, one, nflux, old_flu,ifl)
! This routine maps to a given grid
implicit none
integer :: i, j, k, one, nne, bmin, bmax,ifl
double precision :: new_en(0:nne)
double precision :: old_en(0:one), old_flu(one)
double precision :: stemp,etemp, s, etemp2
real :: nflux(nne)
integer,parameter :: out_unit=20
do i=1,nne
nflux(i)=real(0.d0)
call dbinsrch_olivine(new_en(i-1),bmin,old_en,one+1)
call dbinsrch_olivine(new_en(i),bmax,old_en,one+1)
bmin = bmin-1
bmax = bmax-1
! Linear interpolation
if (bmin.eq.bmax) then
if(new_en(i).le.old_en(1))then
s=real(old_flu(1))
else if(new_en(i).gt.old_en(one))then
s=real(old_flu(one))
else
do j=2,one
if(new_en(i).gt.old_en(j-1).and.new_en(i).le.old_en(j))then
etemp2=(new_en(i)+new_en(i-1))/2
s=old_flu(j-1)+(old_flu(j)-old_flu(j-1))*(etemp2-old_en(j-1))/(old_en(j)-old_en(j-1))
endif
enddo
endif
! Average (integral)
else
stemp=0.d0
etemp=0.d0
do k=bmin,bmax
stemp=stemp+(old_flu(k))*(old_en(k)-old_en(k-1))
etemp=etemp+(old_en(k)-old_en(k-1))
enddo
s=real(stemp/etemp)
endif
nflux(i)=real(s)
enddo
end subroutine map_to_grid_olivine
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dbinsrch_olivine(e,k,ener,n)
!
! search for energy e in array ener(1:n) and return bin k for
! which ener(k).le.e.lt.ener(k+1)
! adapted from J. Wilms's tbvabs_new.f routine
!
implicit none
double precision :: e
integer :: n,k,klo,khi
double precision :: ener(n)
klo=1
khi=n
1 if (khi-klo.gt.1) then
k=(khi+klo)/2
if(ener(k).gt.e)then
khi=k
else
klo=k
endif
goto 1
endif
if (klo.eq.0) then
print *,'Energy out of bounds. Should not happen'
stop
endif
k=klo
end subroutine dbinsrch_olivine
