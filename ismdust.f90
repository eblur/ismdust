! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ISMDUST
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! XSPEC local model for dust absorption edge structure
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ismdust(ear, ne, param, ifl, photar)
!
! The main routine to call all subroutines
!
implicit none
integer,parameter :: num_param = 3
integer,parameter :: ngrain=1
integer,parameter :: nemod=25530 !Number of elements for each cross section.
integer :: ne, ifl, a
double precision :: msil, mgra, rshift, emod(nemod), coemod(nemod)
double precision :: bxs(0:ngrain,nemod), bener(nemod)
double precision :: zfac
real :: ear(0:ne), param(num_param), photar(ne)
logical :: startup=.true.
character (len=40) version
version='2.6'
 if(startup)then
  print *, ' '
  print *, 'ISMdust: High resolution XAFS model Version',version
  print *, 'Corrales, Garcia, Wilms, Nowak, & Baganoff (2015)'
  print *, 'Optical constants come from Draine 2003 (ApJ, 598, 1026)'
  print *, 'WARNING: If used in conjunction with neutral metal absorption models'
  print *, '(e.g. TBabs, TBnew), be sure to change abundances to stop'
  print *, 'from overestimating the ISM metal absorption component.'
  print *, ' '
  call read_cross_sections_ismdust(nemod,bxs,ifl,bener)
  startup=.false.
 endif
! Model parameters
msil = param(1)
mgra = param(2)
rshift = param(3)
zfac = 1.d0/(1.d0+dble(rshift))

call extinction_ismdust(msil, mgra, zfac, emod, nemod, coemod,bxs,ifl,bener)
!
call map_to_grid_ismdust(dble(ear),ne,emod,nemod,photar,coemod,ifl)
return
end subroutine ismdust
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_cross_sections_ismdust(bnene,xs,ifl,ener)
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
integer,parameter :: ngrain=1, out_unit=20
integer :: bnene, ifl, i, j, status
integer :: nemax
double precision :: ener(bnene), xs(0:ngrain,bnene)
character (*), parameter :: fileloc = '/edge_files/xs_ext_grid.fits'
character (*), parameter :: ismreadchat = 'ismdust: reading from '
character (len=255 + 29) :: filename2 ! len(fileloc)
character (len=240) :: local_dir = '.'
character (len=255) :: ismdust_root = ''
character (len=len(ismreadchat)+len(filename2)) :: chatmsg = ''
integer inunit,readwrite,blocksize
integer :: hdutype,colnum
integer :: felem=1, nulld=0
logical :: anynull
!Number of elements for each grain type cross section.
nemax=25530
! Where do we look for the data?
call getenv('ISMDUSTROOT', ismdust_root)
if (ismdust_root .EQ. '') then
ismdust_root = local_dir
print *, 'cannot find any ISMDUSTROOT environment'
endif
! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ismdust_root) // fileloc
chatmsg=ismreadchat // filename2
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

!Read in the cross section information
do i=0,ngrain
  colnum=i+2
  do j=1,nemax
    call ftgcvd(inunit,colnum,j,felem,1,nulld,xs(i,j),anynull,status)
    enddo
enddo

! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the X-Spec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used. An alternative would be to use xwrite()
! with a low chatter level.
!
! This message could be displayed only once, but it is probably worth
! repeating each time it is used.
if (status .ne. 0) then
write (*,*) 'ERROR: unable to read cross sections from ', filename2
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_cross_sections_ismdust
! ======================================= !
subroutine extinction_ismdust(msil, mgra, zfac, e1, bnene, coeff, bxs2,ifl,bener)
!
! This is routine that calculates the optical depth given the column densities
! Finally returns the absorption coefficient exp(-tau)
!
implicit none
integer,parameter :: ngrain=1, out_unit=20
integer :: bnene, ifl
integer :: i, j
double precision :: msil, mgra
double precision :: bener(bnene), bxs2(0:ngrain,bnene), e1(bnene)
double precision :: tau, coeff(bnene)
double precision :: zfac
real hphoto, gphoto
external hphoto, gphoto

! Calculates the optical depth and the extinction coefficient exp(-tau)
do i=1,bnene
  e1(i) = (bener(i)*zfac)/1.d3
  tau   = msil * bxs2(0,i) + mgra * bxs2(1,i)
  coeff(i) = dexp(-tau)
enddo

end subroutine extinction_ismdust
! ======================================= !
subroutine map_to_grid_ismdust(new_en,nne,old_en, one, nflux, old_flu,ifl)
! This routine maps to a given grid
implicit none
integer :: i, j, k, one, nne, bmin, bmax, ifl, btemp
double precision :: new_en(0:nne)
double precision :: old_en(one), old_flu(one)
double precision :: etemp, s
real :: nflux(nne)
integer,parameter :: out_unit=20
do i=1,nne
  nflux(i)=0.
  call dbinsrch_ismdust(new_en(i-1),bmin,old_en,one+1)
  call dbinsrch_ismdust(new_en(i),bmax,old_en,one+1)
  etemp = (new_en(i)+new_en(i-1))/2
  ! Linear interpolation
  if (bmin.eq.bmax) then
    if(new_en(i).le.old_en(0)) then
      s = real(old_flu(1))
    else if(new_en(i).gt.old_en(one)) then
      s = real(old_flu(one))
    else
      s = old_flu(bmin) + (old_flu(bmax+1)-old_flu(bmin)) * (etemp-old_en(bmin)) / (old_en(bmax+1)-old_en(bmin))
      endif
  else
    call dbinsrch_olivine(etemp,btemp,old_en,one+1)
    s = old_flu(btemp) + (old_flu(btemp+1)-old_flu(btemp)) * (etemp-old_en(btemp)) / (old_en(btemp+1)-old_en(btemp))
    endif
  nflux(i) = real(s)
  enddo
end subroutine map_to_grid_ismdust
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dbinsrch_ismdust(e,k,ener,n)
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
end subroutine dbinsrch_ismdust
