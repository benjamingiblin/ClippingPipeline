	program cat2grid_frombin
        real, dimension(:), allocatable :: x1,x2,e1,e2,weight,rm
        real, dimension(:), allocatable :: x1boots,x2boots,e1boots,e2boots,weightboots,rmboots
        real, dimension(:,:), allocatable :: data
        real, dimension(:,:), allocatable :: galdens,eiso1,eiso2,calib
        character*80 filein,filterstring,sweight,efield,xfield,extrafield,fileout,setw
        character*150 arg,filename,rotate
	character*1 boots
        character*150 opt
        integer iargc,narg,filtertest,igraine,igraineboots
        real*8 :: xmin,xmax,ymin,ymax

! ifort cat2grid_frombin.f90 -o cat2grid_frombin.a -L/usr/local/share/cfitsio/ -lcfitsio

	setw='yes'
        rmfileout=.false.
        xfield="x = %x"
        efield="e = %eiso"
        extrafield=''
        filtertest=0
        nx=2048
        ny=2048
	xminin=0.
	xmaxin=21000.
	yminin=0.
	ymaxin=21000.
	igraine=0
	rotate="n"
	igraineboots=-1
	boots='n'
	zbin=0.

987     format(a112)
988     format(a101)

        narg=iargc()

        do i=1,narg
           call getarg(i,opt)
           call getarg(i+1,arg)
           select case (opt)
             case ('-delete')
                rmfileout=.true.
           end select
        enddo

        do i=1,narg
           call getarg(i,opt)
           call getarg(i+1,arg)
           select case (opt)
             case ('-h')
         print *,''

         print 987, 'USAGE: cat2grid.Linux'
         print *,'Format: x,y,e1,e2,weight'
         print *,''
         STOP
             case ('-in')
                filein=arg
             case ('-boots')
                read(arg,*) igraineboots
                boots='y'
             case('-nx')
                read(arg,*) nx
             case('-ny')
                read(arg,*) ny
             case('-xmin')
                read(arg,*) xminin
             case('-xmax')
                read(arg,*) xmaxin
             case('-ymin')
                read(arg,*) yminin
             case('-ymax')
                read(arg,*) ymaxin
             case ('-rot')
                rotate="y"
             case ('-ran')
                read(arg,*) igraine
             case ('-noweight')
                setw='no'
           end select
        enddo

	write(*,*) "read input catalogue...", filein
	open(12,file=filein,form='unformatted')
          read(12) nligne,ncol
	  allocate(data(nligne,ncol))
          read(12) data
        close(12)
	ngal=nligne
	write(*,*) nligne,ncol
        allocate (x1(nligne),x2(nligne),e1(nligne),e2(nligne),weight(nligne),rm(nligne))
        allocate (x1boots(nligne),x2boots(nligne),e1boots(nligne),e2boots(nligne),weightboots(nligne),rmboots(nligne))
	write(*,*) 'filter catalogue...'

!  this is for the older catalogue version before e2 correction has been fixed
!	  xpos=data(i,1)
!	  ypos=data(i,2)
!	  wcs1=data(i,3)
!	  wcs2=data(i,4)
!	  mask=data(i,5)
!	  weight=data(i,6)
!	  e1=data(i,7)
!	  e2=data(i,8)
!	  zl=data(i,9)
!	  type=data(i,10)
!	  odds=data(i,11)
!	  fluxradius=data(i,12)
!	  rmag=data(i,13)
!	  starclass=data(i,14)
!	  starmasslog=data(i,15)

!  this is for the newer catalogue version after e2 correction has been fixed
!	  xpos=data(i,1)
!	  ypos=data(i,2)
!	  wcs1=data(i,3)
!	  wcs2=data(i,4)
!	  mask=data(i,5)
!	  weight=data(i,6)
!	  e1=data(i,7)
!	  e2=data(i,8)
!	  c2=data(i,9)   !  for the simulated catalogues this is equal to the noisefree kappa value
!	  rm=data(i,10)
!	  zl=data(i,11)
!	  type=data(i,12)
!	  odds=data(i,13)
!	  fluxradius=data(i,14)
!	  rmag=data(i,15)
!	  starclass=data(i,16)
!	  starflag=data(i,17)
!	  starmasslog=data(i,18)
!	  flagexlcude=data(i,19)

!  selection is done here
	icpt=1
	do i=1,nligne
!	write(*,*) data(i,:)
!  old catlaogue version
!	if(data(i,6).gt.0..and.data(i,5).le.1..and.data(i,5).ge.0..and.data(i,9).ge.0.5.and.data(i,9).le.1.2) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.0.and.data(i,9).le.10.) then ! set 1
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.5.and.data(i,9).le.1.2) then ! set 2
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 3
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.2.and.data(i,9).le.0.9) then ! set 4
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,5).ge.0..and.data(i,13).ge.21.5.and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 5
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,5).ge.0..and.data(i,13).le.24..and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 6
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,13).le.24.5.and.data(i,13).gt.21.5.and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 7
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,13).le.24.5.and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 8
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,13).le.24.0.and.data(i,13).gt.21.5.and.data(i,9).ge.0.0.and.data(i,9).le.1.2) then ! set 9
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,5).ge.0..and.data(i,13).le.24..and.data(i,9).ge.0.0.and.data(i,9).le.4.2) then ! set 10
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,5).ge.0..and.data(i,13).le.24..and.data(i,9).ge.0.0.and.data(i,9).le.1.0) then ! set 11
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,5).ge.0..and.data(i,13).le.24..and.data(i,9).ge.0.5.and.data(i,9).le.1.2) then ! set 12
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.1.and.data(i,9).le.1.2) then ! set 13
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.3.and.data(i,9).le.1.2) then ! set 14
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.3.and.data(i,9).le.1.1) then ! set 15
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.3.and.data(i,9).le.1.0) then ! set 16
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,9).ge.0.4.and.data(i,9).le.1.1) then ! set 17 this was the last used
!	if(data(i,6).gt.0..and.data(i,5).lt.2..and.data(i,9).ge.0.1.and.data(i,9).le.1.3) then
!  new catalogue version
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.1.and.data(i,11).le.1.3.and.data(i,19).eq.0.) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.1.and.data(i,11).le.1.1.and.data(i,19).eq.0.) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.1.and.data(i,11).le.1.1) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.4.and.data(i,11).le.1.1) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.4.and.data(i,11).le.1.1) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.4.and.data(i,11).le.1.1.and.data(i,19).eq.0.) then
!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.1.and.data(i,11).le.1.1.and.data(i,19).eq.0.) then

	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,11).ge.0.4.and.data(i,11).le.1.1.and.data(i,19).eq.0.) then  !  data and most tests have been run with this

!	if(data(i,6).gt.0..and.data(i,5).lt.1..and.data(i,19).eq.0.) then  !  used for simul (redshift cut out)

	iz=1
	do j=1,Nz
	  if(data(i,11).ge.zslice(j)) iz=j
	enddo
	zbin(iz)=zbin(iz)+1.

	  x1(icpt)=data(i,1)
	  x2(icpt)=data(i,2)
	  e1(icpt)=data(i,7)
	  e2(icpt)=data(i,8)
	  weight(icpt)=data(i,6)
	  rm(icpt)=data(i,10)
	    if(rotate.eq."y") then
	      esave=e1(icpt)
	      e1(icpt)=-e2(icpt)
	      e2(icpt)=esave
	    endif
            if(igraine.ne.0) then
              e=sqrt(e1(icpt)**2+e2(icpt)**2)
              rantheta=ran2(igraine)*3.14159
              e1(icpt)=e*cos(2.*rantheta)
              e2(icpt)=e*sin(2.*rantheta)
            endif
!	write(*,*) x1(icpt),x2(icpt),e1(icpt),e2(icpt),weight(icpt)
	  icpt=icpt+1
	endif
	enddo
	icptmax=icpt-1
	write(*,*) icptmax





	allocate (galdens(nx,ny),eiso1(nx,ny),eiso2(nx,ny),calib(nx,ny))
!	xmin=minval(x1)
!	xmax=maxval(x1)
!	ymin=minval(x2)
!	ymax=maxval(x2)
	xmin=xminin
	xmax=xmaxin
	ymin=yminin
	ymax=ymaxin
!       write(*,*)xmin,ymin
!
!stop

	x1=x1-xmin
	x2=x2-ymin
	dx=xmax/float(nx)
	dy=ymax/float(ny)
!	write(*,*) int(8000./dx)+1,int(8000./dy)+1,dx,dy

!  generates bootstrap sample if -boots keyword is present
        if(boots.eq.'y') then
          write(*,*) 'start bootstrap with seed ',igraineboots
          x1boots=x1
          x2boots=x2
          e1boots=e1
          e2boots=e2
          weightboots=weight
          rmboots=rm
          do igalid=1,icptmax
!       write(*,*) int(ran2(igraine)*float(ngalmax-1))
            igalpick=1+int(ran2(igraineboots)*float(icptmax-1))
!	write(*,*) igalid,igalpick
            x1boots(igalid)=x1(igalpick)
            x2boots(igalid)=x2(igalpick)
            e1boots(igalid)=e1(igalpick)
            e2boots(igalid)=e2(igalpick)
            weightboots(igalid)=weight(igalpick)
            rmboots(igalid)=rm(igalpick)
!	    e=sqrt(e1(igalid)**2+e2(igalid)**2)
!            e1boots(igalid)=e*cos(2.*3.14159*ran2(igraineboots))
!            e2boots(igalid)=e*sin(2.*3.14159*ran2(igraineboots))
!	    x1boots(igalid)=x1(igalid)
!	    x2boots(igalid)=x2(igalid)
!	    weightboots(igalid)=weight(igalid)
          enddo
          x1=x1boots
          x2=x2boots
          e1=e1boots
          e2=e2boots
          weight=weightboots
          rm=rmboots
          write(*,*) 'bootstrap done'
        endif   



	galdens=0.
	eiso1=0.
	eiso2=0.
	calib=0.
	i=1
	if(setw.eq.'yes') then
	do i=1,icptmax
	  index1=maxval((/int(x1(i)/dx)+1,1/))
	  index1=minval((/int(x1(i)/dx)+1,nx/))
	  index2=maxval((/int(x2(i)/dy)+1,1/))
	  index2=minval((/int(x2(i)/dy)+1,ny/))
	  galdens(index1,index2)=galdens(index1,index2)+weight(i)*(1.+rm(i))
	  eiso1(index1,index2)=eiso1(index1,index2)+e1(i)*weight(i)
	  eiso2(index1,index2)=eiso2(index1,index2)+e2(i)*weight(i)
	  calib(index1,index2)=calib(index1,index2)+weight(i)*(1.+rm(i))
	enddo
	endif
	if(setw.eq.'no') then
	do i=1,icptmax
	  index1=maxval((/int(x1(i)/dx)+1,1/))
	  index1=minval((/int(x1(i)/dx)+1,nx/))
	  index2=maxval((/int(x2(i)/dy)+1,1/))
	  index2=minval((/int(x2(i)/dy)+1,ny/))
	  galdens(index1,index2)=galdens(index1,index2)+1.
	  eiso1(index1,index2)=eiso1(index1,index2)+e1(i)
	  eiso2(index1,index2)=eiso2(index1,index2)+e2(i)
	  calib(index1,index2)=calib(index1,index2)+1.
	enddo
	endif

	open(11,file='zbin.dat')
	  do i=1,Nz
	    write(11,*) redshift(i),zbin(i)
	  enddo
	close(11)

	ibitpix=-32
	filename='galdens.fits'
	call writeimage(galdens,nx,ny,filename,ibitpix)
	filename='eiso1.fits'
	call writeimage(eiso1,nx,ny,filename,ibitpix)
	filename='eiso2.fits'
	call writeimage(eiso2,nx,ny,filename,ibitpix)
	filename='calib.fits'
	call writeimage(calib,nx,ny,filename,ibitpix)

	write(*,*) "...Done."
	write(*,*) 'Now run mask2grid with options:'
	write(*,888) 'mask2grid.Linux -dx ',dx,' -dy ',dy,' -nx ',nx,' -ny ',ny,' -xmin ',xmin,' -ymin ',ymin
888	format(a20,f6.2,a5,f6.2,a5,i4,a5,i4,a7,f10.2,a7,f10.2)
!	call system('mask2grid.Linux -h')

	end

        include '../numrec/ran2.for'
        include '/home/cech/cookbook.f90'
