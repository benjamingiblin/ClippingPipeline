        program cat2grid_frombin
!integer, parameter :: dp = selected_real_kind(8, 307)
          real, dimension(:), allocatable :: x1,x2,e1,e2,weight,rm
        real, dimension(:), allocatable :: x1boots,x2boots,e1boots,e2boots,weightboots,rmboots
        real, dimension(:,:), allocatable :: data
        real, dimension(:,:), allocatable :: galdens,eiso1,eiso2,calib
        character*200 filein,LOS,filterstring,sweight,efield,xfield,extrafield,fileout,setw
        character*200 arg,rotate
        character*500 inoutdir,filename
	character*1 boots
        character*200 opt
        integer iargc,narg,filtertest,igraine,igraineboots
        real*8 :: xmin,xmax,ymin,ymax

! ifort cat2grid_fromasc.f90 -o cat2grid_fromasc.a -L/usr/local/share/cfitsio/ -lcfitsio -O0 -CB

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
			 case ('-directory')
				read(arg,'(a)') inoutdir	
			 case ('-LOS')
				read(arg,*) LOS	
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
        open(12,file=filein,status='old')
        read(12,*) nligne,ncol
        allocate(data(nligne,ncol))
        do i = 1,nligne
           read(12,*) (data(i,j),j=1,ncol)
        end do
        close(12)
        
	ngal=nligne
	write(*,*) 'Analysing this many galaxies', nligne
        allocate (x1(nligne),x2(nligne),e1(nligne),e2(nligne),weight(nligne),rm(nligne))
        allocate (x1boots(nligne),x2boots(nligne),e1boots(nligne),e2boots(nligne),weightboots(nligne),rmboots(nligne))
	write(*,*) 'filter catalogue...'

!  selection is done here
	icpt=1
	do i=1,nligne
	  x1(icpt)=data(i,1)
	  x2(icpt)=data(i,2)
	  e1(icpt)=data(i,3)
	  e2(icpt)=data(i,4)
	  weight(icpt)=data(i,5)
	 ! rm(icpt)=data(i,10)
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
	  icpt=icpt+1
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
!        write(*,*) 'xmax is ',xmax
!        write(*,*) 'ymax is ',ymax
!        write(*,*) 'xmin is ',xmin
!        write(*,*) 'ymin is ',ymin

	x1=x1-xmin
	x2=x2-ymin
	dx=(xmax-xmin)/float(nx)
 dy=(ymax-ymin)/float(ny)
!write(*,*) 'ymin, ymax, (ymax-ymin), float(ny) and (ymax-ymin)/float(ny) are: ',ymin,ymax,(ymax-ymin),float(ny),dy
 
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
    ! 25/05/2019 - Giblin, these lines had a bug in them that caused a seg fault
    ! in cases where x1,x2 were too close to the edge. The 1st arg in the minval() 
    ! lines matched the 1st arg in the maxval() lines. It needs to be index1/2
    do i=1,icptmax
	  index1=maxval((/int(x1(i)/dx)+1,1/))
	  index1=minval((/index1,nx/))
	  index2=maxval((/int(x2(i)/dy)+1,1/))
   index2=minval((/index2,ny/))
   ! Confirmed with this statement, seg-fault comes from y_min too close to low-limit
   ! i.e. index2 comes out as 0, which is <1.
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

        ibitpix=-32
	filename=trim(inoutdir) // 'galdens' // trim(LOS) // '.fits'
	call writeimage(galdens,nx,ny,filename,ibitpix)
	filename=trim(inoutdir) // 'eiso1' // trim(LOS) // '.fits'
	call writeimage(eiso1,nx,ny,filename,ibitpix)
	filename=trim(inoutdir) // 'eiso2' // trim(LOS) // '.fits'
	call writeimage(eiso2,nx,ny,filename,ibitpix)
	filename=trim(inoutdir) // 'calib' // trim(LOS) // '.fits'
	call writeimage(calib,nx,ny,filename,ibitpix)

	write(*,*) "...Done."
!	write(*,*) 'Now run mask2grid with options:'
!	write(*,888) 'mask2grid.Linux -dx ',dx,' -dy ',dy,' -nx ',nx,' -ny ',ny,' -xmin ',xmin,' -ymin ',ymin
!888	format(a20,f6.2,a5,f6.2,a5,i4,a5,i4,a7,f10.2,a7,f10.2)
!	call system('mask2grid.Linux -h')

	end

        include '../numrec/ran2.for'
        include 'cookbook.f90'
