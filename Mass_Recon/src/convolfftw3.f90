        subroutine convolfftw3(twoD_DATA_IN_no1,twoD_DATA_IN_no2,nx,ny,flagpsf)

!        real(8) :: twoD_DATA_IN_no1(nx,ny)
!        real(8) :: twoD_DATA_IN_no2(nx,ny)
        real :: twoD_DATA_IN_no1(nx,ny)
        real :: twoD_DATA_IN_no2(nx,ny)

	include '/usr/include/fftw3.f'

! FFTW parameters _ small array
	double precision array_in
	dimension array_in(1:nx,1:ny)
	double complex array_out
	dimension array_out(1:(nx/2+1),1:ny)

! FFTW parameters _ interpolated array
	double precision larray_in
	dimension larray_in(1:nx,1:ny)
!	double complex,save::larray_out
!	dimension larray_out(1:(nx/2+1),1:ny)
	double complex, dimension(:,:),allocatable, save :: larray_out

	double complex llarray_in
	dimension llarray_in(1:(nx/2+1),1:ny)
	double precision llarray_out
	dimension llarray_out(1:nx,1:ny)

	integer*8 :: plan,lplan,llplan
	integer len(2),flagpsf
	integer, save:: internalflagpsf

	if(internalflagpsf.ne.1) then
          allocate (larray_out(nx/2+1,ny))
        endif
	array_in=dble(twoD_DATA_IN_no1)
	larray_in=dble(twoD_DATA_IN_no2)

!	call dfftw_plan_dft_r2c_2d(plan,nx,ny,twoD_DATA_IN_no1,array_out,&
	call dfftw_plan_dft_r2c_2d(plan,nx,ny,array_in,array_out,&
                        FFTW_ESTIMATE)
!	if(internalflagpsf.ne.1) call dfftw_plan_dft_r2c_2d(lplan,nx,ny,twoD_DATA_IN_no2,larray_out,&
	if(internalflagpsf.ne.1) call dfftw_plan_dft_r2c_2d(lplan,nx,ny,larray_in,larray_out,&
                        FFTW_ESTIMATE)
!	call dfftw_plan_dft_c2r_2d(llplan,nx,ny,llarray_in,twoD_DATA_IN_no1,&
	call dfftw_plan_dft_c2r_2d(llplan,nx,ny,llarray_in,llarray_out,&
                        FFTW_ESTIMATE)

	call dfftw_execute(plan)
	if(internalflagpsf.ne.1) call dfftw_execute(lplan)
!	call dfftw_execute(lplan)
	llarray_in=larray_out*array_out
	call dfftw_execute(llplan)
	llarray_out=llarray_out/float(nx*ny)

	twoD_DATA_IN_no1=real(cshift(array=cshift(array=llarray_out,shift=nx/2,dim=1),shift=ny/2,dim=2))

!	write(*,*) internalflagpsf,flagpsf,nx,ny
	if(flagpsf.ne.1) then
	  flagpsf=1
	  internalflagpsf=1
	endif

        end

