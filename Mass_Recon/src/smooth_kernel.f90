        subroutine smooth_kernel(twoD_DATA_IN,nx,ny,func)

!  smooth_kernel is the Fourier convolution kernel. twoD_DATA_IN is the 2D array
!  input data to be smoothed and func if the 1D fourier transform of the
!  convolution filter. func HAS TO BE REAL.

!  ifort -Vaxlib -w smooth_kernel.f90 -o smooth_kernel.Linux -L/home/waerbeke/cfitsio/ -lcfitsio -L/usr/opt/intel_fc_80/mkl72/lib/32 -lmkl_ia32 -lguide -lpthread -lm
!  ./smooth_kernel.Linux

        use MKL_DFTI

        real :: twoD_DATA_IN(nx,ny)
!        real, intent(in), dimension(:,:) :: twoD_DATA_IN

        real, dimension(:,:), allocatable :: Y_IN_2D
        real, dimension(:), allocatable :: Y_IN
        real, dimension(:,:), allocatable :: Y_OUT_2D
        real, dimension(:), allocatable :: Y_OUT
        real k1,k2

        external func

        type (DFTI_DESCRIPTOR), POINTER :: my_descl_handle
        integer status,L(2)
        integer strides_in(3)
        integer strides_out(3)

!  allocate arrays from input and output arrays
!        n1=2048 !nx
!        n2=2048 !ny
        n1=nx
        n2=ny
        allocate (Y_IN_2D(n1,n2))
        allocate (Y_IN(n1*n2))
        allocate (Y_OUT_2D(n1+2,n2+2))
        allocate (Y_OUT((n1+2)*(n2+2)))

!  load input arrays for Fourier transform
        Y_IN_2D=twoD_DATA_IN
	Y_IN=reshape(Y_IN_2D,(/n1*n2/))

        L(1)=n1
        L(2)=n2

        strides_in(1) = 0
        strides_in(2) = 1
        strides_in(3) = n1

        strides_out(1) = 0
        strides_out(2) = 1
        strides_out(3) = n1+2

        status=DftiCreateDescriptor(my_descl_handle,DFTI_SINGLE,DFTI_REAL,2,L)
        Status = DftiSetValue(my_descl_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        Status = DftiSetValue(my_descl_handle, DFTI_INPUT_STRIDES, strides_in)
        Status = DftiSetValue(my_descl_handle, DFTI_OUTPUT_STRIDES, strides_out)

        status=DftiCommitDescriptor(my_descl_handle)
        status=DftiComputeForward(my_descl_handle,Y_IN,Y_OUT)
        Y_OUT_2D=reshape(Y_OUT,(/n1+2,n2+2/))

!  smoothing in Fourier space ; f2 is the square wavenumber module
!  the wavenumber is (k1,k2)
        do j=1,n2/2+1
          k1=0.
          if(j.eq.1.) k2=0
          if(j.gt.1..and.j.lt.(n2/2+1)) k2=float(j-1)/float(n2)
          if(j.eq.(n2/2+1)) k2=1./2.
          f2=k1*k1+k2*k2
          Wsmooth=func(f2)
          Y_OUT_2D(1,2*j-1)=Y_OUT_2D(1,2*j-1)*Wsmooth
          Y_OUT_2D(1,2*j)=Y_OUT_2D(1,2*j)*Wsmooth
        enddo
        do i=2,n1/2
          do j=1,n2
            k1=float(i-1)/float(n1)
            if(j.eq.1) k2=0.
            if(j.gt.1.and.j.lt.(n2/2+1)) k2=float(j-1)/float(n2)
            if(j.eq.(n2/2+1)) k2=1./2.
            if(j.gt.(n2/2+1)) k2=float(j-n2-1)/float(n2)
            f2=k1*k1+k2*k2
            Wsmooth=func(f2)
            Y_OUT_2D(2*i-1,j)=Y_OUT_2D(2*i-1,j)*Wsmooth
            Y_OUT_2D(2*i,j)=Y_OUT_2D(2*i,j)*Wsmooth
          enddo
        enddo
        do j=1,n2/2+1
          k1=1./2.
          if(j.eq.1.) k2=0
          if(j.ne.1..and.j.ne.(n2/2+1)) k2=float(j-1)/float(n2)
          if(j.eq.(n2/2+1)) k2=1./2.
          f2=k1*k1+k2*k2
          Wsmooth=func(f2)
          Y_OUT_2D(n1+1,2*j-1)=Y_OUT_2D(n1+1,2*j-1)*Wsmooth
          Y_OUT_2D(n1+1,2*j)=Y_OUT_2D(n1+1,2*j)*Wsmooth
        enddo
        Y_OUT=reshape(Y_OUT_2D,(/(n1+2)*(n2+2)/))

        Status = DftiSetValue(my_descl_handle, DFTI_INPUT_STRIDES, strides_out)
        Status = DftiSetValue(my_descl_handle, DFTI_OUTPUT_STRIDES, strides_in)
        Status = DftiCommitDescriptor(my_descl_handle)
        Scale = 1.0/real(n1*n2, KIND=4)
!        Status = DftiSetValue(my_descl_handle, DFTI_BACKWARD_SCALE, Scale)
        Status = DftiComputeBackward(my_descl_handle, Y_OUT, Y_IN)
        status=DftiFreeDescriptor(my_descl_handle)

        Y_IN_2D=reshape(Y_IN,(/n1,n2/))
        twoD_DATA_IN=Y_IN_2D*Scale

        end

      !  include '/home/cech/cookbook.f90'
        include '../numrec/gamma.for'
        include '../numrec/bessj.for'
        include '../numrec/bessj0.for'
        include '../numrec/bessi1.for'
        include '../numrec/bessj1.for'
        include '../numrec/factrl.for'

