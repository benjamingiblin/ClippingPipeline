      module param
        real :: fscale
      end module

      program massrecon_new
      use param
! to run it for cosmos 
! ./massrecon_smoothKorig.Linux -n1 1638 -n2 1588 -tophat -scale 40 

! THE POINT OF THIS IS TO GAUSSIAN SMOOTH THE ORIGINAL KAPPA IMAGE AND THAT
! IS ALL. IT THINKS THE ORIGINAL KAPPA IMAGE IS A MASK. LOLz. 

!  ifort massrecon_smoothKorig.f90 -o massrecon_smoothKorig.a -L/usr/local/share/cfitsio/ -lcfitsio -L/roe/intel/Compiler/11.1/072/mkl/lib/em64t/ -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread -L/usr/lib -lm -lfftw3

!      parameter(n1=2048,n2=2048)
      real, dimension(:,:), allocatable :: rkappa,rphi,gamma1,gamma2,galdens,storegaldens,wgaussian
      real, dimension(:,:), allocatable :: rkappa_sto,gamma1_sto,gamma2_sto,rawmask,mask,smoothmask
      character*80 filename,kappaID(1),recon_type,pathname,maskflag
        character*150 arg
        character*150 opt
        integer iargc,narg,filtertest

      external sgaussian,saperture,stophat

      kappaID=(/'kappa93o.fits'/)
      recon_type='gaussian'
      n1=2048
      n2=2048
      fscale=10.
	pathname=''
	maskflag='yes'
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


         print *, 'USAGE: massrecon_smoothKorig.a '

         print *, '       Defaults:'
         STOP
             case ('-gaussian')
                recon_type='gaussian'
             case ('-aperture')
                recon_type='aperture'
             case ('-tophat')
                recon_type='tophat'
             case ('-nomask')
                maskflag='no'
             case ('-scale')
                read(arg,*) fscale
             case('-n1')
                read(arg,*) n1
             case('-n2')
                read(arg,*) n2
               case ('-path')
                read (arg,*) pathname     
     end select
        enddo
	write(*,*) recon_type,fscale,n1,n2
      allocate (rkappa(n1,n2),rphi(n1,n2),gamma1(n1,n2),gamma2(n1,n2),galdens(n1,n2),storegaldens(n1,n2),wgaussian(n1,n2))
      allocate (rkappa_sto(n1,n2),gamma1_sto(n1,n2),gamma2_sto(n1,n2),rawmask(n1,n2),mask(n1,n2),smoothmask(n1,n2))
      ibitpix=-32

!------------------------------------------------------kapparecon.fits

      filename='galdens.fits'
 !       write(*,*)pathname,filename
!stop
      call readimage(galdens,n1,n2,filename)
        storegaldens=galdens
!	write(*,*)'galdense read', filename
!stop
      filename='eiso1.fits'
      call readimage(gamma1_sto,n1,n2,filename)

      filename='eiso2.fits'
      call readimage(gamma2_sto,n1,n2,filename)
      if(maskflag.eq.'yes') then
        filename='./newmask.fits'
        call readimage(rawmask,n1,n2,filename)
!	where(rawmask.eq.0.) rawmask=-1.
!	where(rawmask.eq.1.) rawmask=-1.
!	where(rawmask.eq.32.) rawmask=-1.
!	where(rawmask.eq.256.) rawmask=-1.
!	where(rawmask.eq.1024.) rawmask=-1.
!	where(rawmask.eq.2048.) rawmask=-1.
!	where(rawmask.ge.0.) rawmask=0.
!	rawmask=-rawmask
        filename='mask.fits'
        call writeimage(rawmask,n1,n2,filename,ibitpix)
      endif
      if(maskflag.eq.'no') rawmask=1. 

      gamma1=gamma1_sto!*mask
      gamma2=gamma2_sto!*mask

	  
      select case (recon_type)
        case ('gaussian')
	  
	  smoothmask=rawmask
        flagconvo=0.
        xcenter=float(n1)/2.
        ycenter=float(n2)/2.
	    
        do i=1,n1
          do j=1,n2
            arg2=(float(i)-xcenter)**2+(float(j)-ycenter)**2
            wgaussian(i,j)=exp(-arg2/fscale**2*1.)*(1./3.14159/fscale**2)
          enddo
        enddo
		

          call convolfftw3(smoothmask,wgaussian,n1,n2,flagconvo)
		print *,'Smoothed the orig Kappa map successfully, I hope!'
        filename='smoothmask.fits'
        call writeimage(smoothmask,n1,n2,filename,ibitpix)
	
          call convolfftw3(gamma1,wgaussian,n1,n2,flagconvo)
          call convolfftw3(gamma2,wgaussian,n1,n2,flagconvo)
          call convolfftw3(galdens,wgaussian,n1,n2,flagconvo)
	where(galdens.gt.0.01) gamma1=gamma1/galdens
	where(galdens.gt.0.01) gamma2=gamma2/galdens

!         filename='g1presmooth.fits'
!          call writeimage(gamma1,n1,n2,filename,ibitpix)
!         filename='g2presmooth.fits'
!          call writeimage(gamma2,n1,n2,filename,ibitpix)
!         filename='galdenspresmooth.fits'
!          call writeimage(galdens,n1,n2,filename,ibitpix)

         filename='g1smooth.fits'
          call writeimage(gamma1,n1,n2,filename,ibitpix)
         filename='g2smooth.fits'
          call writeimage(gamma2,n1,n2,filename,ibitpix)
	  gamma1_sto=gamma1
	  gamma2_sto=gamma2

          call recon_kernel_g1(gamma1,n1,n2)
          call recon_kernel_g2(gamma2,n1,n2)
          rkappa=gamma1+gamma2

          call potential_kernel_g1(gamma1_sto,n1,n2)
          call potential_kernel_g2(gamma2_sto,n1,n2)
	  rphi=-(gamma1_sto+gamma2_sto)/2.
	
	where(smoothmask.le..5) rkappa=0.
          filename='galdenssmooth.fits'
          call writeimage(galdens,n1,n2,filename,ibitpix)
         filename='kapparenorm.fits'
          call writeimage(rkappa,n1,n2,filename,ibitpix)
         filename='potential.fits'
          call writeimage(rphi,n1,n2,filename,ibitpix)

        case ('aperture')
          call smooth_kernel(gamma1,n1,n2,saperture)
          call smooth_kernel(gamma2,n1,n2,saperture)
          call smooth_kernel(galdens,n1,n2,stophat)
	where(mask.eq.1.) galdens=0.
	where(mask.eq.1.) gamma1=0.
	where(mask.eq.0.) gamma1=gamma1/galdens
	where(mask.eq.1.) gamma2=0.
	where(mask.eq.0.) gamma2=gamma2/galdens

          call recon_kernel_g1(gamma1,n1,n2)
          call recon_kernel_g2(gamma2,n1,n2)
          rkappa=gamma1+gamma2

	where(mask.eq.1.) rkappa=0.
          filename='galdenssmooth.fits'
          call writeimage(galdens,n1,n2,filename,ibitpix)
!          filename='kapparecon.fits'
!          call writeimage(rkappa,n1,n2,filename,ibitpix)
          filename='kapparenorm.fits'
          call writeimage(rkappa,n1,n2,filename,ibitpix)


        case ('tophat')
	  smoothmask=rawmask

        flagconvo=0.
        xcenter=float(n1)/2.
        ycenter=float(n2)/2.
	wgaussian=0.
        do i=1,n1
          do j=1,n2
            arg2=(float(i)-xcenter)**2+(float(j)-ycenter)**2
            if(arg2.le.fscale) wgaussian(i,j)=1./(3.14159*fscale**2)
          enddo
        enddo

          call convolfftw3(smoothmask,wgaussian,n1,n2,flagconvo)

        filename='smoothmask.fits'
        call writeimage(smoothmask,n1,n2,filename,ibitpix)
	
          call convolfftw3(gamma1,wgaussian,n1,n2,flagconvo)
          call convolfftw3(gamma2,wgaussian,n1,n2,flagconvo)
          call convolfftw3(galdens,wgaussian,n1,n2,flagconvo)

         filename='g1presmooth.fits'
          call writeimage(gamma1,n1,n2,filename,ibitpix)
         filename='g2presmooth.fits'
          call writeimage(gamma2,n1,n2,filename,ibitpix)
         filename='galdenspresmooth.fits'
          call writeimage(galdens,n1,n2,filename,ibitpix)
	where(galdens.ne.0.) gamma1=gamma1/galdens
	where(galdens.ne.0.) gamma2=gamma2/galdens
         filename='g1smooth.fits'
          call writeimage(gamma1,n1,n2,filename,ibitpix)
         filename='g2smooth.fits'
          call writeimage(gamma2,n1,n2,filename,ibitpix)
	  gamma1_sto=gamma1
	  gamma2_sto=gamma2

          call recon_kernel_g1(gamma1,n1,n2)
          call recon_kernel_g2(gamma2,n1,n2)
          rkappa=gamma1+gamma2

          call potential_kernel_g1(gamma1_sto,n1,n2)
          call potential_kernel_g2(gamma2_sto,n1,n2)
	  rphi=-(gamma1_sto+gamma2_sto)/2.

	where(smoothmask.le..5) rkappa=0.
          filename='galdenssmooth.fits'
          call writeimage(galdens,n1,n2,filename,ibitpix)
         filename='kapparenorm.fits'
          call writeimage(rkappa,n1,n2,filename,ibitpix)
         filename='potential.fits'
          call writeimage(rphi,n1,n2,filename,ibitpix)

      end select

	stop

      end

        function sgaussian(f2)
        use param
!  1D profile of the smoothing window in Fourier space
!       gaussian
          sigma2=fscale*fscale
          sgaussian=exp(-f2*3.14159*3.14159*sigma2)
        end

        function saperture(f2)
        use param
!       aperture mass
          l=1
          f=sqrt(f2)*fscale
          saperture=2**l*float(l+3)/3.14159*bessj(3+l,f)/(f+1e-10)**(l+1)
        end

        function stophat(f2)
        use param
!       top-hat
          f=sqrt(f2)*fscale
          stophat=bessj1(f*(2.*3.14159))/(f+1e-10)/3.14159
       if(f.eq.0) stophat=1.
!!      if(f.lt.30.) write(*,*) f,sgaussian
        end

       include 'recon_kernel.f90'
       include 'potential_kernel.f90'
       include 'smooth_kernel.f90'
       include 'convolfftw3.f90'
       include 'cookbook.f90'

