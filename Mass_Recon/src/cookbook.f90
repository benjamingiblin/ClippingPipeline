     subroutine writeimage(data,n1,n2,filename,ibitpix)

!     Create a FITS primary array containing a 2-D image

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer i,j,group,fpixel,nelements
	real data(n1,n2)
      character filename*500
      logical simple,extend

 1    status=0

!     Delete the file if it already exists, so we can then recreate it
 2    call deletefile(filename,status)

!     Get an unused Logical Unit Number to use to open the FITS file
 3    call ftgiou(unit,status)

!     create the new empty FITS file
      blocksize=1
 4    call ftinit(unit,filename,blocksize,status)

!     initialize parameters about the FITS image (300 x 200 16-bit integers)
      simple=.true.
      bitpix=ibitpix
      naxis=2
      naxes(1)=n1
      naxes(2)=n2
      extend=.true.

!     write the required header keywords
 5    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!     initialize the values in the image with a linear ramp function
!     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)

	if(ibitpix.eq.8) call ftpprb(unit,group,fpixel,nelements,int(data,8),status)
	if(ibitpix.eq.16) call ftppri(unit,group,fpixel,nelements,int(data),status)
	if(ibitpix.eq.32) call ftpprj(unit,group,fpixel,nelements,int(data),status)
	if(ibitpix.eq.-32) call ftppre(unit,group,fpixel,nelements,data,status)
	if(ibitpix.eq.-64) call ftpprd(unit,group,fpixel,nelements,data,status)

!     close the file and free the unit number
 8    call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 9    if (status .gt. 0)call printerror(status)
      end


      subroutine writeascii

!     Create an ASCII table containing 3 columns and 6 rows

      integer status,unit,readwrite,blocksize,tfields,nrows,rowlen
      integer nspace,tbcol(3),diameter(6), colnum,frow,felem
      real density(6)
      character filename*40,extname*16
      character*16 ttype(3),tform(3),tunit(3),name(6)
      data ttype/'Name','Diameter','Density'/
      data tform/'A8','I6','F4.2'/
      data tunit/' ','km','g/cm'/
      data name/'Mercury','Venus','Earth','Mars','Jupiter','Saturn'/
      data diameter/4880,12112,12742,6800,143000,121000/
      data density/5.1,5.3,5.52,3.94,1.33,0.69/

 1    status=0

!     Get an unused Logical Unit Number to use to open the FITS file
 2    call ftgiou(unit,status)

!     open the FITS file, with write access
 3    readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)

!     append a new empty extension onto the end of the primary array
 4    call ftcrhd(unit,status)

!     define parameters for the ASCII table (see the above data statements)
      tfields=3
      nrows=6
      extname='PLANETS_ASCII'
      
!     calculate the starting position of each column, and the total row length
      nspace=1
 5    call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)

!     write the required header parameters for the ASCII table
 6    call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit,extname,status)

!     write names to the first column, diameters to 2nd col., and density to 3rd
      frow=1
      felem=1
      colnum=1
 7    call ftpcls(unit,colnum,frow,felem,nrows,name,status)
      colnum=2
      call ftpclj(unit,colnum,frow,felem,nrows,diameter,status)  
      colnum=3
      call ftpcle(unit,colnum,frow,felem,nrows,density,status)  

!     close the FITS file and free the unit number
 8    call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 9    if (status .gt. 0)call printerror(status)
      end


      subroutine writebintable

!     Create a binary table containing 3 columns and 6 rows

      integer status,unit,readwrite,blocksize,hdutype,tfields,nrows
      integer varidat,diameter(6), colnum,frow,felem
      real density(6)
      character filename*40,extname*16
      character*16 ttype(3),tform(3),tunit(3),name(6)
      data ttype/'Name','Diameter','Density'/
      data tform/'8A','1J','1E'/
      data tunit/' ','km','g/cm'/
      data name/'Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'/
      data diameter/6800,143000,121000,47000,45000,6000/
      data density/3.94,1.33,0.69,1.56,2.27,1.0/

 1    status=0
!     Name of the FITS file to append the ASCII table to:
      filename='ATESTFILEZ.FITS'

!     Get an unused Logical Unit Number to use to open the FITS file
 2    call ftgiou(unit,status)

!     open the FITS file, with write access
 3    readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)

!     move to the last (2nd) HDU in the file
 4    call ftmahd(unit,2,hdutype,status)

!     append/create a new empty HDU onto the end of the file and move to it
 5    call ftcrhd(unit,status)

!     define parameters for the binary table (see the above data statements)
      tfields=3
      nrows=6
      extname='PLANETS_BINARY'
      varidat=0
      
!     write the required header parameters for the binary table
 6    call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

!     write names to the first column, diameters to 2nd col., and density to 3rd
      frow=1
      felem=1
      colnum=1
 7    call ftpcls(unit,colnum,frow,felem,nrows,name,status)
      colnum=2
      call ftpclj(unit,colnum,frow,felem,nrows,diameter,status)  
      colnum=3
      call ftpcle(unit,colnum,frow,felem,nrows,density,status)  

!     close the FITS file and free the unit number
 8    call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 9    if (status .gt. 0)call printerror(status)
      end



      subroutine copyhdu(infilename)

!     copy the 1st and 3rd HDUs from the input file to a new FITS file

      integer status,inunit,outunit,readwrite,blocksize,morekeys,hdutype
      character infilename*500,outfilename*500

 1    status=0
!     Name of the FITS files:
!      infilename='ATESTFILEZ.FITS'
      outfilename='BTESTFILEZ.FITS'

!     Delete the file if it already exists, so we can then recreate it
 2    call deletefile(outfilename,status)

!     Get  unused Logical Unit Numbers to use to open the FITS files
 3    call ftgiou(inunit,status)
      call ftgiou(outunit,status)

!     open the input FITS file, with readonly access
      readwrite=0
 4    call ftopen(inunit,infilename,readwrite,blocksize,status)

!     create the new empty FITS file with the standard block size
      blocksize=1
 5    call ftinit(outunit,outfilename,blocksize,status)

!     copy the primary array from the input file to the output file
      morekeys=0
 6    call ftcopy(inunit,outunit,morekeys,status)

!     append/create a new empty extension on the end of the output file
! 7    call ftcrhd(outunit,status)

!     skip to the 3rd extension in the input file
! 8    call ftmahd(inunit,3,hdutype,status)

!     copy this extension from the input file to the output file
! 9    call ftcopy(inunit,outunit,morekeys,status)  

!     close the FITS file and free the unit numbers
 10   call ftclos(inunit, status)
      call ftclos(outunit, status)
 11   call ftfiou(-1, status)

!     check for any error, and if so print out error messages
 12   if (status .gt. 0)call printerror(status)
      end



      subroutine selectrows

!     select rows from an input table and copy them to the output table

      integer status,inunit,outunit,readwrite,blocksize,hdutype
      integer nkeys,nspace,naxes(2),nfound,colnum,frow,felem
      integer noutrows,irow,temp(100),i
      real nullval,density(6)
      character infilename*40,outfilename*40,record*80
      logical exact,anynulls

 1    status=0
!     Names of the FITS files:
      infilename='ATESTFILEZ.FITS'
      outfilename='BTESTFILEZ.FITS'

!     Get  unused Logical Unit Numbers to use to open the FITS files
 2    call ftgiou(inunit,status)
      call ftgiou(outunit,status)

!     open the FITS files, with the appropriate read/write access
      readwrite=0
 3    call ftopen(inunit,infilename,readwrite,blocksize,status)
      readwrite=1
      call ftopen(outunit,outfilename,readwrite,blocksize,status)

!     move to the 3rd HDU in the input file (a binary table in this case)
 4    call ftmahd(inunit,3,hdutype,status)

!     move to the last extension in the output file
 5    do while (status .eq. 0)
          call ftmrhd(outunit,1,hdutype,status)
      end do

      if (status .eq. 107)then
!         this is normal; it just means we hit the end of file
          status=0
 6        call ftcmsg
      end if

!     create a new empty extension in the output file
 7    call ftcrhd(outunit,status)

!     find the number of keywords in the input table header
 8    call ftghsp(inunit,nkeys,nspace,status)

!     copy all the keywords from the input to the output extension
 9    do i=1,nkeys
          call ftgrec(inunit,i,record,status)
          call ftprec(outunit,record,status)
      end do

!     force FITSIO to read the output file keywords to define the data structure
 10   call ftrdef(outunit,status)
!     get the width of the table (in bytes) and the number of rows
 11   call ftgknj(inunit,'NAXIS',1,2,naxes,nfound,status)

!     find which column contains the DENSITY values
      exact=.false.
 12   call ftgcno(inunit,exact,'DENSITY',colnum,status)

!     read the DENSITY column values
      frow=1
      felem=1
      nullval=-99.
 13   call ftgcve(inunit,colnum,frow,felem,naxes(2),nullval,density,anynulls,status)

!     If the density is less than 3.0, copy the row to the output table
      noutrows=0
 14   do irow=1,naxes(2)
          if (density(irow) .lt. 3.0)then
              noutrows=noutrows+1
 15           call ftgtbb(inunit,irow,1,naxes(1),temp,status)
              call ftptbb(outunit,noutrows,1,naxes(1),temp,status)
          end if
      end do

!     update the NAXIS2 keyword with the correct no. of rows in the output file
 16   call ftmkyj(outunit,'NAXIS2',noutrows,'&',status)

!     close the FITS file and free the unit numbers
 17   call ftclos(inunit, status)
      call ftclos(outunit, status)
      call ftfiou(-1, status)

!     check for any error, and if so print out error messages
 18   if (status .gt. 0)call printerror(status)
      end




      subroutine readheader(filename)

!     Print out all the header keywords in all extensions of a FITS file

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i
      character filename*500,record*80

 1    status=0

!     Get an unused Logical Unit Number to use to open the FITS file
 2    call ftgiou(unit,status)

!     name of FITS file 
!      filename='ATESTFILEZ.FITS'

!     open the FITS file, with read-only access
      readwrite=0
 3    call ftopen(unit,filename,readwrite,blocksize,status)

100   continue

!     Determine the number of keywords in the header
 4    call ftghsp(unit,nkeys,nspace,status)

!     Read each 80-character keyword record, and print it out
      do i = 1, nkeys
 5        call ftgrec(unit,i,record,status)
          print *,record
      end do

!     Print out and END record, and a blank line to mark the end of the header
      if (status .eq. 0)then
          print *,'END'
          print *,' '
      end if

!     try moving to the next extension in the FITS file, if it exists
 6    call ftmrhd(unit,1,hdutype,status)

      if (status .eq. 0)then
!         success, so loop back and print out keywords in this extension
 7        go to 100

      else if (status .eq. 107)then
!         hit end of file, so quit
 8        print *,'***** END OF FILE *****'
          status=0
          call ftcmsg
      end if

!     close the file, free the unit number, and exit
 9    call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 10   if (status .gt. 0)call printerror(status)
      end

      subroutine readimage(data_lin,n1,n2,filename)

!     Read a FITS image and determine the minimum and maximum pixel value

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real datamin,datamax,nullval,buffer(100)
	real data_lin(n1*n2)
      logical anynull
      character filename*500

 1    status=0

!     Get an unused Logical Unit Number to use to open the FITS file
 2    call ftgiou(unit,status)

!     open the FITS file previously created by WRITEIMAGE
      readwrite=0
 3    call ftopen(unit,filename,readwrite,blocksize,status)

!     determine the size of the image
 4    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!     check that it found both NAXIS1 and NAXIS2 keywords
 5    if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
       end if

!     check if the on-hand size of the image is not larger than
!     the real image size
	if(naxes(1)*naxes(2).lt.n1*n2) then
	  write(*,*) 'ERROR in the on-hand image size'
	  stop
	endif

!     initialize variables
      npixels=naxes(1)*naxes(2)
      group=1
      firstpix=1
      nullval=-999
      datamin=1.0E30
      datamax=-1.0E30

 6        call ftgpve(unit,group,firstpix,npixels,nullval,data_lin,anynull,status)

!     close the file and free the unit number
 7    call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 8    if (status .gt. 0)call printerror(status)
      end
      subroutine readtable

!     read and print data values from an ASCII or binary table

      integer status,unit,readwrite,blocksize,hdutype,ntable
      integer felem,nelems,nullj,diameter,nfound,irow,colnum
      real nulle,density
      character filename*40,nullstr*1,name*8,ttype(3)*10
      logical anynull

 1    status=0

!     Get an unused Logical Unit Number to use to open the FITS file
 2    call ftgiou(unit,status)

!     open the FITS file previously created by WRITEIMAGE
      filename='ATESTFILEZ.FITS'
      readwrite=0
 3    call ftopen(unit,filename,readwrite,blocksize,status)

!     loop twice, first reading the ASCII table, then the binary table
 4    do ntable=1,2

!         move to the next extension
 5        call ftmrhd(unit,1,hdutype,status)

          print *,' '
          if (hdutype .eq. 1)then
              print *,'Extension ',ntable,' is an ASCII table.'
          else if (hdutype .eq. 2)then
              print *,'Extension ',ntable,' is a binary table.'
          end if

!         read the TTYPEn keywords, which give the names of the columns
 6        call ftgkns(unit,'TTYPE',1,3,ttype,nfound,status)
          write(*,2000)ttype
2000      format(8x,3a10)

!         read the data, one row at a time, and print them out
          felem=1
          nelems=1
          nullstr=' '
          nullj=0
          nulle=0.
          do irow=1,6
              colnum=1
 7            call ftgcvs(unit,colnum,irow,felem,nelems,nullstr,name,anynull,status)
              colnum=2
 8            call ftgcvj(unit,colnum,irow,felem,nelems,nullj,diameter,anynull,status)
              colnum=3
 9            call ftgcve(unit,colnum,irow,felem,nelems,nulle,density,anynull,status)
              write(*,2001)irow,name,diameter,density
2001          format(i4,a10,i10,f10.2)
          end do
      end do

!     close the file and free the unit number
 10   call ftclos(unit, status)
      call ftfiou(unit, status)

!     check for any error, and if so print out error messages
 11   if (status .gt. 0)call printerror(status)
      end
      subroutine printerror(status)

!     Print out the FITSIO error messages to the user

      integer status
      character errtext*30,errmessage*80

!     check if status is OK (no error); if so, simply return
      if (status .le. 0)return

!     get the text string which describes the error
 1    call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

!     read and print out all the error messages on the FITSIO stack
 2    call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
      subroutine deletefile(filename,status)

!     A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!     simply return if status is greater than zero
      if (status .gt. 0)return

!     Get an unused Logical Unit Number to use to open the FITS file
 1    call ftgiou(unit,status)

!     try to open the file, to see if it exists
 2    call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
 3        call ftdelt(unit,status)
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
 4        call ftcmsg
      else
!         there was some other error opening the file; delete the file anyway
          status=0
 5        call ftcmsg
          call ftdelt(unit,status)
      end if

!     free the unit number for later reuse
 6    call ftfiou(unit, status)
      end
