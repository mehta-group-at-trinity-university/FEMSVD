      subroutine CalcEigenErrors(info,iparam,MatrixDim,H,LeadDim,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

      integer info,iparam(*),MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
      double precision H(LeadDim,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

      integer j
      double precision dnrm2
      double precision, allocatable :: HPsi(:),SPsi(:)

      if ( info .eq. 0) then

       if (iparam(5) .lt. NumStates) write(6,*) 'Not all states found'

! Compute the residual norm: ||  A*x - lambda*x ||

       allocate(HPsi(MatrixDim))
       allocate(SPsi(MatrixDim))
       do j = 1,NumStates
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,LeadDim,Psi(1,j),1,0.0d0,HPsi,1)
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
        call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
        Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
        Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
       enddo
       deallocate(HPsi)
       deallocate(SPsi)
      else
       write(6,*) ' '
       write(6,*) ' Error with _sband, info= ', info
       write(6,*) ' Check the documentation of _sband '
       write(6,*) ' '
      end if

      return
      end
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  Matrices A and B are stored in LAPACK-style band form.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost 
!  of n*nev if both are requested (in this case a separate array Z must be 
!  supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are called Ritz values and Ritz vectors respectively.  They are referred 
!  to as such in the comments that follow.  The computed orthonormal basis 
!  for the invariant subspace corresponding to these Ritz values is referred 
!  to as a Lanczos basis.
!
! \Usage
!   call dsband
!      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA, 
!        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V, 
!        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
!
! \Arguments
!
!  RVEC    Logical (INPUT)
!          Specifies whether Ritz vectors corresponding to the Ritz value 
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the associated Ritz vectors. 
!
!  HOWMNY  Character*1  (INPUT) 
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute all Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' , SELECT is not referenced.
!
!  D       Double precision array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by dsaupd transformed to
!          those of the original eigensystem A*z = lambda*B*z. If 
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
!          as the those of A*z = lambda*B*z.
!
!  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Lanczos basis array V computed by DSAUPD.
!
!  LDZ     Integer.  (INPUT) 
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   Double precision  (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
! 
!  N       Integer.  (INPUT) 
!          Dimension of the eigenproblem.  
! 
!  AB      Double precision array of dimension LDA by N. (INPUT)
!          The matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!  MB      Double precision array of dimension LDA by N. (INPUT)
!          The matrix M in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set. 
!          The j-th column of M is stored in the j-th column of the
!          array AB as follows:
!          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!          Not referenced if IPARAM(7) = 1
!
!  LDA     Integer. (INPUT)
!          Leading dimension of AB, MB, RFAC.
!
!  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
!          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
!          is invoked.  It is used to store the LU factors of
!          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
!          It is not referenced when IPARAM(7) = 1.
!
!  KL      Integer. (INPUT)
!          Max(number of subdiagonals of A, number of subdiagonals of M)
!
!  KU      Integer. (OUTPUT)
!          Max(number of superdiagonals of A, number of superdiagonals of M)
!
!  WHICH   Character*2.  (INPUT)
!          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
!          the following.
!  
!            'LM' -> want the NEV eigenvalues of largest magnitude.
!            'SM' -> want the NEV eigenvalues of smallest magnitude.
!            'LA' -> want the NEV eigenvalues of largest REAL part.
!            'SA' -> want the NEV eigenvalues of smallest REAL part.
!            'BE' -> Compute NEV eigenvalues, half from each end of the 
!                    spectrum.  When NEV is odd, compute one more from 
!                    the high end than from the low end. 
!
!          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
!          
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

!  NEV     Integer. (INPUT)
!          Number of eigenvalues of OP to be computed.
!   
!  TOL     Double precision scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = DLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH).
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          Represents the dimension of the Lanczos basis constructed
!          by dsaupd for OP.
!
!  V       Double precision array N by NCV.  (OUTPUT)
!          Upon INPUT: the NCV columns of V contain the Lanczos basis 
!                      vectors as constructed by dsaupd for OP.
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
!                       represent the Ritz vectors that span the desired 
!                       invariant subspace.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Lanczos basis vector array V computed by dsaupd. In this case
!          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
!          the desired Ritz vectors. 
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: 
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          It is set to 1 in this subroutine.  The user do not need
!          to set this parameter.
!          ------------------------------------------------------------
!          ISHIFT = 1: exact shifts with respect to the reduced 
!                      tridiagonal matrix T.  This is equivalent to 
!                      restarting the iteration with a starting vector 
!                      that is a linear combination of Ritz vectors 
!                      associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced. 
!
!          IPARAM(3) = MXITER
!          On INPUT:  max number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" eigenvalues.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used. 
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4,5; See under \Description of dsband for the 
!          five modes available.
!
!          IPARAM(8) = NP
!          Not referenced.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
! WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
!
! WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
!
! LWORKL   Integer.  (INPUT)
!          LWORKL must be at least NCV**2 + 8*NCV.
!
! IWORK    Integer array of dimension at least N. (WORKSPACE)
!          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
!          factorization of M or (A-SIGMA*M).
!            
! INFO     Integer.  (INPUT/OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)  
!                returns the number of wanted converged Ritz values.
!          =  3: No shifts could be applied during a cycle of the 
!                Implicitly restarted Arnoldi iteration. One possibility 
!                is to increase the size of NCV relative to NEV. 
!                See remark 4 in DSAUPD.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Informational error from LAPACK routine dsteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -13: HOWMNY must be one of 'A' or 'P'
!          = -14: DSAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current
!                   Arnoldi factorization.
!
!\Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dgbtrf  LAPACK band matrix factorization routine.
!     dgbtrs  LAPACK band linear system solve routine. 
!     dlacpy  LAPACK matrix copy routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the dot product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dgbmv   Level 2 BLAS that computes the band matrix vector product.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing 
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it.
!
!\Author    
!     Danny Sorensen
!     Richard Lehoucq
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
      subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,&
                       tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)
               implicit none

      character        which*2, bmat, howmny
      integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info  
      Double precision tol, sigma
      logical          rvec

      integer          iparam(*), iwork(*)
      logical          select(*)
      Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*)

      integer          ipntr(14)

      integer          ido, i, j, Row, Col, type, ierr

      Double precision one, zero
      parameter        (one = 1.0, zero = 0.0)

      Double precision ddot, dnrm2, dlapy2
      external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy

! iparam(3) : Max number of Arnoldi iterations
      iparam(3) = 10000
      iparam(7) = 3
      rvec = .TRUE.
      howmny = 'A'
      which = 'LM'
      bmat = 'G'
      type = 4  
      ido = 0
      iparam(1) = 1

      !print *, ncv,nev,n

      rfac = 0.0d0
      do i = 1,n
       do j = i,min(i+k,n)
        Row = k+1+i-j
        Col = j
        rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
       enddo
       do j = max(1,i-k),i-1
        Row = 2*k+1
        Col = j
        rfac(Row+i-j,j) = rfac(Row+j-i,i)
       enddo
      enddo

      call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
      if ( ierr .ne. 0 )  then
       print*, ' '
       print*, '_SBAND: Error with _gbtrf:',ierr
       print*, ' '
       go to 9000
      end if

  90  continue
     
      call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)

      if (ido .eq. -1) then
       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
       if (ierr .ne. 0) then
        print*, ' ' 
        print*, '_SBAND: Error with _gbtrs.'
        print*, ' ' 
        go to 9000
       end if
      else if (ido .eq. 1) then
       call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
       if (ierr .ne. 0) then 
        print*, ' '
        print*, '_SBAND: Error with _gbtrs.' 
        print*, ' '
        go to 9000
       end if
      else if (ido .eq. 2) then
       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
      else 
       if ( info .lt. 0) then
        print *, ' '
        print *, ' Error with _saupd info = ',info
        print *, ' Check the documentation of _saupd '
        print *, ' '
        go to 9000
       else 
        if ( info .eq. 1) then
         print *, ' '
         print *, ' Maximum number of iterations reached.'
         print *, ' '
        else if ( info .eq. 3) then
         print *, ' '
         print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
         print *, ' '
        end if
        if (iparam(5) .gt. 0) then
         call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
         if ( info .ne. 0) then
          print *, ' ' 
          print *, ' Error with _neupd = ', info
          print *, ' Check the documentation of _neupd '
          print *, ' ' 
          go to 9000
         endif
        endif
       endif
       go to 9000
      endif

      go to 90 

 9000 continue

    end subroutine MyDsband
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Mydggev(N,H,LDH,L,LDL,eval,evec)
  implicit none
  integer LDH,N,LDL,info
  double precision H(LDH,N),L(LDL,N),eval(N),evec(N,N)
  double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:)
  integer lwork,i,im,in,j

  allocate(alphar(N),alphai(N),beta(N))


  info = 0
  lwork = -1
  allocate(work(1))
  call dggev('N','V',N,H,LDH,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 
  do im = 1,N
     alphar(im)=0.0d0
     alphai(im)=0.0d0
     beta(im)=0.0d0
     eval(im)=0.0d0
     do in = 1,N
        evec(im,in)=0.0d0
     enddo
  enddo

  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dggev('N','V',N,H,LDH,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 

  do i = 1, N
     if (abs(alphai(i)).ge.1e-12) then
        print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i) ! 
     endif
     if(abs(beta(i)).ge.1e-12) then
        eval(i) = -alphar(i)/beta(i)
     endif
  enddo

  call deigsrt(eval,evec,N,N)

  do i = 1,N
     eval(i)=-eval(i)
  enddo

  deallocate(alphar,alphai,beta)
  deallocate(work)

end subroutine Mydggev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Mydgeev(N,H,LDH,eval,evec)
  implicit none
  integer LDH,N,info
  double precision H(LDH,N),eval(N),evec(N,N)
  double precision, allocatable :: alphar(:),alphai(:),work(:),VL(:,:)
  integer lwork,i,im,in,j

  allocate(alphar(N),alphai(N))

  info = 0
  lwork = -1
  allocate(work(1))
  write(6,*) "getting info"
  call dgeev('N','V',N,H,LDH,alphar,alphai,VL,1,evec,N,work,lwork,info) ! 
  write(6,*) "done getting info, info = ", info
  alphar(:)=0d0
  alphai(:)=0d0

  eval(:)=0d0
  evec(:,:)=0d0

  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))

  write(6,*) "diagonalizing now"
  call dgeev('N','V',N,H,LDH,alphar,alphai,VL,1,evec,N,work,lwork,info) ! 

  write(6,*) "checking for complex eigenvalues"
  do i = 1, N
     if (abs(alphai(i)).ge.1e-12) then
        print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i) ! 
     endif
     
     eval(i) = -alphar(i)

  enddo
  
  write(6,*) "sorting"
  call deigsrt(eval,evec,N,N)

  do i = 1,N
     eval(i)=-eval(i)
  enddo

  deallocate(alphar,alphai)
  deallocate(work)

end subroutine Mydgeev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CompSqrMatInv(A, N)
  implicit none
  integer N,info,lwk
  integer, allocatable :: ipiv(:)
  double precision, allocatable :: work(:)
  double complex A(N,N)
  allocate(ipiv(N))
  call zgetrf(N, N, A, N, ipiv, info)
  allocate(work(1))
  lwk = -1
  call zgetri(N, A, N, ipiv, work, lwk, info)
  lwk = work(1)
  deallocate(work)
  allocate(work(lwk))
  call zgetri(N, A, N, ipiv, work, lwk, info)
  deallocate(ipiv,work)
end subroutine CompSqrMatInv
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE deigsrt(d,v,n,np)
  INTEGER n,np
  double precision d(np),v(np,np)
  INTEGER i,j,k
  double precision p
  do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
        if(d(j).ge.p)then
           k=j
           p=d(j)
        endif
     enddo
     if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do j=1,n
           p=v(j,i)
           v(j,i)=v(j,k)
           v(j,k)=p
        enddo
     endif
  enddo
  return
END SUBROUTINE deigsrt
