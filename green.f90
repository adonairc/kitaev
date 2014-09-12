program green
	implicit none

	external ZGETRF
        external ZGETRI

	! Variables definitions
	real*16, parameter :: t = 1.0
	real*16, parameter :: gammal = 0.004*t
	real*16, parameter :: t0 = 20*gammal
	real*16,parameter :: mu0 = 5*gammal !mu0 = -edot
	real*16,parameter :: mu = 0
	real*16, parameter :: delta = 0.2*t
	
	real*16, parameter :: eta = 10.0D-6
	real*16, parameter :: min = -20*gammal
	real*16, parameter :: max = 20*gammal
	real*16, parameter :: PI = 3.141592D+00
	real*16 :: omega,step

	complex*16 :: imag
	complex*16,dimension(2,2) :: ID, V, W, V0, W0, W0T, WT, g, gd, gtr, gtl, gdd, gnn, gb, inv
	complex*16,dimension(2) :: work
	
	integer,parameter :: N = 30000 !Number of sites
	!integer,parameter :: B = 2 ! Bulk site
	integer,parameter :: qtd = 2000 ! Discretization points
	integer,dimension(2) :: ipiv
	integer :: i,j,info

	imag = (0.0D+00,1.0D+00)

	ID(1,1) = (1,0)
	ID(1,2) = (0,0)
	ID(2,1) = (0,0)
	ID(2,2) = (1,0)	

	W(1,1) = (0,0)
	W(1,2) = cmplx(0.0,0.5*(delta+t))
	W(2,1) = cmplx(0.0,0.5*(delta-t))
	W(2,2) = (0,0)

	WT = transpose(conjg(W))

	W0(1,1) = (0,0)
	W0(1,2) = cmplx(0.0,0.5*t0)
	W0(2,1) = cmplx(0.0,-0.5*t0)
	W0(2,2) = (0,0)

	W0T = transpose(conjg(W0))

	V(1,1) = (0,0)
	V(1,2) = cmplx(0.0,0.5*mu)
	V(2,1) = cmplx(0.0,-0.5*mu)
	V(2,2) = (0,0)
	
	V0(1,1) = (0,0)
	V0(1,2) = cmplx(0.0,0.5*mu0)
	V0(2,1) = cmplx(0.0,-0.5*mu0)
	V0(2,2) = (0,0)

	step = (max-min)/qtd

	open(unit=1,file="dados")

	omega = min
	do while (omega < (max+step))

		g(1,1) = 2.0/cmplx(omega,eta)
		g(1,2) = (0,0)
		g(2,1) = (0,0)
		g(2,2) = 2.0/cmplx(omega,eta)
		
		gd(1,1) = 2.0/cmplx(omega,eta+gammal)
		gd(1,2) = (0,0)
		gd(2,1) = (0,0)
		gd(2,2) = 2.0/cmplx(omega,eta+gammal)

		inv = ID-matmul(g,V)

		call ZGETRF(2,2,inv,2,ipiv,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if

		call ZGETRI(2,inv,2,ipiv,work,2,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if
			
		gb = matmul(inv,g)
		gtr = gb
		gtl = gb


		! Right branch
		do j=1,N
			inv = ID-matmul(matmul(gb,WT),matmul(gtr,W))

			call ZGETRF(2,2,inv,2,ipiv,info)
			if(info .ne. 0) then
				write(*,*)"failed",info
			end if

			call ZGETRI(2,inv,2,ipiv,work,2,info)
			if(info .ne. 0) then
				write(*,*)"failed",info
			end if

			gtr = matmul(inv,gb)
		enddo

		! Green's function of the dot
		inv = ID-matmul(gd,V0)
		call ZGETRF(2,2,inv,2,ipiv,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if

		call ZGETRI(2,inv,2,ipiv,work,2,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if
		
		gdd = matmul(inv,gd)
		

		! Combination with first site
! 		inv = ID-matmul(matmul(gb,W0),matmul(gtl,W0T))

! 		call ZGETRF(2,2,inv,2,ipiv,info)
! 		if(info .ne. 0) then
! 			write(*,*)"failed",info
! 		end if

! 		call ZGETRI(2,inv,2,ipiv,work,2,info)
! 		if(info .ne. 0) then
! 			write(*,*)"failed",info
! 		end if

! 		gtl = matmul(inv,gb)

! 		! Left branch
! 		do j=2,B
! 			inv = ID-matmul(matmul(gb,W),matmul(gtl,WT))

! 			call ZGETRF(2,2,inv,2,ipiv,info)
! 			if(info .ne. 0) then
! 				write(*,*)"failed",info
! 			end if

! 			call ZGETRI(2,inv,2,ipiv,work,2,info)
! 			if(info .ne. 0) then
! 				write(*,*)"failed",info
! 			end if

! 			gtl = matmul(inv,gb)
! 		enddo

		! Combination at bulk
		inv = ID-matmul(matmul(gtr,W0),matmul(gdd,W0T))
		call ZGETRF(2,2,inv,2,ipiv,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if

		call ZGETRI(2,inv,2,ipiv,work,2,info)
		if(info .ne. 0) then
			write(*,*)"failed",info
		end if

		gnn = matmul(inv,gtr)

		write(1,*) omega/gammal,PI*gammal*(-1.0D+00/PI)*aimag(0.25*(gnn(1,1)+gnn(2,2)+imag*gnn(1,2)-imag*gnn(2,1)))
		!write(1,*) omega,(-1.0D+00/PI)*(eta*PI)*aimag(gt(2,2))

		omega = omega + step
	enddo	
	close(1)
end program green
