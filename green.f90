subroutine print_matrix(M)
	! Print 2x2 matrices
	complex*16,dimension(2,2),intent(in) :: M
	integer :: i,j
	write(*,*)
	do, i=1,2
		!write(*,'(50g15.10)') ('(',real(M(i,j)),',',aimag(M(i,j)),')',j=1,2)
		write(*,*) ('(',real(M(i,j)),',',aimag(M(i,j)),')',j=1,2)
	enddo

end subroutine print_matrix

program green
	implicit none

	external ZGETRF
    external ZGETRI

	! Definição das variáveis
	real*16,parameter :: mu = 0
	real*16, parameter :: eta = 10.0D-5
	real*16, parameter :: t = 0.5
	real*16, parameter :: delta = 0.2*t
	real*16, parameter :: min = -2.0D+00
	real*16, parameter :: max = 2.0D+00
	real*16, parameter :: PI = 3.141592D+00
	real*16 :: omega,step

	complex*16 :: imag
	complex*16,dimension(2,2) :: ID, V, W, WT, g, gt, gb, inv
	complex*16,dimension(2) :: work
	
	integer,parameter :: N = 30000 ! Number of sites
	integer,parameter :: qtd = 1000 ! Discretization
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

	V(1,1) = (0,0)
	V(1,2) = cmplx(0.0,mu)
	V(2,1) = cmplx(0.0,-mu)
	V(2,2) = (0,0)
	
	step = (max-min)/qtd

	open(unit=1,file="dados")

	omega = min
	do while (omega < (max+step))

		g(1,1) = 2.0/cmplx(omega,eta)
		g(1,2) = (0,0)
		g(2,1) = (0,0)
		g(2,2) = 2.0/cmplx(omega,eta)
		
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
		gt = gb

		do j=1,N
			inv = ID-matmul(matmul(gb,W),matmul(gt,WT))
			
			!call print_matrix(inv)

			call ZGETRF(2,2,inv,2,ipiv,info)
			if(info .ne. 0) then
				write(*,*)"failed",info
			end if

			call ZGETRI(2,inv,2,ipiv,work,2,info)
			if(info .ne. 0) then
				write(*,*)"failed",info
			end if

			gt = matmul(inv,gb)
		enddo

		write(1,*) omega,(-1.0D+00/PI)*aimag(0.25*(gt(1,1)+gt(2,2)+imag*gt(1,2)-imag*gt(2,1)))
		!write(1,*) omega,(-1.0D+00/PI)*(eta*PI)*aimag(gt(2,2))

		omega = omega + step
	enddo	
	close(1)
end program green