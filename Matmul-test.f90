program main
   real(8),allocatable  ::  M1(:,:),M2(:,:)
   real(8)              ::  Sum_T1,sum_T2
   integer              ::  n
   real              ::  timer1,timer2,timer3
   
   do n = 1000,5000,200
      allocate(M1(n,n))
      allocate(M2(n,n))
      call  RANDOM_NUMBER(M1)
      call  RANDOM_NUMBER(M2)

      call  cpu_time(timer1)
      call  myWritesum(M1,M2,n,sum_T1)
      call  cpu_time(timer2)
      call  interBLAS(M1,M2,n,sum_T2)
      call  cpu_time(timer3)

      print '("data n*n:   ",I5,"    muSum:",f10.5,"   BLAS3:",f10.5  )',n,timer2-time1,timer3-time2

      deallocate(M1)
      deallocate(M2)
   enddo

end program

!=======
! Using BLAS3
!=======
subroutine interBLAS(M1,M2,n,sum_T)
integer :: n,i,j
real(8) :: M1(n,n)
real(8) :: M2(n,n)
real(8) :: sum_T
real(8) :: alpha,beta
real(8),allocatable :: R(:,:)

alpha = 1
beta = 0

allocate(R(n,n))
call DGEMM('N','T',n,n,n,alpha,M1,n,M2,n,beta,R,n)
sum_T = 0
do i = 1,n
       sum_T = sum_T + R(i,i)
enddo

deallocate(R)



end subroutine


!=====
!  My subroutine
!=====
subroutine myWritesum(M1,M2,n,sum_T)

integer :: i,j,n
real(8) :: M1(n,n)
real(8) :: M2(n,n)
real(8) :: sum_T

sum_T = 0
do i = 1,n
   do j = 1,n
       sum_T = sum_T + M1(i,j)*M2(i,j)
   enddo
enddo

end subroutine
