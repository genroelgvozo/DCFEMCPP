function calc_approx3(f,t)
    real(8)::calc_approx
    real(8),dimension(8)::f
    real(8),dimension(3)::t
    
    calc_approx3 = f(1) * (t(1)-1)*(t(2)-1)*(t(3)-1)/(-8)
    calc_approx3 = calc_approx3 + f(2) * (t(1)+1)*(t(2)-1)*(t(3)-1)/8
    calc_approx3 = calc_approx3 + f(3) * (t(1)-1)*(t(2)+1)*(t(3)-1)/8
    calc_approx3 = calc_approx3 + f(4) * (t(1)+1)*(t(2)+1)*(t(3)-1)/(-8)
    calc_approx3 = calc_approx3 + f(5) * (t(1)-1)*(t(2)-1)*(t(3)+1)/8
    calc_approx3 = calc_approx3 + f(6) * (t(1)+1)*(t(2)-1)*(t(3)+1)/(-8)
    calc_approx3 = calc_approx3 + f(7) * (t(1)-1)*(t(2)+1)*(t(3)+1)/(-8)
    calc_approx3 = calc_approx3 + f(8) * (t(1)+1)*(t(2)+1)*(t(3)+1)/8
end function

!> Calculate approximation of divariation of some function
!! @param f values of function in 4 points
!! @param t coordinates of point for calculating
!! @param n number of variable, which is divariative
function calc_approx_diff3(f,t,n)
    real(8)::calc_approx_diff
    real(8),dimension(8)::f
    real(8),dimension(3)::t
    integer(4)::n
    if(n==1) then
    calc_approx_diff3 = f(1) * (t(2)-1)*(t(3)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(2) *(t(2)-1)*(t(3)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(3) *(t(2)+1)*(t(3)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(4) *(t(2)+1)*(t(3)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(5) *(t(2)-1)*(t(3)+1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(6) *(t(2)-1)*(t(3)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(7) *(t(2)+1)*(t(3)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(8) *(t(2)+1)*(t(3)+1)/8
    end if
    if(n==2) then
    calc_approx_diff3 = f(1) * (t(1)-1)*(t(3)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(2) * (t(1)+1)*(t(3)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(3) * (t(1)-1)*(t(3)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(4) * (t(1)+1)*(t(3)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(5) * (t(1)-1)*(t(3)+1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(6) * (t(1)+1)*(t(3)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(7) * (t(1)-1)*(t(3)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(8) * (t(1)+1)*(t(3)+1)/8
    end if
    if(n==3) then
    calc_approx_diff3 = f(1) * (t(1)-1)*(t(2)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(2) * (t(1)+1)*(t(2)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(3) * (t(1)-1)*(t(2)+1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(4) * (t(1)+1)*(t(2)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(5) * (t(1)-1)*(t(2)-1)/8
    calc_approx_diff3 = calc_approx_diff3 + f(6) * (t(1)+1)*(t(2)-1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(7) * (t(1)-1)*(t(2)+1)/(-8)
    calc_approx_diff3 = calc_approx_diff3 + f(8) * (t(1)+1)*(t(2)+1)/8
    end if
end function
    
subroutine get_beta_brick(t, x, b)
    real(8),dimension(3)::t
    real(8),dimension(8,3)::x
    real(8),dimension(3,3)::b
    
    do i=1,3
        do j=1,3
            b(i,j)=calc_approx_diff3(x(:,i),t,j)
        end do
    end do
end subroutine


!> Calculate Jacobian of the Quad Element
!! @param t coordinates of point for calculate Jacobian (2)
!! @param x coordinates of element's vertexes (4,2)	
!! @param alpha Result matrix with Jacobian
subroutine get_jacoby_brick(t, x, alpha)
    use lapack95
    use f95_precision
    real(8),dimension(3)::t
    real(8),dimension(8,3)::x
    real(8),dimension(3,3)::alpha
    
    real(8),dimension(3,3)::d
    integer(4),dimension(3)::ipiv
    real(8)::det

    call get_beta_brick(t,x,d)

    
    call getrf(d,ipiv)
    call getri(d,ipiv)
    
    alpha = d
end subroutine    
    
    
subroutine calc_stress_brick(e,lam,mu,t,x,eps,stress)

    real(8),dimension(24)::e
    real(8)::lam,mu
    real(8),dimension(3)::t
	real(8),dimension(8,3)::x
    real(8),dimension(3,3)::eps,stress
    
        
    real(8),dimension(3,3)::alpha

    real(8),dimension(8,3)::u
    
    do i=1,8
        u(i,:)=e(3*i-2:3*i)
    end do
    
    
    call get_jacoby_brick(t,x,alpha)

    
    do i=1,3
        eps(i,i)=0
        do j=1,3
            eps(i,i)=eps(i,i) + calc_approx_diff3(u(:,i),t,j) * alpha(j,i)
        end do
        
    end do

        

    eps(2,1) = calc_approx_diff3(u(:,1),t,1)*alpha(1,2)+calc_approx_diff3(u(:,1),t,2)*alpha(2,2)+calc_approx_diff3(u(:,1),t,3)*alpha(3,2)
    eps(2,1) = eps(2,1)+calc_approx_diff3(u(:,2),t,1)*alpha(1,1)+calc_approx_diff3(u(:,2),t,2)*alpha(2,1)+calc_approx_diff3(u(:,2),t,3)*alpha(3,1)
    eps(2,1) = eps(2,1)/2
    eps(1,2) = eps(2,1)
    
    eps(3,1) = calc_approx_diff3(u(:,1),t,1)*alpha(1,3)+calc_approx_diff3(u(:,1),t,2)*alpha(2,3)+calc_approx_diff3(u(:,1),t,3)*alpha(3,3)
    eps(3,1) = eps(3,1)+calc_approx_diff3(u(:,3),t,1)*alpha(1,1)+calc_approx_diff3(u(:,3),t,2)*alpha(2,1)+calc_approx_diff3(u(:,3),t,3)*alpha(3,1)
    eps(3,1) = eps(3,1)/2
    eps(1,3) = eps(3,1)
    
    eps(3,2) = calc_approx_diff3(u(:,2),t,1)*alpha(1,3)+calc_approx_diff3(u(:,2),t,2)*alpha(2,3)+calc_approx_diff3(u(:,2),t,3)*alpha(3,3)
    eps(3,2) = eps(3,2)+calc_approx_diff3(u(:,3),t,1)*alpha(1,2)+calc_approx_diff3(u(:,3),t,2)*alpha(2,2)+calc_approx_diff3(u(:,3),t,3)*alpha(3,2)
    eps(3,2) = eps(3,2)/2
    eps(2,3) = eps(3,2)
    
    epse = eps(1,1) + eps(2,2) + eps(3,3)
    do i=1,3
        do j=1,3
            stress(i,j) = 2 * mu * eps(i,j)
            if(i==j) then
                stress(i,j) = stress(i,j) + lam * epse
            end if
        end do
    end do
    
    end subroutine
    
function calc_func_brick(e,lam,mu,x)
    
    real(8),dimension(24)::e
    real(8),dimension(8,3)::x
    real(8)::lam,mu
    real(8),dimension(3,3)::eps,stress
    real(8),dimension(3)::t
    real(8),dimension(3,3)::d
    real(8)::detb
    integer(4) :: Nint
    real(8),dimension(3,14) :: coord
    real(8),dimension(14) :: weight
    real(8) :: sq
    ! 1 - trap by vertexes
    ! 2 - int by Gauss
    integer(4) :: integ_scheme = 2
        if(integ_scheme == 0) then
        Nint = 1
        coord(1,1) = 0

        coord(2,1) = 0

        
        coord(3,1) = 0
        weight=8
    end if
    if(integ_scheme == 1) then
        Nint = 8
        coord(1,1) = - 1.0D+00
        coord(1,2) = - 1.0D+00
        coord(1,3) = - 1.0D+00
        coord(1,4) = - 1.0D+00
        coord(1,5) = + 1.0D+00
        coord(1,6) = + 1.0D+00
        coord(1,7) = + 1.0D+00
        coord(1,8) = + 1.0D+00
        
        coord(2,1) = - 1.0D+00
        coord(2,2) = - 1.0D+00
        coord(2,3) = + 1.0D+00
        coord(2,4) = + 1.0D+00
        coord(2,5) = - 1.0D+00
        coord(2,6) = - 1.0D+00
        coord(2,7) = + 1.0D+00
        coord(2,8) = + 1.0D+00
        
        coord(3,1) = - 1.0D+00
        coord(3,2) = + 1.0D+00
        coord(3,3) = - 1.0D+00
        coord(3,4) = + 1.0D+00
        coord(3,5) = - 1.0D+00
        coord(3,6) = + 1.0D+00
        coord(3,7) = - 1.0D+00
        coord(3,8) = + 1.0D+00
        weight=1.0/8
    end if
    
    if(integ_scheme == 2) then
        
        Nint = 8
        
        coord(1,1) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(1,2) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(1,3) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(1,4) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(1,5) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(1,6) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(1,7) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(1,8) = + 1.0D+00 / dsqrt(3.0D+00)
        
        coord(2,1) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(2,2) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(2,3) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(2,4) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(2,5) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(2,6) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(2,7) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(2,8) = + 1.0D+00 / dsqrt(3.0D+00)
        
        coord(3,1) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(3,2) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(3,3) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(3,4) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(3,5) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(3,6) = + 1.0D+00 / dsqrt(3.0D+00)
        coord(3,7) = - 1.0D+00 / dsqrt(3.0D+00)
        coord(3,8) = + 1.0D+00 / dsqrt(3.0D+00)

        weight = 1.0D+00
    end if
    
    if(integ_scheme == 3) then
        sq = dsqrt(19.0D+00/33.0D+00)
        
        Nint = 14
        coord= 0
        do i=1,3
            coord(i,2*i-1) = -sq
            coord(i,2*i) = sq
        end do
        weight(1:6) = 320.0/361.0
        
        coord(7:10,1)=-sq
        coord(11:14,1)=sq
        
        do i=1,2
            coord(6+4*(i-1)+1:6+4*(i-1)+2,2)=-sq
            coord(6+4*(i-1)+3:6+4*(i-1)+4,2)=sq
        end do
        do i=1,4
            coord(6+2*i-1,3)=-sq
            coord(6+2*i,3)=sq
        end do
        
        weight(7:14)=121.0/361.0
    end if
    
    calc_func_brick = 0
    
    do i=1,Nint
        t(1)=coord(1,i)
        t(2)=coord(2,i)
        t(3)=coord(3,i)
        call calc_stress_brick(e,lam,mu,t,x,eps,stress)
        call get_beta_brick(t,x,d)
        detb=d(1,1)*d(2,2)*d(3,3)+d(1,2)*d(2,3)*d(3,1)+d(1,3)*d(2,1)*d(3,2)-d(1,3)*d(2,2)*d(3,1)-d(1,2)*d(2,1)*d(3,3)-d(1,1)*d(2,3)*d(3,2)
        
        do m = 1,3
            do n= 1,3
                calc_func_brick = calc_func_brick + eps(m,n)*stress(m,n)*detb*weight(i)
            end do
        end do
            
    end do
    calc_func_brick = calc_func_brick/2
end function

    
subroutine get_kelem_shell(h1,h2, E,nu, t, K)
    real(8)::h1,h2
    real(8)::E,nu
    real(8),dimension(12,12),intent(out)::K
    real(8),dimension(12,12)::K1,K2,K3,K4,L
    real(8)::a,b,ro,Dx,Dy,Dxy,D1
    a = h1/2
    b = h2/2
    Dx = E*t**3/12/(1-nu*nu)
    Dy = E*t**3/12/(1-nu*nu)
    D1 = E*t**3/12/(1-nu*nu)*nu
    Dxy = E*t**3/12/(1-nu*nu)*(1-nu)/2
    K1(1,:) = (/     60,     0,  30,  30,   0,  15, -60,  0,  30,  -30,  0, 15 /)
    K1(2,:) = (/    0,     0,     0,   0,   0,   0,   0,  0,   0,    0,  0,  0 /)
    K1(3,:) = (/   30,   0,    20,    15,   0,  10, -30,  0,  10,  -15,  0,  5 /)
    K1(4,:) = (/   30,   0,  15,    60,     0,  30, -30,  0,  15,  -60,  0, 30 /)
    K1(5,:) = (/    0,   0,   0,   0,     0,     0,   0,  0,   0,    0,  0,  0 /)
    K1(6,:) = (/   15,   0,  10,  30,   0,    20,   -15,  0,   5,  -30,  0, 10 /)
    K1(7,:) = (/  -60,   0, -30, -30,   0, -15,    60,    0, -30,   30,  0, -15 /)
    K1(8,:) = (/    0,   0,   0,   0,   0,   0,   0,    0,     0,    0,  0,  0 /)
    K1(9,:) = (/   30,   0,  10,  15,   0,   5, -30,  0,    20,    -15,  0, 10 /)
    K1(10,:) = (/ -30,   0, -15, -60,   0, -30,  30,  0, -15,     60,    0, -30 /)
    K1(11,:) = (/   0,   0,   0,   0,   0,   0,   0,  0,   0,    0,    0,    0 /)
    K1(12,:) = (/  15,   0,   5,  30,   0,  10, -15,  0,  10,  -30,  0,   20 /)
    K1 = K1 * b**2/a**2
    
    K2(1,:) = (/     60,   -30,   0, -60, -30,   0,  30, -15,   0,  -30, -15,   0 /)
    K2(2,:) = (/  -30,    20,     0,  30,  10,   0,  15,   5,   0,   15,   5,   0 /)
    K2(3,:) = (/    0,   0,     0,     0,   0,   0,   0,   0,   0,    0,   0,   0 /)
    K2(4,:) = (/  -60,  30,   0,    60,    30,   0, -30,  15,   0,   30,  15,   0 /)
    K2(5,:) = (/  -30,  10,   0,  30,    20,     0, -15,   5,   0,   15,  10,   0 /)
    K2(6,:) = (/    0,   0,   0,   0,   0,     0,     0,   0,   0,    0,   0,   0 /)
    K2(7,:) = (/   30, -15,   0, -30, -15,   0,    60,   -30,   0,  -60, -30,   0 /)
    K2(8,:) = (/  -15,  10,   0,  15,   5,   0, -30,    20,     0,   30,  10,   0 /)
    K2(9,:) = (/    0,   0,   0,   0,   0,   0,   0,   0,     0,      0,   0,   0 /)
    K2(10,:) = (/ -30,  15,   0,  30,  15,   0, -60,  30,   0,    60,     30,   0 /)
    K2(11,:) = (/ -15,   5,   0,  15,  10,   0, -30,  10,   0,  30,    20,      0 /)
    K2(12,:) = (/   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,      0 /)
    K2 = K2*a**2/b**2
    
    L = 0
    L(1,1) = 1
    L(2,2) = 2*b
    L(3,3) = 2*a
    L(4,4) = 1
    L(5,5) = 2*b
    L(6,6) = 2*a
    L(7,7) = 1
    L(8,8) = 2*b
    L(9,9) = 2*a
    L(10,10) = 1
    L(11,11) = 2*b
    L(12,12) = 2*a
    
    K = Dx*K1+Dy*K2+D1*K3+Dxy*K4
    
end subroutine
    
    
subroutine BUILDFEM(N1,N2,l1,l2,E,nu,t,K,R)
    integer(4)::N1,N2
    real(8)::l1,l2,lam,mu,t
    real(8),dimension(12*N1*N2,12*N1*N2)::K
    real(8),dimension(12*N1*N2)::R
    
    real(8),dimension(12,12)::ki
    real(8)::h1,h2

    K=0
    R=0

    h1 = l1 / (N1-1)
    h2 = l2 / (N2-1)


    print*,"lam FEM = ",lam
    print*,"mu FEM = ",mu
    call get_kelem_shell(h1,h2,E,nu,t,ki)
    open(2,file="Ki.txt")
    write(2,'(24F18.8)')ki
    close(2)
    do i=1,N1-1
        do j=1,N2-1
            do k1=1,N3-1
            
            iloc(1) = (k1-1)*N1*N2+(j-1)*N1 + i
            iloc(2) = iloc(1) + 1
            iloc(3) = iloc(1)+ N1
            iloc(4) = iloc(2)+ N1
            iloc(5) = iloc(1) + N1*N2
            iloc(6) = iloc(2) + N1*N2
            iloc(7) = iloc(3) + N1*N2
            iloc(8) = iloc(4) + N1*N2

            do m=1,8
                do n=1,8
                    K((iloc(m)-1)*3+1:(iloc(m)-1)*3+3,(iloc(n)-1)*3+1:(iloc(n)-1)*3+3) =K((iloc(m)-1)*3+1:(iloc(m)-1)*3+3,(iloc(n)-1)*3+1:(iloc(n)-1)*3+3) + ki(3*m-2:3*m,3*n-2:3*n)
                end do
            end do
            !k((i11-1)*2+1:(i11-1)*2+4,(i11-1)*2+1:(i11-1)*2+4) = k((i11-1)*2+1:(i11-1)*2+4,(i11-1)*2+1:(i11-1)*2+4) + ki(1:4,1:4)
            !k((i11-1)*2+1:(i11-1)*2+4,(i21-1)*2+1:(i21-1)*2+4) = k((i11-1)*2+1:(i11-1)*2+4,(i21-1)*2+1:(i21-1)*2+4) + ki(1:4,5:8)
            !k((i21-1)*2+1:(i21-1)*2+4,(i11-1)*2+1:(i11-1)*2+4) = k((i21-1)*2+1:(i21-1)*2+4,(i11-1)*2+1:(i11-1)*2+4) + ki(5:8,1:4)
            !k((i21-1)*2+1:(i21-1)*2+4,(i21-1)*2+1:(i21-1)*2+4) = k((i21-1)*2+1:(i21-1)*2+4,(i21-1)*2+1:(i21-1)*2+4) + ki(5:8,5:8)
            end do
        end do
    end do
    
    end subroutine
    
    
subroutine get_stress_elem(N1,N2,N3,h1,h2,h3,y,eni,enj,enk,numn,lam,mu,eps1,stress1)
    integer(4)::eni,enj,enk,numn,N1,N2,N3
    real(8)::lam,mu,h1,h2,h3
    real(8),dimension(3*N1*N2*N3)::y
    real(8),dimension(3,3)::eps1,stress1
    real(8),dimension(8,3)::x
    real(8),dimension(24)::etemp
    real(8),dimension(3)::tcoord
    integer(4),dimension(8)::iglob
    x(1,1)=0
    x(1,2)=0
    x(1,3)=0
    
    x(2,1)=h1
    x(2,2)=0
    x(2,3)=0
    
    x(3,1)=0
    x(3,2)=h2
    x(3,3)=0
    
    x(4,1)=h1
    x(4,2)=h2
    x(4,3)=0
    
    x(5,1)=0
    x(5,2)=0
    x(5,3)=h3
    
    x(6,1)=h1
    x(6,2)=0
    x(6,3)=h3
    
    x(7,1)=0
    x(7,2)=h2
    x(7,3)=h3
    
    x(8,1)=h1
    x(8,2)=h2
    x(8,3)=h3
    
    
    if(numn == 1) then
        tcoord(1)=-1
        tcoord(2)=-1
        tcoord(3)=-1
    end if
    
    if(numn == 2) then
        tcoord(1)=1
        tcoord(2)=-1
        tcoord(3)=-1
    end if
    
    if(numn == 3) then
        tcoord(1)=-1
        tcoord(2)=1
        tcoord(3)=-1
    end if
    
    if(numn == 4) then
        tcoord(1)=1
        tcoord(2)=1
        tcoord(3)=-1
    end if
    
    if(numn == 5) then
        tcoord(1)=-1
        tcoord(2)=-1
        tcoord(3)=1
    end if
    
    if(numn == 6) then
        tcoord(1)=1
        tcoord(2)=-1
        tcoord(3)=1
    end if
    
    if(numn == 7) then
        tcoord(1)=-1
        tcoord(2)=1
        tcoord(3)=1
    end if
    
    if(numn == 8) then
        tcoord(1)=1
        tcoord(2)=1
        tcoord(3)=1
    end if
    iglob(1) = (enk-1)*N1*N2+(enj-1)*N1 + eni
    iglob(2) = iglob(1) +1
    iglob(3) = iglob(1)+N1
    iglob(4) = iglob(3)+1 
    iglob(5) = iglob(1)+N1*N2
    iglob(6) = iglob(2)+N1*N2
    iglob(7) = iglob(3)+N1*N2
    iglob(8) = iglob(4)+N1*N2
    do i=1,8
        etemp(3*i-2:3*i) = y(3*iglob(i)-2:3*iglob(i))
    end do
    
    call calc_stress_brick(etemp,lam,mu,tcoord,x,eps1,stress1)
end subroutine
    
subroutine get_stress_node(N1,N2,N3,h1,h2,h3,y,ni,nj,nk,lam,mu,eps,stress)
    integer(4)::N1,N2,N3,ni,nj,nk
    
    real(8),dimension(3*N1*N2*N3)::y
    real(8)::lam,mu,h1,h2,h3
    real(8),dimension(3,3)::eps,stress,etemp,stemp
    integer(4)::count
    eps = 0
    stress = 0
    etemp=0
    stemp=0
    count=0

    if(ni>1) then
        if(nj>1) then
            if(nk>1) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni-1,nj-1,nk-1,8,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
            if(nk<N3) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni-1,nj-1,nk,4,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
        end if
        if(nj<N2) then
            if(nk>1) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni-1,nj,nk-1,6,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
            if(nk<N3) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni-1,nj,nk,2,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
        end if
    end if
        
    if(ni<N1) then
        if(nj>1) then
            if(nk>1) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni,nj-1,nk-1,7,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
            if(nk<N3) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni,nj-1,nk,3,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
        end if
        if(nj<N2) then
            if(nk>1) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni,nj,nk-1,5,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
            if(nk<N3) then
                call get_stress_elem(N1,N2,N3,h1,h2,h3,y,ni,nj,nk,1,lam,mu,etemp,stemp)
                count=count+1
                eps = eps + etemp
                stress = stress + stemp
            end if
        end if
    end if
        
    
    eps = eps / count
    stress = stress / count

end subroutine
