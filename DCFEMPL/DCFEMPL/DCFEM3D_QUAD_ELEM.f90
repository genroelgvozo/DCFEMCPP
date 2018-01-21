function calc_appro(f,t)
    real(8)::calc_appro
    real(8),dimension(4)::f
    real(8),dimension(2)::t
    
    calc_appro = f(1) * (t(1)-1)*(t(2)-1)/4
    calc_appro = calc_appro + f(2) * (t(1)+1)*(t(2)-1)/(-4)
    calc_appro = calc_appro + f(3) * (t(1)-1)*(t(2)+1)/(-4)
    calc_appro = calc_appro + f(4) * (t(1)+1)*(t(2)+1)/4
end function

!> Calculate approximation of divariation of some function
!! @param f values of function in 4 points
!! @param t coordinates of point for calculating
!! @param n number of variable, which is divariative
function calc_approx_dif(f,t,n)
    real(8)::calc_approx_dif
    real(8),dimension(4)::f
    real(8),dimension(2)::t
    integer(4)::n
    
    calc_approx_dif = f(1) * (t(3-n)-1)/4
    calc_approx_dif = calc_approx_dif + f(2) * (t(3-n)+2*n-3)/(-4)
    calc_approx_dif = calc_approx_dif + f(3) * (t(3-n)+3-2*n)/(-4)
    calc_approx_dif = calc_approx_dif + f(4) * (t(3-n)+1)/4
end function

!> Calculate Beta of the Quad Element
!! @param t coordinates of point for calculate Beta (2)
!! @param x coordinates of element's vertexes (4,2)	
!! @param b Result matrix with Beta (2,2)
subroutine get_beta_quad(t, x, b)
    real(8),dimension(2)::t
    real(8),dimension(4,2)::x
    real(8),dimension(2,2)::b
    real(8)::det
    integer(4)::i,j
    
    do i=1,2
        do j=1,2
            b(i,j)=calc_approx_dif(x(:,i),t,j)
        end do
    end do
end subroutine


!> Calculate Jacobian of the Quad Element
!! @param t coordinates of point for calculate Jacobian (2)
!! @param x coordinates of element's vertexes (4,2)	
!! @param alpha Result matrix with Jacobian
subroutine get_jacoby_quad(t, x, alpha)
	use DCFEM3D
    real(8),dimension(2)::t
    real(8),dimension(4,2)::x
    real(8),dimension(2,2)::alpha
    
    real(8),dimension(2,2)::b
    real(8)::det
    
    call get_beta_quad(t,x,b)
	
    det=b(1,1)*b(2,2)-b(1,2)*b(2,1)

    do i=1,2
        do j=1,2
            b(i,j)=b(i,j)/det
        end do
    end do


    alpha(1,1)=b(2,2)
    alpha(2,2)=b(1,1)
    alpha(2,1)=-b(2,1)
    alpha(1,2)=-b(1,2)
end subroutine


!> Calculate deformation and stresses at quad element
!! @param e vector of solution at element (24)
!! @param lam value of the first parameter Lame
!! @param mu value of the second parameter Lame
!! @param t coordinates of point, in which we calculate stress (2)
!! @param x coordinate of element's vertexes (4,2)
!! @param eps array with result deformation (3,3)
!! @param stress array with result stress (3,3)
subroutine calc_stress_quad(e,lam,mu,t,x,eps,stress)
	
    real(8),dimension(24)::e
    real(8)::lam,mu
    real(8),dimension(2)::t
	real(8),dimension(4,2)::x
    real(8),dimension(3,3)::eps,stress
    real(8)::epse
        
    real(8),dimension(2,2)::alpha

    real(8),dimension(4,3)::u,v
    
    eps=0
    stress=0
    
    do i=1,4
        u(i,:)=e((i-1)*6+1:(i-1)*6+3)
        v(i,:)=e((i-1)*6+4:(i-1)*6+6)
    end do
    
    
    call get_jacoby_quad(t,x,alpha)

    
    do i=1,2
        eps(i,i) = calc_approx_dif(u(:,i),t,1) * alpha(1,i) + calc_approx_dif(u(:,i),t,2) * alpha(2,i)
    end do

    eps(3,3) = calc_appro(v(:,3),t)
    

    eps(2,1) = calc_approx_dif(u(:,1),t,1)*alpha(1,2)+calc_approx_dif(u(:,1),t,2)*alpha(2,2)
    eps(2,1) = eps(2,1)+calc_approx_dif(u(:,2),t,1)*alpha(1,1)+calc_approx_dif(u(:,2),t,2)*alpha(2,1)
    eps(2,1) = eps(2,1)/2
    eps(1,2) = eps(2,1)
    
    do i=1,2
        eps(3,i) = (calc_approx_dif(u(:,3),t,1)*alpha(1,i)+calc_approx_dif(u(:,3),t,2)*alpha(2,i)+calc_appro(v(:,i),t))/2
        eps(i,3) = eps(3,i)
    end do
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

!> Calculate functional Lagrange. Use quadtratures of Gauss to integrate
!! In this suite, use four points in vertex of quad with weights are equal 1/4
!! @param e vector with solution in the vertex
!! @param lam value of the first parameter Lame
!! @param mu value of the second parameter Lame
!! @param x coordinate of element's vertexes (4,2)
function calc_func_quad(e,lam,mu,x)
    
    real(8),dimension(24)::e
    real(8),dimension(4,2)::x
    real(8)::lam,mu
    real(8),dimension(3,3)::eps,stress
    real(8),dimension(2)::t
    real(8),dimension(2,2)::b
    real(8)::detb
    integer(4) :: Nint
    real(8),dimension(6) :: coord1
    real(8),dimension(6) :: coord2
    real(8) :: weight
    ! 1 - trap by vertexes
    ! 2 - int by Gauss
    integer(4) :: integ_scheme
    integ_scheme = 2
    
    
    if(integ_scheme == 1) then 
        
        Nint = 4
        
        coord1(1) = - 1.0D+00
        coord1(2) = - 1.0D+00
        coord1(3) = + 1.0D+00
        coord1(4) = + 1.0D+00
        coord2(1) = - 1.0D+00
        coord2(2) = + 1.0D+00
        coord2(3) = - 1.0D+00
        coord2(4) = + 1.0D+00
        
        weight =0.25D+00
    end if
    
    if(integ_scheme == 2) then
        
        Nint = 4
        
        coord1(1) = - 1.0D+00 / dsqrt(3.0D+00)
        coord1(2) = - 1.0D+00 / dsqrt(3.0D+00)
        coord1(3) = + 1.0D+00 / dsqrt(3.0D+00)
        coord1(4) = + 1.0D+00 / dsqrt(3.0D+00)
        coord2(1) = - 1.0D+00 / dsqrt(3.0D+00)
        coord2(2) = + 1.0D+00 / dsqrt(3.0D+00)
        coord2(3) = - 1.0D+00 / dsqrt(3.0D+00)
        coord2(4) = + 1.0D+00 / dsqrt(3.0D+00)
        
        weight = 1.0D+00
    end if

    
    calc_func_quad = 0
    
    do i=1,Nint
        t(1)=coord1(i)
        t(2)=coord2(i)

        call calc_stress_quad(e,lam,mu,t,x,eps,stress)
        call get_beta_quad(t,x,b)
        detb = b(1,1)*b(2,2)-b(1,2)*b(2,1)
        do m = 1,3
            do n= 1,3
                calc_func_quad = calc_func_quad + eps(m,n)*stress(m,n)*detb
            end do
        end do
            
    end do
    calc_func_quad = calc_func_quad*weight/2
end function


!> Calculate K matrix of the quad element
!! @param x coordinates of element's vertexes
!! @param lam value of the first parameter Lame
!! @param mu value of the second parameter Lame
!! @param K array with result
subroutine get_kelem_quad(x, lam, mu, K)
    real(8),dimension(4,2)::x
    real(8)::lam,mu
    real(8),dimension(24,24),intent(out)::K
    
    real(8),dimension(24)::e
    
    do i=1,24
        do j=1,24
            e=0
            e(i)=1
            e(j)=e(j)+1
            K(i,j)=calc_func_quad(e,lam,mu,x)
            
            e=0
            e(i)=1
            K(i,j)=K(i,j)-calc_func_quad(e,lam,mu,x)
            
            e=0
            e(j)=1
            K(i,j)=K(i,j)-calc_func_quad(e,lam,mu,x)
            
        end do
    end do
end subroutine
    
subroutine get_stress_elem_dcfem(N1,N2,h1,h2,U,eni,enj,numn,lam,mu,eps,stress)
    integer(4)::eni,enj,numn,N1,N2
    real(8)::lam,mu,h1,h2
    real(8),dimension(6*N1*N2)::U
    real(8),dimension(3,3)::eps,stress
    real(8),dimension(4,2)::x
    real(8),dimension(24)::etemp
    real(8),dimension(2)::tcoord
    integer(4),dimension(4)::iglob
    
    x(1,1)=0
    x(1,2)=0
    x(2,1)=h1
    x(2,2)=0
    x(3,1)=0
    x(3,2)=h2
    x(4,1)=h1
    x(4,2)=h2
    
    if(numn == 1) then
        tcoord(1)=-1
        tcoord(2)=-1
    end if
    
    if(numn == 2) then
        tcoord(1)=1
        tcoord(2)=-1
    end if
    
    if(numn == 3) then
        tcoord(1)=-1
        tcoord(2)=1
    end if
    
    if(numn == 4) then
        tcoord(1)=1
        tcoord(2)=1
    end if    
    
    iglob(1) = (enj-1)*N1 + eni
    iglob(2) = iglob(1) + 1
    iglob(3) = iglob(1) + N1
    iglob(4) = iglob(3) + 1 
    
    do i=1,4
        etemp(6*(i-1)+1:6*(i-1)+3) = U(3*iglob(i)-2:3*iglob(i))
        etemp(6*(i-1)+4:6*(i-1)+6) = U(3*N1*N2+3*iglob(i)-2:3*N1*N2+3*iglob(i))
    end do
    
    call calc_stress_quad(etemp,lam,mu,tcoord,x,eps,stress)

end subroutine
    
subroutine get_stress(N1,N2,numi,numj,y,h1,h2,lam,mu,eps1,stress1)
    integer(4)::numi,numj,N1,N2
    real(8)::h1,h2,lam,mu
    real(8),dimension(3,3)::eps1,stress1,etemp,stemp
    real(8),dimension(6*N1*N2)::y
    real(8),dimension(2) :: t
    real(8),dimension(4,2) :: x
    integer(4)::count
    eps1=0
    stress1=0
    count=0

    num = (numj-1)*N1 + numi
    if(numi>1) then
        if(numj>1) then
      
            call get_stress_elem_dcfem(N1,N2,h1,h2,y,numi-1,numj-1,4,lam,mu,etemp,stemp)
            
            eps1 = eps1 + etemp
            stress1 = stress1 + stemp
            count=count+1
        end if
        
        if(numj<N2) then
            
            call get_stress_elem_dcfem(N1,N2,h1,h2,y,numi-1,numj,2,lam,mu,etemp,stemp)
            
            eps1 = eps1 + etemp
            stress1 = stress1 + stemp
            count=count+1
        end if
    
    end if
    
    if(numi<N1) then
        if(numj>1) then
            
            call get_stress_elem_dcfem(N1,N2,h1,h2,y,numi,numj-1,3,lam,mu,etemp,stemp)
            
            eps1 = eps1 + etemp
            stress1 = stress1 + stemp
            count=count+1
        end if
        
        if(numj<N2) then
            
            call get_stress_elem_dcfem(N1,N2,h1,h2,y,numi,numj,1,lam,mu,etemp,stemp)
            
            eps1 = eps1 + etemp
            stress1 = stress1 + stemp
            count=count+1
        end if
    
    end if
    eps1 =eps1/count
    stress1 = stress1 / count
end subroutine