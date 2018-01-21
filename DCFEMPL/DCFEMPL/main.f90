

    
program DCFEM2D
    use lapack95
    use f95_precision
    integer(4) :: N1,N2,N3 !Количество узлов вдоль x1 и x2
    real(8)::xb1,xb2,xb3,x
    inTeger(4)::nf
    real(8)::l1,l2,l31,l32,l3
    real(8)::h1,h2,h3
    real(8)::E,nu
    integer(4)::nv,ne
    real(8),dimension(:),allocatable:: lam,mu,x1,x2,x3
    real(8),dimension(:,:),allocatable::coords
    integer(4),dimension(:,:),allocatable::elems
    integer(4),dimension(:),allocatable::nkf
    real(8),dimension(:),allocatable::xf,R,y
    real(8),dimension(:,:),allocatable::K,f1,f2,f3
    real(8),dimension(8,8)::ki
    real(8),dimension(3,3)::eps,stress
    real(8),dimension(2)::t
    integer(4)::err,inum,inn
    real(8)::res
    N1 = 11
    N2 = 11
    N3 = 11
        
    nv = N1*N2
    ne = (N1-1)*(N2-1)
    allocate(x1(N1),x2(N2),x3(N3),lam(ne),mu(ne))

    
    E= 3000
    nu = 0.16
    l1 = 300
    l2 = 300
    l31 = 300
    l32 = 300
    l3 = l31+l32
    
    
    h1 = l1 / (N1-1)
    h2 = l2 / (N2-1)
    h3 = l32 / (N3-1)
    do i=1,N1
        x1(i)=h1*(i-1)
    end do
    do i=1,N2
        x2(i)=h2*(i-1)
    end do
    do i=1,N3
        x3(i)=l31+h3*(i-1)
    end do
    
    allocate(coords(2,nv),elems(4,ne))
    
    do i=1,N1
        do j=1,N2
            coords(1,(j-1)*N1+i) = x1(i)
            coords(2,(j-1)*N1+i) = x2(j)
        end do
    end do
    
    do i=1,N1-1
        do j=1,N2-1
            elems(1,(j-1)*(N1-1)+i)=i+(j-1)*N1
            elems(2,(j-1)*(N1-1)+i)=i+(j-1)*N1+1
            elems(3,(j-1)*(N1-1)+i)=i+j*N1
            elems(4,(j-1)*(N1-1)+i)=i+j*N1+1
        end do
    end do
        
    print*,"lam=",E*nu/(1+nu)/(1-2*nu)
    print*,"mu=",E/(2+2*nu)
    do i=1,ne
        lam(i)=E*nu/(1+nu)/(1-2*nu)
        mu(i)=E/(2+2*nu)
    end do
    
    nf = 1
    allocate(xf(nf),f1(nf,nv),f2(nf,nv),f3(nf,nv))
    f1= 0
    f2= 0
    f3= 0
    
    xf(1)=150
    
    f2(1,nv-N1+1)=-1500
    do i=2,N1-1
        f2(1,nv-N1+i)=-3000
    end do
    f2(1,nv)=-1500

    call INIT()
    
    call BUILDSDE2(1, nv , ne, coords, elems, lam, mu, nf, xf, f1, f2, f3, err)
    
    allocate(K(3*N1*N2*N3,3*N1*N2*N3),R(3*N1*N2*N3))
    R=0
    call BUILDFEM(N1,N2,N3,l1,l2,l32,lam(1),mu(1),K,R)
    call BUILDLEFTBOUND(N1*N2,0.0D+00)
    !
    call BUILDCONNECTBOUND(N1,N2,N3,h1,h2,h3,l31,lam(1),mu(1))
    !call BUILDCONNECTBOUND2(N1,N2,N3,l31)
    call FINISH()


    do i=1,3*N1*N2
		K(i,:)=0
        K(i,i)=1000
        K(3*N1*N2*N3-i+1,:)=0
        K(3*N1*N2*N3-i+1,3*N1*N2*N3-i+1)=1000
    end do
    R=0
    numn = 5*N1*N2 + (N2-1)*N1 + 1
    R(3*numn-1) = -100*15
    do i=2,N1-1
        numn = 5*N1*N2 + (N2-1)*N1 + i
        R(3*numn-1) = -100*30
    end do
    numn = 5*N1*N2 + (N2-1)*N1 + N1
    R(3*numn-1) = -100*15
    !!R(3*(5)*N1*N2+3*(nv-(N1-1)/2)-1) = -1000
    !
    !open(2,file="Kfem.txt")
    !do i=1,3*N1*N2*N3
    !    write(2,'(3993F16.8)')K(i,:)
    !end do
    !
    !call gesv(K,R)
    !
    !open(2,file="FEM.txt")
    !write(2,*)R
    !close(2)
    !
    !
    !open(2,file="U2_x1=150_x2=300.txt")
    !do i=1,N3
    !    inum = 3*((i-1)*N1*N2+(N2-1)*N1+6-N1)-1
    !    write(2,*)x3(i),R(inum)
    !end do
    !do i=1,N1
    !    R(2*N1*N2-2*i+2)=1000
    !end do
    !
    !
    open(2,file="temp/matK.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(K)/4)
    
    write(2,rec=1)K
    close(2)
    open(2,file="temp/vecRF.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(R)/4)
    
    write(2,rec=1)R
    close(2)
    deallocate(K,R)
    
    call LOADSOLVE()
    
    call SOLVE(3*N1*N2*N3,3*N1*N2)
    
    allocate(y(6*N1*N2))
    !!
    open(2,file="U2_x1=150_x2=270.txt")
    open(12,file="S33_x1=150_x2=270.txt")
    open(13,file="resDCFEM.txt")
    numi = 6
    numj = 2
    inum = (N2-numj)*N1+numi
    h=10
    do i=1,60
        x=(i-1)*5
        if(abs(x-150)<1) then
            call calcY(x,-1,y)            
            call get_stress(N1,N2,numi,N2-numj+1,y,h1,h2,lam(1),mu(1),eps,stress)
            write(12,'(6F15.7)')x,stress(1,1),stress(2,2),stress(3,1),stress(3,2),stress(3,3)
        end if
        call calcY(x,1,y)
        write(13,'(726F10.5)')y
        write(2,'(4F15.7)')x,y(3*(inum-1)+1),y(3*(inum-1)+2),y(3*(inum-1)+3)
        call get_stress(N1,N2,numi,N2-numj+1,y,h1,h2,lam(1),mu(1),eps,stress)
        write(12,'(6F15.7)')x,stress(1,1),stress(2,2),stress(3,1),stress(3,2),stress(3,3)
    end do
    
    x=300
    call calcY(x,-1,y)
    write(13,'(726F10.5)')y
    write(2,'(4F15.7)')x,y(3*(inum-1)+1),y(3*(inum-1)+2),y(3*(inum-1)+3)
    call get_stress(N1,N2,numi,N2-numj+1,y,h1,h2,lam(1),mu(1),eps,stress)
    write(12,'(6F15.7)')x,stress(1,1),stress(2,2),stress(3,1),stress(3,2),stress(3,3)
    
    do i=1,N3
        inum = (i-1)*N1*N2+(N2-numj)*N1+numi
        inn = 2
    
        res1 = getSol(inum,1)
        res2 = getSol(inum,2)
        res3 = getSol(inum,3)
        !res = R(3*(inum-1)+inn)
        write(2,'(4F15.7)')x3(i),res1,res2,res3
    end do
    
    close(2)
    deallocate(y)
    allocate(y(3*N1*N2*N3)) 
    call copy_solve(N1,N2,N3,y)

    do i=1,N3
        
        call get_stress_node(N1,N2,N3,h1,h2,h3,y,numi,N2-numj+1,i,lam(1),mu(1),eps,stress)
        write(12,'(6F15.7)')x3(i),stress(1,1),stress(2,2),stress(3,1),stress(3,2),stress(3,3)
    end do
    
    close(12)
    close(13)
end program
    
    
