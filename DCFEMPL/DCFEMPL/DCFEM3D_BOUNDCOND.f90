subroutine LDISPLBOUNDCOND(k,xb,n)
    use DCFEM3D
    integer(4)::k,n,np
    real(8)::xb
    real(8),dimension(:,:),allocatable :: br
    real(8),dimension(:),allocatable :: g
    
    character(20)::s=""
    write(s,'(I4.4)')k
    

    write(FILEXB,rec=k) xb
    numeq=n
    write(FILEJOB,*) 6*numeq
    
    allocate(br(3*numeq,6*numeq),g(3*numeq))
    
    open(FILEMATBP,file="temp/matB+"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(br)/4)
    open(FILEMATBM,file="temp/matB-"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(br)/4)
    open(FILEVECG,file="temp/vecG"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(g)/4)

    br=0
    do i=1,numeq*3
        br(i,i)=1
    enddo
    
    
    write(FILEMATBP,rec=1)br
    br=0
    write(FILEMATBM,rec=1)br
    g=0
    write(FILEVECG,rec=1)g
    
    close(FILEMATBP)
    close(FILEMATBM)
    close(FILESBK)
    close(FILEVECG)
    deallocate(br,g)
end subroutine


subroutine RDISPLBOUNDCOND(k [value], xb [value], n [value], np [value], u1 [value], u2 [value])
    use DCFEM3D
    !DEC$ ATTRIBUTES DLLEXPORT :: RDISPLBOUNDCOND
    integer(4)::k,n,np
    real(8)::xb
    real(8),dimension(:,:),allocatable :: br
    real(8),dimension(:),allocatable :: g
    real(8)::u1,u2
    
    character(20)::s=""
    
    write(FILEXB,rec=k) xb
    
    numeq=n
    write(FILEJOB,*) 6*numeq
    
    allocate(br(6*numeq,6*numeq),g(6*numeq))
    s=""
    write(s,'(I4.4)')k

    open(FILEMATBP,file="temp/matB+"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(br)/4)
    open(FILEMATBM,file="temp/matB-"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(br)/4)
    open(FILEVECG,file="temp/vecG"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(g)/4)
    

    
    br=0
    do i=1,numeq*3
        br(3*numeq+i,i)=1
    enddo

    
    write(FILEMATBM,rec=1)br
    br=0
    write(FILEMATBP,rec=1)br
    g=0
    if(abs(u1)>epsz) then
        do i=1,3*n-1,2
            g(3*np+i)=u1
        end do
    end if
    if(abs(u2)>epsz) then
        do i=2,3*n,2
            g(3*np+i)=u2
        end do
    end if
    
    write(FILEVECG,rec=1)g
    
    close(FILEMATBP)
    close(FILEMATBM)
    close(FILESBK)
    close(FILEVECG)
    deallocate(br,g)
    end subroutine

    
subroutine BUILDLEFTBOUND(nv,Xb)

    use DCFEM3D
    integer(4)::nv
    real(8)::Xb
    
    real(8),dimension(:,:),allocatable :: b
    real(8),dimension(:),allocatable :: g
        
    character(20)::s=""
    write(s,'(I4.4)')1
    
    write(FILEXB,rec=1) Xb
    write(FILEJOB,*) 6*nv
    
    allocate(b(3*nv,6*nv),g(3*nv))
    
    open(FILEMATBP,file="temp/matB+"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(b)/4)
    open(FILEVECG,file="temp/vecG"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(g)/4)
    
    b = 0
    g = 0
    do i=1,3*nv
        b(i,i) = 1000
    end do
    
    write(FILEMATBP,rec=1)b
    write(FILEVECG,rec=1)g
    
    close(FILEMATBP)
    close(FILEVECG)
    deallocate(b,g)
end subroutine
    
       
subroutine BUILDCONNECTBOUND(N1,N2,N3,h1,h2,h3,Xb,lam,mu)

    use DCFEM3D
    integer(4)::N1,N2,N3
    real(8)::h1,h2,h3,Xb
    real(8):: lam,mu
    
    real(8),dimension(:,:),allocatable :: bm,bp
    real(8),dimension(:),allocatable :: g
    
    
    real(8),dimension(:),allocatable :: y
    real(8),dimension(3,3)::eps,stress
    character(20)::s=""
    write(s,'(I4.4)')2
    
    allocate(bm(6*N1*N2,6*N1*N2),g(6*N1*N2))
    
    bm = 0
    g = 0
    
    
    write(FILEXB,rec=2) xb
    write(FILEJOB,*) 6*N1*N2
    
    open(FILEMATBM,file="temp/matB-"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(bm)/4)
    
    open(FILEVECG,file="temp/vecG"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(g)/4)    
    
    do i=1,3*N1*N2
        bm(i,i) = 1000
    end do
    allocate(y(6*N1*N2))

    do ig=1,N1
        do jg=1,N2
            do i=1,6*N1*N2
                y=0
                y(i)=1
                call get_stress(N1,N2,ig,jg,y,h1,h2,lam,mu,eps,stress)
                bm(3*N1*N2+(jg-1)*N1+ig,i)=stress(3,1)
                bm(4*N1*N2+(jg-1)*N1+ig,i)=stress(3,2)
                bm(5*N1*N2+(jg-1)*N1+ig,i)=stress(3,3)
            end do
        end do
    end do
    
    !!!--------------------------
    !!!Фрагмент, приравнивающий касательные напряжения

    write(FILEMATBM,rec=1)bm
    write(FILEVECG,rec=1)g
    
    deallocate(bm,g)
    
    allocate(bp(6*N1*N2,3*N1*N2*N3))
    open(FILEMATBP,file="temp/matB+"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(bp)/4)
    
    bp=0
    do i=1,3*N1*N2
        bp(i,i) = 1000
    enddo
    deallocate(y)
    allocate(y(3*N1*N2*N3))
    do ig=1,N1
        do jg=1,N2
            
            do i=1,6*N1*N2
                    y=0
                    
                    y(i)=1
                    call get_stress_node(N1,N2,N3,h1,h2,h3,y,ig,jg,1,lam,mu,eps,stress)
                    bp(3*N1*N2+(jg-1)*N1+ig,i)=stress(3,1)
                    bp(4*N1*N2+(jg-1)*N1+ig,i)=stress(3,2)
                    bp(5*N1*N2+(jg-1)*N1+ig,i)=stress(3,3)

            end do

        end do
    end do

    write(FILEMATBP,rec=1)bp
    
    deallocate(bp)
    
    close(FILEMATBP)
    close(FILEVECG)
    close(FILEMATBM)
    end subroutine
    

    subroutine BUILDCONNECTBOUND2(N1,N2,N3,xb)

    use DCFEM3D
    integer(4)::N1,N2,N3
    real(8)::xb
    real(8),dimension(:,:),allocatable :: bm,bp
    real(8),dimension(:),allocatable :: g
    
    character(20)::s=""
    write(s,'(I4.4)')2
    
    allocate(bm(6*N1*N2,6*N1*N2),g(6*N1*N2))
    
    bm = 0
    g = 0
    
    
    write(FILEXB,rec=2) xb
    write(FILEJOB,*) 6*N1*N2
    
    open(FILEMATBM,file="temp/matB-"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(bm)/4)
    
    open(FILEVECG,file="temp/vecG"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(g)/4)    
    
    do i=1,3*N1*N2
        bm(i,i) = 1000
    end do
    
     write(FILEMATBM,rec=1)bm
    write(FILEVECG,rec=1)g
    
    deallocate(bm,g)
    
    allocate(bp(6*N1*N2,3*N1*N2*N3))
    open(FILEMATBP,file="temp/matB+"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(bp)/4)
    
    bp=0
    do i=1,3*N1*N2
        bp(3*N1*N2+i,i) = 1000
    enddo
      
    write(FILEMATBP,rec=1)bp
    
    deallocate(bp)
    
    close(FILEMATBP)
    close(FILEVECG)
    close(FILEMATBM)
    end subroutine