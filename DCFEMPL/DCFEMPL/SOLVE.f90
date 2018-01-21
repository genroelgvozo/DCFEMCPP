module SOLVER
    integer(2)::isOpen=0,isCoeff=0
    integer(2),parameter :: FILEJOB = 1,FILEMATA=2, FILETEMP=9,FILETEMPE=10,FILELOG=14,FILEVECR=3
    integer(2),parameter :: FILETEMPL = 13, FILEMATBP = 4, FILEMATBM = 15, FILEVECG = 12, FILEXB = 8, FILEXF=16
    
    real(8),parameter :: epse = 1e-5
    real(8),parameter :: eps = 0.005
    
    integer(4)::nk
    
    integer(4)::n,nf
    real(8),dimension(:),allocatable :: gk,gk2
    real(8),dimension(:,:),allocatable :: a
    real(8),dimension(:),allocatable :: er,ei,R
    real(8),dimension(:,:),allocatable:: Tmkl,Tmkl1
    complex(8),dimension(:),allocatable :: e,e2
    complex(8),dimension(:,:),allocatable :: T,T1,P2
    real(8),dimension(:,:),allocatable :: B
    complex(8),dimension(:,:),allocatable :: ctemp,ctemp2
    real(8),dimension(:,:),allocatable :: Kf,Km,Coeffs
    complex(8),dimension(:),allocatable ::  Sv
    
    integer(4):: l1,l2
    
    real(8)::ct
    complex(8)::cct
    
    
end module
    
subroutine PROCESS_MATRIX()
    use SOLVER
        use blas95
    use lapack95
    use f95_precision
    character(12)::s
    integer(4)::flagcompl = 0
    real(8)::dott
    read(FILEJOB,*)n
    read(FILEJOB,*)nf
    
    write(s,'(I4.4)') 1
    
    allocate(a(n,n),er(n),ei(n),Tmkl(n,n),Tmkl1(n,n))
    
    open(FILEMATA,file="temp/matA"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(a)/4)
    open(FILETEMPL,file="temp/templ"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(l1)/4)
    open(FILEXB,file="temp/xb.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(eps)/4)
    
    read(FILEMATA,rec=1)a
    
    
	write(FILELOG,*)"---Подсчет собственных значений---"
    call cpu_time(time1)
    call geev(a,er,ei, Tmkl1,Tmkl)
    !call evcrg(a,e,T)
    write(FILELOG,*)"---Матрица T вычислена--"

    allocate(T(n,n),e(n),e2(n))
    open(FILETEMP,file="temp/temp"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(T)/4)
    open(FILETEMPE,file="temp/tempe"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(e)/4)
    
    
    do i=1,n
        e(i)=CMPLX(er(i),ei(i))
    end do
    
    do i=1,n
        if(flagcompl == 1) then
            flagcompl =0
            cycle
        end if
        if(abs(ei(i))>1e-5) then
            do j=1,n
                T(j,i)=CMPLX(Tmkl(j,i),Tmkl(j,i+1))
                T(j,i+1)=CMPLX(Tmkl(j,i),-Tmkl(j,i+1))
            end do
            flagcompl = 1
        else
            T(:,i)=Tmkl(:,i)
        end if
    end do
    write(FILETEMP,rec=1)T
    read(FILEMATA,rec=1)a
    
    write(FILELOG,*)"Транспонирование матрицы A"
    do i1=1,n
        do j1=1,i1-1,1
           ct=a(i1,j1)
            a(i1,j1)=a(j1,i1)
            a(j1,i1)=ct
        end do
    end do
        
    write(FILELOG,*)"---Подсчет матрицы T~--"
    call geev(a,er,ei,Tmkl1, Tmkl)
    
    do i=1,n
        e2(i)=CMPLX(er(i),ei(i))
    end do
    
    flagcompl = 0
    do i=1,n
        if(flagcompl == 1) then
            flagcompl =0
            cycle
        end if
        if(abs(ei(i))>1e-5) then
            do j=1,n
                T(j,i)=CMPLX(Tmkl(j,i),Tmkl(j,i+1))
                T(j,i+1)=CMPLX(Tmkl(j,i),-Tmkl(j,i+1))
            end do
            flagcompl = 1
        else
            T(:,i)=Tmkl(:,i)
        end if
    end do

    !call evcrg(a,e2,T)
    !call cpu_time(time2)
    write(FILELOG,'("Время нахождения собственных значений:",F10.1," millisecs")')(time2-time1)*1000

    read(FILEMATA,rec=1)a
    !a=transpose(a)
    write(FILETEMP,rec=2)T
    
    deallocate(Tmkl,er,ei,Tmkl1)
    
    allocate(T1(n,n))
    
    T1=T
    read(FILETEMP,rec=1)T
    
    call SORTEIG()

    do i1=1,n
        do j1=1,i1-1,1
            cct=T1(i1,j1)
            T1(i1,j1)=T1(j1,i1)
            T1(j1,i1)=cct
        end do
    end do
    write(FILETEMP,rec=1)T
    
    call NORMT()
    
    read(FILETEMP,rec=1)T
    
    do j=1,n
        dott = sum(T(:,j)*T1(j,:))
        write(FILELOG,*)"dot",j,"=",dott
    end do
    !deallocate(ctemp)
    !        iraz=l1
    !    allocate(ctemp(iraz,iraz))
    !    !call mcrcr(T(:,1:l1),T1(1:l1,:),ctemp)
    !    call gemm(T1(1:l1,:),T(:,1:l1),ctemp)
    !    
    !    forall(i=1:iraz) ctemp(i,i)=ctemp(i,i)-1
    !    summ=0
    !    do j=1,iraz
    !        do k=1,iraz
    !            if(summ<abs(ctemp(j,k))) then
    !                summ=abs(ctemp(j,k))
    !            end if
    !        end do
    !    end do
    !    write(FILELOG,*)"max(T*T1 - E)=",summ
    write(FILETEMP,rec=2)T1
    write(FILETEMPE,rec=1)e
    write(FILETEMPL,rec=1)l1
    write(FILETEMPL,rec=2)l2
    

    
    do i=1,n
    write(FILELOG,'(4F15.6)')e(i)
    end do
    write(FILELOG,*)"-------------"
    write(FILELOG,*)l1
    write(FILELOG,*)l2
    
    call CALCP()
    write(FILETEMP,rec=3)P2
    write(FILETEMP,rec=4)T
    do j=2,n-l2+2
        call gemm(T,T1,P2)
        T=P2

        write(FILETEMP,rec=3+j)T
    end do
    allocate(ctemp(n,n),ctemp2(n,n))
        call calcfuncfull(0.0D+00,1,0,ctemp,0)
        call calcfuncfull(0.0D+00,-1,0,ctemp2,0)
        ctemp=ctemp-ctemp2
        forall(i=1:n) ctemp(i,i)=ctemp(i,i)-1
        summ=0
        do j=1,n
            do k=1,n
                summ=summ+abs(ctemp(j,k))
            end do
        end do
        write(FILELOG,*)"sum(E(+0) - E(-0)-E)=",summ
    deallocate(ctemp,ctemp2)
    deallocate(a,e2)
    close(FILEMATA)
    
    end subroutine
    
subroutine calcfuncfull(x,l,s,E1,debug)
    use SOLVER
    real(8),intent(in)::x
    integer(4),intent(in)::l
    integer(4),intent(in)::s
    complex(8),dimension(n,n),intent(inout)::E1
    integer(4)::debug
    integer(4)::i,j,k,p
    real(8)::x1
    
    E1=0
    if(dabs(x)<eps) then
        x1=x+l
    else
        x1=x
    endif
    
    read(FILETEMP,rec=1)T
    read(FILETEMP,rec=2)T1
    
    if(debug==0 .or. debug==1) then
    
    do i=1,n
        do j=1,n
            do k=1,l2
                if( (x1>0 .and. dreal(e(k))<=0) .or. (x1<0 .and. dreal(e(k))>0) ) then
                    E1(i,j)=E1(i,j)+T(i,k)*dsign(1.0,x1) *cdexp(e(k)*x) *T1(k,j) *e(k)**s
                end if
            end do
        end do
    end do
    endif

    if(debug==0 .or. debug==2) then
    tt=1
      if(x1>0) then
 !       tt=0.5
  !  else
!        tt=-0.5
   ! end if
        if(s<=0) then
            read(FILETEMP,rec=3)T
            E1=E1+tt*T*x**(-s)
        end if
        
        x1=x**(1-s)
        
        i=1
        do j=1,1-s
            i=i*j
        end do
        x1=x1/i
        
        read(FILETEMP,rec=4)T
        
        do i=1,n-l2+1
            E1=E1+tt*T*x1
            read(FILETEMP,rec=4+i)T
            x1=x1*x/(i+1-s)
        end do

	endif
    
    end if


end subroutine
    
subroutine svertka(x,l,S1,debug)
    use SOLVER

    real(8)::x
    integer(4)::l
    complex(8),dimension(n)::S1
    integer(4)::debug
    real(8)::xf,xb,xb1
    character(12)::s=""
    write(s,'(I4.4)') 1
    
    S1=0
    
    allocate(gk(n),ctemp(n,n))
    open(FILEVECR,file="temp/vecR"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(gk)/4)
    open(FILEXF,file="temp/xf"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(xf)/4)
    read(FILEXB,rec=1)xb
    read(FILEXB,rec=2)xb1
    
    do i=1,nf
        read(FILEXF,rec=i)xf
        read(FILEVECR,rec=i)gk
        
        if(xf>xb .and. xf<xb1) then
            call calcfuncfull(x-xf,l,0,ctemp,debug)
            
            do j=1,n
                S1(j)=S1(j)+sum(ctemp(j,:)*gk)
            end do
        end if
    end do
    
    deallocate(gk,ctemp)
    close(FILEVECR)
    close(FILEXF)
end subroutine

subroutine SOLVE(Nfem,Nbound)
    
    use SOLVER
    use lapack95
    use blas95
    use f95_precision
    integer(4)::Nfem,Nbound
    real(8)::x,xb,xb1
    x=0.0D+00
    if(isCoeffs) then
        deallocate(Coeffs)
    end if
    print*,"Before Allocating global matrix"
    allocate(Km(n+Nfem,n+Nfem),Kf(Nfem,Nfem),Coeffs(n+Nfem,1),R(Nfem))
    print*,"After Allocating global matrix"
    isCoeffs=1
    
    Km=0
    Coeffs=0
    print*,"Km Coeffs zero"
    
    read(FILEXB,rec=1)xb
    read(FILEXB,rec=2)xb1
    Kf = 0
    print*,"Reading Xb"
    
    open(2,file="temp/matK.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(Kf)/4)
    
    
    
    read(2,rec=1)Kf
    print*,"Reading Kfem"
    do i=1,Nfem
        do j=1,Nfem
            Km(n+i,n+j)=Kf(i,j)
        end do
    end do
        

    deallocate(Kf)
    
    print*,"KFem OK"
    
    close(2)
    
    
    open(2,file="temp/vecRF.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(R)/4)
    read(2,rec=1)R
    Coeffs(n+1:n+Nfem,1)=R
    deallocate(R)
    close(2)
    
    allocate(b(n,Nfem),Sv(n),R(n))
    
    open(FILEMATBP,file="temp/matB+0002.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(b)/4)
    
    read(FILEMATBP,rec=1)b
    
    
    
    Km(1:n+Nbound,:)=0
    Km(Nbound+1:Nbound+n,n+1:n+Nfem)=-b
    print*,"B+ Fem OK"

    deallocate(b)
    close(FILEMATBP)
    
     
    allocate(b(n,n))
    open(FILEMATBP,file="temp/matB-0002.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(b)/4)
    open(FILEVECG,file="temp/vecG0002.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(R)/4)    
    read(FILEMATBP,rec=1)b
    read(FILEVECG,rec=1)R
    
    allocate(ctemp(n,n),ctemp2(n,n),gk2(n))
    
    call calcfuncfull(xb1-xb,-1,0,ctemp,0)
    
    call calcfuncfull(0.0D+00,-1,0,ctemp2,0)
    
    ctemp = ctemp - ctemp2
    
    T = b

    ctemp2=0
    
    call gemm(T,ctemp,ctemp2)
    

    
    
    Km(Nbound+1:Nbound+n,1:n) = real(ctemp2)
    
    deallocate(ctemp)
    
    call svertka(xb1,-1,Sv,0)
    allocate(ctemp(n,n))
    gk2 = 0
    R= real(Sv)
    call gemv(b,R,gk2)
    
    !do i=1,2*N1
    !    Km(2*N1+i,1:4*N1)=real(ctemp2(i,:))
    !    Coeffs(2*N1+i,1)=gk2(i)
    !end do
    !
    !do i=1,N1
    !    j=2*(i-1)*N2+1
    !    Km(4*N1+j,1:4*N1)=real(ctemp2(2*N1+i,:))
    !    Coeffs(4*N1+j,1)=gk2(2*N1+i)
    !    Km(4*N1+j+1,1:4*N1)=real(ctemp2(3*N1+i,:))
    !    Coeffs(4*N1+j+1,1)=gk2(3*N1+i)
    !end do   
    !
    Coeffs(Nbound+1:Nbound+n,1)=-gk2
    
    print*,"B- DCFEM OK"
    close(FILEMATB)
    close(FILEVECG)
    
    deallocate(b,R,gk2)
    allocate(b(Nbound,n),R(Nbound))
    
    open(FILEMATBP,file="temp/matB+0001.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(b)/4)
    open(FILEVECG,file="temp/vecG0001.bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(R)/4)    
    read(FILEMATBP,rec=1)b
    read(FILEVECG,rec=1)R
    
    
    call calcfuncfull(0.0D+00,1,0,ctemp,0)
    
    call calcfuncfull(xb-xb1,1,0,ctemp2,0)
    
    ctemp = ctemp - ctemp2
    
    deallocate(T,ctemp2)
    
    allocate(T(Nbound,n),ctemp2(Nbound,n))
    
    T=b
    
    ctemp2=0
    call gemm(T,ctemp,ctemp2)
    deallocate(T)
    allocate(T(n,n))
    Km(1:Nbound,1:n) = real(ctemp2)
    
    deallocate(ctemp,R)

    call svertka(xb,1,Sv,0)
    b= -b
    allocate(gk2(n),R(Nbound))

    gk2 = real(Sv)
    R=0
    call gemv(b,gk2,R)
    Coeffs(1:Nbound,1)=R
    deallocate(R,gk2,Sv,b,ctemp2)
    open(2,file="RightVector.txt")
    write(2,*)Coeffs
    close(2)
    open(2,file="Km.txt")
    write(2,*)Km
    close(2)
    
    call gesv(Km,Coeffs)
    
    open(2,file="res.txt")

    do i=1,Nfem/3
        write(2,*)Coeffs(n+3*i-2,1),Coeffs(n+3*i-1,1),Coeffs(n+3*i,1)
    end do
    close(2)
    close(FILEVECG)
end subroutine

    
subroutine LOADSOLVE()
    use SOLVER
    
    if(isOpen==1) then
        close(FILEJOB)
        close(FILELOG)
    end if
       
    open(FILEJOB,file="temp/solve.txt")
	open(FILELOG,file="temp/logSOLVER.txt")
    
    isOpen = 1
    
    
    read(FILEJOB,*) nk
    
    call PROCESS_MATRIX()
    
end subroutine
    
    
subroutine calcY(x,l,y)
    use SOLVER
    real(8)::x
    integer(4)::l
    real(8),dimension(n)::y
    integer(4)::k
    real(8)::xb,xb1
    character(12)::s=""
    
    read(FILEXB,rec=1)xb
    read(FILEXB,rec=2)xb1

    
    write(s,'(I4.4)')1

    
    do i=1,n
        y(i)=0
    end do
    
    allocate(ctemp(n,n))
    
    
    call calcfuncfull(x-xb,l,0,ctemp,0)
    
    do i=1,n
        y(i)=y(i)+real( sum(ctemp(i,:)*Coeffs(1:n,1)))
    end do
    call calcfuncfull(x-xb1,l,0,ctemp,0)
    
    do i=1,n
        y(i)=y(i)-real( sum(ctemp(i,:)*Coeffs(1:n,1)))
    end do
    
    deallocate(ctemp)
    
    allocate(Sv(n))
    call svertka(x,l,Sv,0)
    
    do i=1,n
        y(i)=y(i)+real(Sv(i))
    end do
    deallocate(Sv)
    
    write(FILELOG,*)"---------------------"
    write(FILELOG,*)"Значение в точке:",x
    write(FILELOG,*)y(1:n)
    
end subroutine
    
    
function getSol(i,num)
    use SOLVER
    integer(4)::i,num
    real(8)::getSol
    
    getSol = Coeffs(n+(i-1)*3+num,1)
end function
    
    
subroutine copy_solve(N1,N2,N3,y)
    use SOLVER
    integer(4)::N1,N2,N3
    real(8),dimension(3*N1*N2*N3)::y
    
    y = Coeffs(n+1:n+3*N1*N2*N3,1)
end subroutine
