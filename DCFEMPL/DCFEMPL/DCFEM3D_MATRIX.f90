
	
subroutine CLOSEFILES()
    use DCFEM3D
    close(FILEJOB)
    close(FILEXB)
    close(OUTFILE)
end subroutine

subroutine FINISH()
    use DCFEM3D
    

    call CLOSEFILES()
    isOpen=0
end subroutine

subroutine INIT() 
    
    use DCFEM3D


    if(isOpen) then
        call FINISH()
    end if
    call system("mkdir temp")
	call system("mkdir logs")
    
    
    open(filejob,file="temp/solve.txt")
    open(outfile,file="logs/dcfem.log")
    
    open(filexb,file="temp/xb.bin", form = "unformatted", access = "direct", recl = sizeof(epsz)/4)
    write(outfile,*) "инициализация.."


    write(filejob,*) 1

    
    isOpen=1
endsubroutine



subroutine invertKvv(n,a)
	
    use DCFEM3D
    use lapack95

    integer(4)::n
    real(8),dimension(n,n)::a
	
    integer(4),dimension(:),allocatable :: flags,ipiv

    allocate(flags(n),ipiv(n))
    flags=0
    

    
    do i=1,n
        flag1=0
        flag2=0
        do j=1,n
            if(abs(a(i,j))>epsz) flag1=1
            if(abs(a(j,i))>epsz) flag2=1
        end do
        if(flag1==0 .and. flag2==0) then
            a(i,i)=1
            flags(i)=1
        end if
    end do
    
    write(OUTFILE,*)"Procedure invertKvv"
    
    call getrf(a,ipiv)
    write(OUTFILE,*)"Procedure getrf"
    call getri(a,ipiv)
    

    do i=1,n
        if(flags(i)==1) a(i,i)=0
    end do
    
    deallocate(flags)
end subroutine




subroutine CHECK_MESH(nv [value], ne [value], coords, elems)

	use DCFEM3D

	integer(4)::nv
    integer(4)::ne
    real(8),dimension(2,nv)::coords
	integer(4),dimension(4,ne)::elems
	
	open(LOGFILE,file="log/CHECK_MESH.log")
	
	write(LOGFILE,*)""
	write(LOGFILE,*)"---Procedure of checking mesh---"
	write(LOGFILE,*)""
	
	write(LOGFILE,*)"	Number of vertex	:",nv
	write(LOGFILE,*)"	Number of elements	:",ne
	
	write(LOGFILE,*)""
	write(LOGFILE,*)"---Checking of intersection elements---"
	
	do i=1,ne
	end do
	
	write(LOGFILE,*)"	Done : no elements are intersect"
	
	close(LOGFILE)
	
end subroutine



subroutine BUILDSDE2(k1 , nv , ne , coords, elems, lam, mu, nf , xf, f1, f2, f3, err)
    
    use lapack95
    use blas95
    
    use DCFEM3D
    use f95_precision
    integer(4)::k1
    integer(4)::nv
    integer(4)::ne
    real(8),dimension(2,nv)::coords
	integer(4),dimension(4,ne)::elems
	real(8),dimension(ne)::lam,mu
	integer(4)::nf
	real(8),dimension(nf)::xf
    real(8),dimension(nf,nv)::f1,f2,f3
	integer(4)::err
	
	!LOCAL VARIABLES
	real(8),dimension(:,:),allocatable::ki,kuu,kuv,kvu,kvv	!Локальная матрица жесткости и подматрицы
	real(8),dimension(:,:),allocatable::x					!Coordinates of element's vertexes
	real(8),dimension(:,:),allocatable::a,akv					!Global matrix of system
	real(8),dimension(:),allocatable::Rvec
	character(20)::s=""     
	
	write(FILEJOB,*)6*nv
    write(FILEJOB,*)nf

	!Checking mesh    
	!call CHECK_MESH(nv,ne,coords,elems)
	

	allocate(a(6*nv,6*nv),Rvec(6*nv))
	
	write(OUTFILE,*)"Matrix A is allocated"
		
    write(s,'(I4.4)')k1
    allocate(x(4,2))
    open(FILEMATA,file="temp/matA"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(a)/4)
    open(FILEVECR,file="temp/vecR"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(Rvec)/4)
    open(FILEXF,file="temp/xf"//trim(s)//".bin", FORM = "UNFORMATTED", ACCESS = "DIRECT", RECL = sizeof(x(1,1))/4)
    	
	!Form matrix for all quad elements
	allocate(ki(24,24),kuu(12,12),kuv(12,12),kvu(12,12),kvv(12,12))

	

	numi=3*nv
	
	do i=1,ne


		
		if(elems(4,i) > 0) then
			
			ki = 0
			
			do j=1,4
				x(j,1)=coords(1,elems(j,i))
				x(j,2)=coords(2,elems(j,i))
			end do
			
			call get_kelem_quad(x, lam(i), mu(i), ki)
		
			do ik=1,4
                do jk=1,4
                    kuu((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)=ki((ik-1)*6+1:(ik-1)*6+3,(jk-1)*6+1:(jk-1)*6+3)
                    kuv((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)=ki((ik-1)*6+1:(ik-1)*6+3,(jk-1)*6+4:(jk-1)*6+6)
                    kvu((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)=ki((ik-1)*6+4:(ik-1)*6+6,(jk-1)*6+1:(jk-1)*6+3)
                    kvv((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)=ki((ik-1)*6+4:(ik-1)*6+6,(jk-1)*6+4:(jk-1)*6+6)
                end do
            end do
			
			!write(FILELOG,*)kuu
			
			do ik=1,4
				do jk=1,4
					a((elems(ik,i)-1)*3+1:(elems(ik,i)-1)*3+3,(elems(jk,i)-1)*3+1:(elems(jk,i)-1)*3+3)=a((elems(ik,i)-1)*3+1:(elems(ik,i)-1)*3+3,(elems(jk,i)-1)*3+1:(elems(jk,i)-1)*3+3)+kuu((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)
					a((elems(ik,i)-1)*3+1:(elems(ik,i)-1)*3+3,numi+(elems(jk,i)-1)*3+1:numi+(elems(jk,i)-1)*3+3)=a((elems(ik,i)-1)*3+1:(elems(ik,i)-1)*3+3,numi+(elems(jk,i)-1)*3+1:numi+(elems(jk,i)-1)*3+3)+kuv((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)
					a(numi+(elems(ik,i)-1)*3+1:numi+(elems(ik,i)-1)*3+3,(elems(jk,i)-1)*3+1:(elems(jk,i)-1)*3+3)=a(numi+(elems(ik,i)-1)*3+1:numi+(elems(ik,i)-1)*3+3,(elems(jk,i)-1)*3+1:(elems(jk,i)-1)*3+3)+kvu((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)
					a(numi+(elems(ik,i)-1)*3+1:numi+(elems(ik,i)-1)*3+3,numi+(elems(jk,i)-1)*3+1:numi+(elems(jk,i)-1)*3+3)=a(numi+(elems(ik,i)-1)*3+1:numi+(elems(ik,i)-1)*3+3,numi+(elems(jk,i)-1)*3+1:numi+(elems(jk,i)-1)*3+3)+kvv((ik-1)*3+1:ik*3,(jk-1)*3+1:jk*3)
					
				end do
			end do
			
		end if
	end do

	deallocate(kuu,kuv,kvu,kvv,ki,x)


	do i=1,numi
		do j=1,numi
			a(i,numi+j)=a(i,numi+j)-a(numi+i,j)
		end do
	end do
	write(OUTFILE,*)"Matrix if Formed 1"
    a(numi+1:2*numi,1:numi) = a(1:numi,1:numi)
	write(OUTFILE,*)"Matrix is Formed 2"
    a(1:numi,1:numi) = a(numi+1:2*numi,numi+1:2*numi)
    a(numi+1:2*numi,numi+1:2*numi) = a(1:numi,numi+1:2*numi)

    
    
	write(OUTFILE,*)"Matrix is Formed"
	
    allocate(akv(3*nv,3*nv))
    call invertKvv(3*nv,a(1:3*nv,1:3*nv))

	write(OUTFILE,*)"Matrix Kvv is inverted"
    do i=1,nf
        
        Rvec=0
        
        write(FILEXF,rec=i)xf(i)

        do j=1,nv
            Rvec(numi+3*j-2)=f1(i,j)
            Rvec(numi+3*j-1)=f2(i,j)
            Rvec(numi+3*j)=f3(i,j)
        end do
        do j=1,numi
            do k=1,numi
                Rvec(j)=Rvec(j)+a(j,k)*Rvec(numi+k)
            end do
        end do
        
        Rvec(numi+1:2*numi) = -Rvec(1:numi)
        Rvec(1:numi)=0
        
        write(FILEVECR,rec=i)Rvec
    end do
  
    
    !call mrrrr(a(1:numi,1:numi),a(numi+1:2*numi,1:numi),a(1:numi,numi+1:2*numi))
    !open(12,file="MAta3.txt")
    !write(12,*)a
    !close(12)
    !do i=1,numi
    !    do j=1,numi
    !        akv(i,j) = a(numi+i,j)
    !    end do
    !end do
    !!call dgemm('n','n',numi,numi,numi,1.0D+00,a,6*nv,akv,numi,0.0D+00,akv,numi)
    !call gemm(a(1:numi,1:numi),akv,akv)
    !open(12,file="tempakv.txt")
    !write(12,*)akv
    !a(1:numi,numi+1:2*numi) = akv
    !
    
    a(1:numi,numi+1:2*numi)=matmul(a(1:numi,1:numi),a(numi+1:2*numi,1:numi))
    a(numi+1:2*numi,1:numi)=a(1:numi,numi+1:2*numi)
    

    !call mrrrr(a(1:numi,1:numi),a(numi+1:2*numi,numi+1:2*numi),a(1:numi,numi+1:2*numi))
    

    do i=1,numi
        do j=1,numi
            akv(i,j) = a(numi+i,numi+j)
        end do
    end do
    !call dgemm('n','n',numi,numi,numi,1.0D+00,a,6*nv,akv,numi,0.0D+00,akv,numi)
    !call gemm(a(1:numi,1:numi),akv,akv)
    !write(12,*)"-----------------------------"
    !write(12,*)akv
    !close(12)
    !a(1:numi,numi+1:2*numi) = akv
    !
    a(1:numi,numi+1:2*numi)=matmul(a(1:numi,1:numi),a(numi+1:2*numi,numi+1:2*numi))
    do i=1,numi
        do j=1,numi
            a(numi+i,numi+j)=a(i,numi+j)
        end do
    end do
    

    a(1:numi,:) = 0
    do j=1,numi
        a(j,numi+j)=1
    end do
	
    open(12,file="Mata.txt")
    write(12,*)a
    close(12)

    write(FILEMATA,rec=1)a
	
    close(FILEMATA)
    close(FILEVECR)
    close(FILEXF)
    deallocate(a,Rvec)
end subroutine
