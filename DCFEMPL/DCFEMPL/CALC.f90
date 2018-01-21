subroutine SORTEIG()
    use SOLVER
    complex(8)::c=0

    do i=1,n
        do j=i,n
            if(abs(real(e2(j)-e(i)))<epse .and. abs(aimag(e2(j)-e(i)))<epse) then
                c=e2(j)
                e2(j)=e2(i)
                e2(i)=c
                do k=1,n
                    c=T1(k,i)
                    T1(k,i)=T1(k,j)
                    T1(k,j)=c
                end do
                exit
            end if
        end do
        c=0
    end do
    
    !j=n
    !i=1
    !do while(i<=j)
    !    if(abs(real(e(i)))<eps .and. abs(aimag(e(i)))<eps) then
    !        c=e(j)
    !        e(j)=e(i)
    !        e(i)=c
    !        do k=1,n
    !            c=T(k,i)
    !            T(k,i)=T(k,j)
    !            T(k,j)=c
    !            c=T1(k,i)
    !            T1(k,i)=T1(k,j)
    !            T1(k,j)=c
    !        end do
    !        j=j-1
    !        i=i-1
    !    end if
    !    i=i+1
    !end do
	do i=n,2,-1
		do j=1,i-1
			if(abs(real(e(j)))<abs(real(e(j+1)))) then
				c=e(j)
				e(j)=e(j+1)
				e(j+1)=c
				do k=1,n
					c=T(k,j)
					T(k,j)=T(k,j+1)
					T(k,j+1)=c
					c=T1(k,j)
					T1(k,j)=T1(k,j+1)
					T1(k,j+1)=c
				end do
			end if
		end do
    end do    
    l2 = n-12
      
    
    !!отделение нулей и перенос их в конец
    !j=n
    !i=1
    !do while(i<=j)
    !    if(abs(real(e(i)))<eps .and. abs(aimag(e(i)))<eps) then
    !        c=e(i)
    !        e(i)=e(j)
    !        e(j)=c
    !        do k=1,n
    !            c=T(k,i)
    !            T(k,i)=T(k,j)
    !            T(k,j)=c
    !            c=T1(k,i)
    !            T1(k,i)=T1(k,j)
    !            T1(k,j)=c
    !        end do
    !        j=j-1
    !        i=i-1
    !    end if
    !    i=i+1
    !end do
    !
    !l2=j
    k=1
    
    do i=1,l2
        p=0
        do m=1,l2
            if(abs(real(e(m)-e(i)))<epse .and. abs(aimag(e(m)-e(i)))<epse) then
                p=p+1
            end if
        end do
        if(p==1 .and. i>k) then
            c=e(i)
            e(i)=e(k)
            e(k)=c
            do m=1,n
                c=T(m,i)
                T(m,i)=T(m,k)
                T(m,k)=c
                c=T1(m,i)
                T1(m,i)=T1(m,k)
                T1(m,k)=c
            end do
        end if
        if(p==1) k=k+1
    end do
    
    l1=k-1
    
    do i=l2,l1+1,-1
        do j=l1+1,i-1
            if(real(e(j))<real(e(j+1))) then
                c=e(j)
                e(j)=e(j+1)
                e(j+1)=c
                do m=1,n
                    c=t(m,j)
                    t(m,j)=t(m,j+1)
                    t(m,j+1)=c
                    c=t1(m,j)
                    t1(m,j)=t1(m,j+1)
                    t1(m,j+1)=c
                end do
            end if
        end do
    end do          
      
end subroutine


subroutine NORMT()
    use SOLVER
    use lapack95
    use f95_precision
    complex(8)::c
    integer(4),dimension(:),allocatable :: ipiv
    
    allocate(ctemp(l2-l1,l2-l1),ipiv(l2-l1))
    
    do i=1,l1
        c=0
        do k=1,n
            c=c+T(k,i)*T1(i,k)
        end do
        !if(abs(real(c))<epse .and. abs(aimag(c))<epse) then
        !    c=0
        !else
            c=1.0D+00/c
        !end if
        
        do j=1,n
            T1(i,j)=T1(i,j)*c
        end do
    end do
    
    do i=1,l2-l1
        do j=1,l2-l1
            ctemp(i,j)=0
            do k=1,n
                ctemp(i,j)=ctemp(i,j)+T1(l1+i,k)*T(k,l1+j)
            end do
        end do
    end do
    
    if(l2-l1>0) then
        call getrf(ctemp,ipiv)
        call getri(ctemp,ipiv)
    end if
    
    do i=1,l2-l1
        do j=1,n
            T(i,j)=0
            do k=1,l2-l1
                T(i,j)=T(i,j)+ctemp(i,k)*T1(k+l1,j)
            end do
        end do
    end do
    do i=1,l2-l1
        do j=1,n
            T1(i+l1,j)=T(i,j)
        end do
    end do
    deallocate(ipiv)
end subroutine


subroutine CALCP()
    use SOLVER
    use blas95
    use lapack95
    use f95_precision
    complex(8)::c
    
    allocate(P2(n,n),ctemp2(n,n))
    
    P2=0
    ctemp2=0
    call zgemm('n','n',n,n,l1,1.0D+00,T,n,T1,n,0.0D+00,P2,n)

    !do i=1,n
    !    do j=1,n
    !        P2(i,j)=0
    !          
    !        do k=1,l1
    !            c=0
    !            do m=1,n
    !                c=c+T(m,k)*T1(k,m)
    !            end do
    !            P2(i,j)=P2(i,j)+c*T(i,k)*T1(k,j)
    !        end do
    !    end do
    !end do
    
    flag=0
    call zgemm('n','n',n,n,n,1.0D+00,P2,n,P2,n,0.0D+00,ctemp2,n)
    do i=1,n
        do j=1,n
            if(abs(ctemp2(i,j)-P2(i,j))>1e-3) then
                write(FILELOG,*)"P11 != P11^2 ",abs(ctemp2(i,j)-P2(i,j)),abs(P2(i,j))
            end if
        end do
    end do
    !flag=0
    !do i=1,n
    !    do j=1,n
    !        ctemp2(i,j)=0
    !        do k=1,n
    !            ctemp2(i,j)=ctemp2(i,j)+P2(i,k)*P2(k,j)
    !        end do
    !        if(abs(ctemp2(i,j)-P2(i,j))>1e-3) then
    !            write(FILELOG,*)"P11 != P11^2 ",abs(ctemp2(i,j)-P2(i,j))
    !        end if
    !    end do
    !end do

    
    if(l2-l1>0) then

        
        do i=1,n
            do j=1,n
                do k=1,l2-l1
                    P2(i,j)=P2(i,j)+T(i,l1+k)*T1(l1+k,j)
                end do
            end do
        end do
    end if
    
    P2=-P2
    
    forall(i=1:n) P2(i,i)=1+P2(i,i)
    
    deallocate(ctemp,ctemp2)
    allocate(ctemp(n,n))
    flag=0
    do i=1,n
        do j=1,n
            ctemp(i,j)=0
            do k=1,n
                ctemp(i,j)=ctemp(i,j)+P2(i,k)*P2(k,j)
            end do
            if(abs(ctemp(i,j)-P2(i,j))>epse) then
                flag=1
            end if
        end do
    end do
    if(flag==1) then
        write(FILELOG,*)"P2 != P2^2"
    end if
    deallocate(ctemp)
    
    T=a
    T1=0.0D+00
    !call gemm(P2,T,T1)
    call zgemm('n','n',n,n,n,1.0D+00,P2,n,T,n,0.0D+00,T1,n)
    T=T1
end subroutine