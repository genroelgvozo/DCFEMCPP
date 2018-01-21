
module DCFEM3D
    real(8),parameter :: EPS_0 = 1e-5	!Точность сравнения с нулем
    real(8),parameter :: EPS_EQ = 1e-5	!Точность сравнение друг с другом
	real(8),parameter :: epsz = 1e-5	!Точность сравнения с нулем
    real(8),parameter :: epsq = 1e-5	!Точность сравнение друг с другом
	integer(4),parameter :: LOGFILE = 13 !Файл логов (локальный)
	integer(4),parameter :: OUTFILE = 9 !Глобальный файл логов
    integer(4),parameter :: FILEJOB =1, FILEMATA = 2 , FILEVECR = 3, FILEMATBP = 4, FILEMATBM = 10, FILEVECG = 7, FILEXB = 8, FILEXF =11, FILESBK =12
    integer(2) :: isOpen=0
	
	real(4)::time1,time2 !> For calculate time
end module