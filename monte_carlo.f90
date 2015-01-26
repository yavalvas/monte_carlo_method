!вычисление определённого интеграла простейшим
!методом Монте-Карло
!и с исользованием приёмов, уменьшающих дисперсию
!вычислений: выделение главной части;
!сущетвенная выборка; симметризация
!подынтегральной функци. здесь вычисляется
!двойной интеграл от функции (y^2)*exp(x) в интервале x=(0,1), y=(0,1).

PROGRAM MMK2
INTEGER,PARAMETER ::N=100000
REAL,PARAMETER ::I0=1.-2.*EXP(1.)+EXP(2.), I1=EXP(1.)-1. !точное значение интеграла двойного интеграла I0 и одинарного I1
REAL TETA1, TETA2, TETA3, TETA4,X,y,Z,GAMMA1,GAMMA2,GAMMA3,GAMMA4,ksiG1,IMK,AA,DN
OPEN(1,FILE='RESULTS.DAT') !Открытие файла
                        !для записи результатов
                        !вычислений
OPEN(2,FILE='PRECISIONTETA1.csv') !Открытие файла для записи значений точности при различных N
OPEN(3,FILE='PRECISIONGAMMA1.csv') !Открытие файла для записи значений точности при различных N
OPEN(4,FILE='PRECISIONTETA2.csv') !Открытие файла для записи значений точности при различных N
OPEN(5,FILE='PRECISIONGAMMA2.csv') !Открытие файла для записи значений точности при различных N
OPEN(6,FILE='PRECISIONTETA3.csv') !Открытие файла для записи значений точности при различных N
OPEN(7,FILE='PRECISIONGAMMA3.csv') !Открытие файла для записи значений точности при различных N
OPEN(8,FILE='PRECISIONTETA4.csv') !Открытие файла для записи значений точности при различных N
OPEN(9,FILE='PRECISIONGAMMA4.csv') !Открытие файла для записи значений точности при различных N
OPEN(10,FILE='PRECISIONGAMMA.csv') !Точность при различных N для обычного метода GAMMA
TETA1=0.;TETA2=0.;TETA3=0.;TETA4=0.
GAMMA=0.;GAMMA1=0.;GAMMA2=0.;GAMMA3=0.;GAMMA4=0.

DN=0. !DN-число точек, попадающих под график ф-ции Z=EXP(X+Y)
AA=EXP(1.) !AA-максимальное значение подынтегральной ф-ции в интервале
			!интегрирования (0,1)
DO I=1,N
    X=CSI()
	y=CSI()

    TETA1 = TETA1+EXP(X) !простейший метод вычисления
                         !среднего значения функции
    if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=TETA1/I				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I1)/I1)
		WRITE(2,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
    GAMMA1 = GAMMA1+EXP(X)*EXP(y) !простейший метод вычисления двойного интеграла
    
    if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=GAMMA1/I				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I0)/I0)
		WRITE(3,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
	Y=AA*CSI()
	if(Y<EXP(X)) then
		DN=DN+1.
	ELSE
		DN=DN+0.
	END if
	
	
    TETA2= TETA2+EXP(X)-X    !выделение главной
                             !части G(X)=1+X
    if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=0.5+TETA2/I				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I1)/I1)
		WRITE(4,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
	y=CSI()
	GAMMA2 = GAMMA2+(EXP(X+y)-(X+y))
	if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=1.+GAMMA2/I				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I0)/I0)
		WRITE(5,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
    
    Y=SQRT(1.+3.*X)-1.
    TETA3=TETA3+EXP(Y)/(1.+Y) !существенная выборка с плотностью
                              !P(X)=2*(1+X)/3
    if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=3.*TETA3/I/2.				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I1)/I1)
		WRITE(6,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
	Z=(SQRT(1.+3.*(y+X))-1.)*(SQRT(1.+3.*(y+X))-1.)
	GAMMA3=GAMMA3+(EXP(Z)/(1.+Z))*(EXP(Z)/(1.+Z))
	if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=3.*GAMMA3/I/2./2.				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I0)/I0)
		WRITE(7,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
    TETA4=TETA4+EXP(X)+EXP(1.-X) !симметризация подынтегральной функции
	if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=TETA4/I/2.				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I1)/I1)
		WRITE(8,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
	
	GAMMA4=GAMMA4+(EXP(X)+EXP(1.-X))*(EXP(y)+EXP(1.-y))
	if (mod(I,100)==0) then !чтобы не было переполнения excel файла
		IMK=GAMMA4/I/2./2.				!для вычисления точности при текущем шаге
		ksiG1 = abs((IMK-I0)/I0)
		WRITE(9,*)    I,';',ksiG1	!запись в файл csv. ; - для разделения столбцов
	end if
	
END DO

TETA1=TETA1/N		  !вычисление среднего значения функции
GAMMA1 = GAMMA1/N

TETA2=0.5+TETA2/N	  !выделение главной части
GAMMA2=1.+GAMMA2/N

TETA3=3.*TETA3/N/2.   !существенная выборка
GAMMA3=3.*GAMMA3/N/2./2.

TETA4=TETA4/N/2.      !симметризация подынтегральной функции
GAMMA4=GAMMA4/N/2./2.

GAMMA=DN/N*AA !вычисление GAMMA стандартным метдом
DN=0.
AA=exp(1.)
DO I=1,N
	y=CSI()
	Z=AA*CSI()
	if(Z<EXP(y)) then
		DN=DN+1.
	ELSE
		DN=DN+0.
	END if
	if (mod(I,100)==0) then
		IMK=GAMMA*DN*AA/I
		ksiG1 = abs((IMK-I1)/I1)
		WRITE(10,*)		I,';',ksiG1
	end if
END DO
GAMMA=GAMMA*DN*AA/N


WRITE(1,*)	'Действительное значение интеграла: ',I1
WRITE(1,*) 'Одинарный: ',TETA1,'С выделением главной части: ',TETA2    !запись результатов
WRITE(1,*) 'С существенной выборкой: ',TETA3,'С симметризацией: ',TETA4      !вычислений
WRITE(1,*)	'Действительное значение двойного интеграла', I0
WRITE(1,*) 'Двойной: ',GAMMA1,'С выделением главной части: ', GAMMA2
WRITE(1,*) 'С существенной выборкой: ', GAMMA3,'С симметризацией: ', GAMMA4
WRITE(1,*) 'Обычным методом:', GAMMA

END PROGRAM MMK2


!подпрограмма для генерации случайных чисел

REAL FUNCTION CSI()
   REAL CI
   INTEGER M       !M=1 интервал рандомных чисел (0,1)
   CALL RANDOM_SEED(M)
   CALL RANDOM_NUMBER(CI)
   CSI=CI
END FUNCTION CSI
