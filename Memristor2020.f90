
 
module vars
integer,parameter:: Ni=3
integer,parameter:: Nj=120
integer,parameter:: Nk=Ni
real(8),dimension(Ni)::Vi
real(8),dimension(Nj)::Vj
real(8),dimension(Nk)::Vk
real(8),dimension(Ni,Nj)::Rij
real(8),dimension(Nj,Nk)::Rjk

real(8),dimension(Ni,Nj)::CondRij
real(8),dimension(Nj,Nk)::CondRjk


real(8),dimension(Ni,Nj)::Iij
real(8),dimension(Nj,Nk)::Ijk
real(8),dimension(Nk)::TotCorK
integer,parameter:: VariablesVoltaje=Ni+Nj+Nk
!integer,parameter::Variables=VariablesVoltaje+VariablesCorrienteA+VariablesCorrienteB
!real(8),dimension(Variables):: VectorIncognitas
real(8),dimension(VariablesVoltaje,variablesVoltaje):: MatrizEq
real(8),dimension(VariablesVoltaje,1):: VectorSoluciones
INTEGER,DIMENSION(Ni):: TABLAElegidos
real(8),dimension(Nk):: CurMed


integer:: PuntaI=2
integer:: maxcorr,maxcorr2
integer:: PuntaO
integer nroeq
integer:: i,j,k
real(8):: V0=1.00
real(8):: Vescala=1.0
integer:: pivot(VariablesVoltaje), ok

real(8), parameter::  dt=1 !dt=0.0001

INTEGER BIGstep

!Modelo De Memristor
!*******************

real(8),dimension(Ni,Nj):: beta_ij 
real(8),dimension(Nj,Nk):: beta_jk 
real(8),dimension(Ni,Nj):: RON_ij 
real(8),dimension(Nj,Nk):: RON_jk 
real(8),dimension(Ni,Nj):: ROFF_ij=10.**6
real(8),dimension(Nj,Nk):: ROFF_jk=10.**6
real(8),dimension(Ni,Nj):: VT_ij
real(8),dimension(Nj,Nk):: VT_jk



integer experimento
integer Nroexperimento
end module

program C1
use vars; 
implicit none
real(8) p
integer semilla



!***********************************************************************
!								Cond. Iniciales						   !
!***********************************************************************



call condicionesiniciales


open (11,file="ResultadosF/salidaXXX.dat",status="unknown")
open (12,file="ResultadosF/TiemposXXX.dat",status="unknown")
open (55,file="ResultadosF/AllXXX.dat",status="unknown")
open (101,file="ResultadosF/ResXXX.dat",status="unknown")

!***********************************************************************
!								DINAMICA							   !
!***********************************************************************


do NroExperimento=1,1

do i=1,ni;
p=rand()
TABLAElegidos(i)=1+int(p*Nk)
write(*,*) i, TABLAElegidos(i)
enddo
write(*,*) "*********************"


call BIGPROCESS
write(12,*) Experimento,NroExperimento
call flush(12)
enddo



End Program



subroutine BIGPROCESS
use vars; 
implicit none
real(8) p
integer nerror
integer nerror2



do experimento=1,5000

p=rand()
PuntaI=1+int(p*ni)
DO Bigstep=1,80

call leo(puntaI)

if ((maxcorr-tablaelegidos(PuntaI))**2>0) then
call corrijo(maxcorr)
!elseif (maxcorr<1.05*maxcorr2) then
!call corrijo(maxcorr2)
else

goto 1234
endif

ENDDO !Bigstep

 
1234 continue


Nerror=0
Nerror2=0
do PuntaI=1,ni
        call leo(PuntaI)
if (maxcorr /= tablaelegidos(PuntaI)) Nerror=Nerror+10
Nerror2=Nerror2+(maxcorr -tablaelegidos(PuntaI))**2
!if (maxcorr<1.01*maxcorr2) Nerror=Nerror+1
write(*,*) "Boton:", PuntaI,"Target:", tablaelegidos(PuntaI), "Maxcorr=", maxcorr  ,"<>" ,  real(curmed)
write(55,*) "Boton:", PuntaI,"Target:", tablaelegidos(PuntaI), "Maxcorr=", maxcorr  ,"<>" ,  real(curmed)
!read(*,*)
enddo
write(11,*) Experimento, Nerror,Nerror2, "Semilla XXX ", "Nro", NroExperimento, "BS",Bigstep
write(*,*) Experimento, Nerror,Nerror2, "Semilla XXX ", "Nro", NroExperimento, "BS",Bigstep
!write(101,*) Experimento, real(Rij(1,:)), "Rij1"
!write(101,*) Experimento, real(Rij(2,:)), "Rij2"
!write(101,*) Experimento, real(Rij(3,:)), "Rij3"
!write(101,*) Experimento, real(Rjk(:,1)), "Rjk1"
!write(101,*) Experimento, real(Rjk(:,2)), "Rjk2"
!write(101,*) Experimento, real(Rjk(:,3)), "Rjk3"



if (nerror==0) goto 5678


ENDDO !Experiment
call flush(11)
5678 continue

write(101,*) Experimento, real(Rij(1,:)), "Rij1"
write(101,*) Experimento, real(Rij(2,:)), "Rij2"
write(101,*) Experimento, real(Rij(3,:)), "Rij3"
write(101,*) Experimento, real(Rjk(:,1)), "Rjk1"
write(101,*) Experimento, real(Rjk(:,2)), "Rjk2"
write(101,*) Experimento, real(Rjk(:,3)), "Rjk3"

end subroutine


subroutine leo(pi)
use vars 
implicit none
integer pi
PuntaI=pi
do PuntaO=1,nk; 
	V0=0.001; call Proceso
	CurMed(PuntaO)=totcork(PuntaO)
enddo
maxcorr=1;      do k=1,nk;      if (CurMed(k)>CurMed(maxcorr)) maxcorr=k; enddo

maxcorr2=1; if (maxcorr==1) maxcorr2=2; do k=1,nk; 	if (maxcorr==K) CYCLE
if (CurMed(k)>CurMed(maxcorr2)) maxcorr2=k; enddo


end subroutine



subroutine corrijo(po)
use vars 
implicit none
integer po
integer step
PuntaO=po
V0=0.2
do step=1,5; call Proceso; enddo
end subroutine



subroutine condicionesiniciales
use vars;
implicit none
real(8) p


do i=1,ni; do j=1,nj
p=rand()
beta_ij(i,j)=1.0-p*0.2 
p=rand()
p=rand()
VT_ij(i,j)= 10.**(-1.)*(1+p)/2.d0
p=rand()
Ron_ij(i,j)=100*(1+p)/2.d0
enddo;enddo

do j=1,nj; do k=1,nk
p=rand()
beta_jk(j,k)=1.0-p*0.2 
p=rand()
p=rand()
VT_jk(j,k)= 10.**(-1.)*(1+p)/2.d0
p=rand()
Ron_jk(j,k)=100*(1+p)/2.d0
enddo;enddo

call actualizaResVolt

end subroutine





subroutine Proceso
use vars!Mi sistema
implicit none
	call ArmoMatriz
	call DGESV(variablesVoltaje, 1, MatrizEq, variablesVoltaje, pivot, VectorSoluciones, variablesVoltaje, ok)
	call MuestroResults
	call actualizaResvolt
    call Midecurr
end subroutine



subroutine Midecurr
use vars!Mi sistema
implicit none
	TotCorK=0
	do k=1,nk;	do j=1,nj;	TotCorK(k)=TotCorK(k)+Ijk(j,k); 	enddo;!write(*,*) "Corriente en el nodo de salida numero", k, TotCorK(k)
	enddo; 

	
end subroutine	





subroutine actualizaResVolt
use vars!Mi sistema
implicit none
real(8) vm
real(8) f
do i=1,ni; do j=1,nj

Vm=Vi(i)-Vj(j)
f= beta_ij(i,j)*(Vm-0.5*(abs(Vm+VT_ij(i,j))-abs(Vm-VT_ij(i,j)) ))
Rij(i,j)=Rij(i,j)+f*dt
if (Rij(i,j)> ROFF_ij(i,j)) Rij(i,j)= ROFF_ij(i,j)
if (Rij(i,j)< RON_ij(i,j)) Rij(i,j)= RON_ij(i,j)
enddo; enddo

do j=1,nj; do k=1,nk
Vm=Vj(j)-Vk(k)
f= beta_jk(j,k)*(Vm-0.5*(abs(Vm+VT_jk(j,k))-abs(Vm-VT_jk(j,k)) ))
Rjk(j,k)=Rjk(j,k)+f*dt
if (Rjk(j,k)> ROFF_jk(j,k)) Rjk(j,k)= ROFF_jk(j,k)
if (Rjk(j,k)< RON_jk(j,k)) Rjk(j,k)= RON_jk(j,k)
enddo; enddo
CondRij=1.d0/Rij
CondRjk=1.d0/Rjk

end subroutine






subroutine ArmoMatriz
use vars!Mi sistema
implicit none


MatrizEq=0
vectorsoluciones=0

!***********************************************************************
!Primeras i ecuaciones: !suma sobre j de I1(i,j) es 0 , salvo para i=PuntaI
do i=1,Ni; do j=1,Nj
NroEq=i
MatrizEq(NroEq,i)=MatrizEq(NroEq,i)+CondRij(i,j)
MatrizEq(NroEq,j+Ni)=-CondRij(i,j)
enddo; enddo


i=PuntaI
MatrizEq(i,:)=0
Matrizeq(i,i)=1
VectorSoluciones(i,1)=V0

!***********************************************************************
!siguentes j ecuaciones (Ni+1 a Ni+nj)
 !Suma sobre i de Iij mas la suma sobre k de Ijk (un 1 en cada elemenro
! la variable Iij es la VariablesVoltaje+j+(i-1)*nj
do j=1,nj
NroEq=j+ni

do i=1,ni
Matrizeq(NroEq,i)=Condrij(i,j)
Matrizeq(NroEq,j+ni)=Matrizeq(NroEq,j+ni)-CondRij(i,j)
enddo

do k=1,nk
Matrizeq(NroEq,j+ni)=Matrizeq(NroEq,j+ni)-CondRjk(j,k)
Matrizeq(NroEq,k+ni+nj)=CondRjk(j,k)
enddo

enddo
!***********************************************************************

do k=1,nk
NroEq=ni+nj+k
do j=1,nj
Matrizeq(Nroeq,j+ni)=CondRjk(j,k)
Matrizeq(Nroeq,K+ni+NJ)=Matrizeq(Nroeq,K+ni+NJ)-CondRjk(j,k)
enddo
enddo

k=PuntaO
NroEq=ni+nj+k
MatrizEq(NroEq,:)=0
MatrizEq(NroEq,Ni+Nj+k)=1


end subroutine


subroutine muestroresults
use vars!Mi sistema
implicit none

do i=1,ni; Vi(i)=VectorSoluciones(i,1); enddo
do j=1,nj; Vj(j)=VectorSoluciones(ni+j,1); enddo
do k=1,nk; Vk(k)=VectorSoluciones(ni+nj+k,1); enddo

do i=1,ni; do j=1,nj
Iij(i,j)= (Vi(i)-Vj(j))*condRij(i,j) 
enddo; enddo

do j=1,nj; do k=1,nk
Ijk(j,k)= (Vj(j)-Vk(k))*condRjk(j,k) 
enddo; enddo

end subroutine
