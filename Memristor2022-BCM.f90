!***********************************************************************
!	Simple Code for Simulating The Memristor Metwork
!   Similar to "Memristor2020.f90" 
!    * Now the variable names are is as the paper,
!    * Explanations are now in English.

! The program is organized as Follows:
! 1-Module Variables (includes All Variables needed for computation)
! 2-  Program MemNetwork:
!	calls "InitialConditions" and for each task,
!	generates a random input-output-map and calls MAINPROCESS 
! 3-Subroutine InitialConditions
! 4- Subroutine MAINPROCESS: Most important one.
! it calls READroutine, WRITEroutine and COMPUTEERROR
! 5- Subroutine READroutine, calls CircuitProcess
! 6- Subroutine WRITEroutine,  calls CircuitProcess
! 7- Subroutine  COMPUTEERROR, calls READroutine
! 8- CircuitProcess
!	Writes electric equations for the nwtwork (call MatixMake)
!   Solves de equiations (call DGESV)
!	Recovers Voltage and current values (call ShowResults)
!	Updates Memristor state (call UpdateResistances)
!   Meassures the current flow (call MeassureCurr)
! 9- Subroutine MeassureCurr Meassures the current through each output 
!   node
!10- Subroutine UpdateResistances: Uses Memristor's equations and 
!   calculated voltagesto update their resistances
!11- Subroutine MatrixMake (makes the system of equations)
!12- Subroutine ShowResults Retrieves voltage and current results 
!   from the solution to equations
!***********************************************************************
! The program solves a linear equation problem using  DGESV 
!From LAPACK (linear Algebra Package). LAPACK should be installed first. 
!Recommended operations From Terminal
! $ mkdir ResultsFolder
! $ gfortran ExampleCode.f90 -o Program -llapack
! $./Program
! You should change the random number generator to the desired one.
! RECOMENDED TO OPEN IN GEANY
!***********************************************************************



!***********************************************************************
! 1- Module Variables
!***********************************************************************
module variables
!Network Parameters
integer,parameter:: N_in=3,  N_bulk=200, N_out=3 
integer,parameter:: VoltageVariables=N_in+N_bulk+N_out

!Node Voltages and total Current
real(8),dimension(N_in)::V_in
real(8),dimension(N_bulk)::V_bulk
real(8),dimension(N_out)::V_out

real(8),dimension(N_out)::TotCorK

!Edge  Values (Resistance, Conductance and Current)
real(8),dimension(N_in,N_bulk)::R_InBulk,  CondR_InBulk, I_inBulk
real(8),dimension(N_bulk,N_out)::R_BulkOut, CondR_BulkOut, I_BulkOut


! Memristor Model (beta, Rmin/max and V_threshold)

!real(8),dimension(N_in,N_bulk):: beta_InBulk !Not Used Anymore
!real(8),dimension(N_bulk,N_out):: beta_BulkOut !Not Used Anymore 
real(8),dimension(N_in,N_bulk):: Rmin_InBulk 
real(8),dimension(N_bulk,N_out):: Rmin_BulkOut 
real(8),dimension(N_in,N_bulk):: Rmax_InBulk=10000
real(8),dimension(N_bulk,N_out):: Rmax_BulkOut=10000
!real(8),dimension(N_in,N_bulk):: VT_ij !Not Used Anymore
!real(8),dimension(N_bulk,N_out):: VT_jk !Not Used Anymore
real(8),parameter:: vt0=0.915,vt1=1.3048, vthr0=4.7404,vthr1=2.4629 !Threshold values in volts
real(8),parameter:: avalue=0.1494,bvalue=1.6182


! Electrical Circuit equations and variables

real(8),dimension(VoltageVariables,VoltageVariables):: MatrixEq
real(8),dimension(VoltageVariables,1):: SolutionVector
integer:: pivot(VoltageVariables), ok


INTEGER,DIMENSION(N_in):: InputOutputMAP
real(8),dimension(N_out):: CurMed




integer Error !Hamming Distance Error
integer:: chosenInNode, MaxCurrNode, TestOutNode

integer:: i,j,k
real(8):: V0
real(8),parameter:: vwrite=-5, vread=0.0001 !Write and Read Voltages
real(8), parameter::  dt0=10.**(-3)  !Correction time in seconds
real(8), parameter::dt=dt0/10.

integer correct

integer TrainingStep, TaskNumber 
integer,parameter:: NumberOfTasks=1 !Change this number to learn several input-output maps
end module

!***********************************************************************
! 2- Program
!***********************************************************************

program MemNetwork
use variables; 
implicit none
real(8) p !will be used as random number

!Initial Conditions					   !
call InitialConditions

open (11,file="ResultsFolder/Output.dat",status="unknown")
open (12,file="ResultsFolder/Time.dat",status="unknown")
open (55,file="ResultsFolder/AllData.dat",status="unknown")
open (101,file="ResultsFolder/Resistances.dat",status="unknown")
!								DINAMICS							   !

do TaskNumber=1,NumberOfTasks
!Generate a Random input-output MAP
	do i=1,N_in
		p=rand()
		InputOutputMAP(i)=1+int(p*N_out)
		write(*,*) i, InputOutputMAP(i)
	enddo
! Solve The map
	call MAINPROCESS
	!Report The Result
	write(12,*) TrainingStep,TaskNumber
	call flush(12)
enddo

End Program

!***********************************************************************
! 3- INITIAL CONDITIONS
!***********************************************************************
subroutine InitialConditions
use variables;
implicit none
real(8) p

do i=1,N_in; do j=1,N_bulk

p=rand()
Rmin_InBulk(i,j)=1000*(0.5+p)
p=rand()
R_InBulk(i,j)=Rmin_InBulk(i,j)*(1+p)
enddo;enddo

do j=1,N_bulk; do k=1,N_out

p=rand()
Rmin_BulkOut(j,k)=1000*(0.5+p)
p=rand()
R_BulkOut(j,k)=Rmin_BulkOut(j,k)*(1+p)
enddo;enddo




call UpdateResistances

end subroutine

!***********************************************************************
! 4- MAIN PROCESS
!***********************************************************************

subroutine MAINPROCESS
use variables; 
implicit none
real(8) p

do TrainingStep=1,1000

p=rand()
chosenInNode=1+int(p*N_in) !select an input node randomly
DO correct=1,80 !Correction step within the traning step
	call READroutine(chosenInNode)
	if ((MaxCurrNode-InputOutputMAP(chosenInNode))**2>0) then
		call WRITEroutune(chosenInNode,MaxCurrNode) !Correct the network
	else
		goto 1234 !No need to correct this input node
	endif
ENDDO !correct
1234 continue

call COMPUTEERROR
if (Error==0) goto 5678 ! SUCCESS, input-output map has been read.

ENDDO !TrainingStep
call flush(11)
5678 continue
!Write Final Resistance Values
write(101,*) TrainingStep, real(R_InBulk(1,:)), "R_InBulk1"
write(101,*) TrainingStep, real(R_InBulk(2,:)), "R_InBulk2"
write(101,*) TrainingStep, real(R_InBulk(3,:)), "R_InBulk3"
write(101,*) TrainingStep, real(R_BulkOut(:,1)), "R_BulkOut1"
write(101,*) TrainingStep, real(R_BulkOut(:,2)), "R_BulkOut2"
write(101,*) TrainingStep, real(R_BulkOut(:,3)), "R_BulkOut3"

end subroutine

!***********************************************************************
! 5- READ
!***********************************************************************

subroutine READroutine(pi)
use variables 
implicit none
integer pi
chosenInNode=pi
do TestOutNode=1,N_out; 
	V0=Vread; call CircuitProcess
	CurMed(TestOutNode)=totcork(TestOutNode)
enddo
MaxCurrNode=1;      do k=1,N_out;      if (CurMed(k)>CurMed(MaxCurrNode)) MaxCurrNode=k; enddo

end subroutine

!***********************************************************************
! 6- WRITE
!***********************************************************************

subroutine WRITEroutune(pi,po)
use variables 
implicit none
integer:: pi,po
integer step
chosenInNode=pi
TestOutNode=po
V0=Vwrite
do step=1,10; call CircuitProcess; enddo
end subroutine

!***********************************************************************
! 7- Compute Error
!***********************************************************************

subroutine COMPUTEERROR
use variables 
implicit none
!End of the training step, calculate Error.
Error=0
do chosenInNode=1,N_in
    call READroutine(chosenInNode)
	if (MaxCurrNode /= InputOutputMAP(chosenInNode)) Error=Error+1

write(*,*) "Input Node:", chosenInNode,"Target Output Node:", InputOutputMAP(chosenInNode), &
& "Max Current Node:", MaxCurrNode  ,"All Currents" ,  real(curmed)
write(55,*) "Input Node:", chosenInNode,"Target Output Node:", InputOutputMAP(chosenInNode),  &
& "Max Current Node:", MaxCurrNode  ,"All Currents" ,  real(curmed)
enddo
!Report Error
write(11,*) "Training Step:", TrainingStep,"Hamming Error:", Error, "Map#", TaskNumber, "correction#",correct
write(*,*)  "Training Step:", TrainingStep,"Hamming Error:", Error, "Map#", TaskNumber, "correction#",correct

end subroutine


!***********************************************************************
! 8- Circuit Process
!***********************************************************************

subroutine CircuitProcess
use variables!Mi sistema
implicit none
	call MatixMake !Make the equation matrix
	call DGESV(VoltageVariables, 1, MatrixEq, VoltageVariables, pivot, SolutionVector, VoltageVariables, ok)
	call ShowResults
	call UpdateResistances
    call MeassureCurr
end subroutine

!***********************************************************************
! 9- Meassure Current
!***********************************************************************
subroutine MeassureCurr
use variables
implicit none
	TotCorK=0
	do k=1,N_out;	do j=1,N_bulk;	TotCorK(k)=TotCorK(k)+I_BulkOut(j,k); 	enddo;
	enddo; 	
end subroutine	

!***********************************************************************
! 10- Update Resistances
!***********************************************************************

subroutine UpdateResistances
use variables
implicit none
real(8) Vm !memristor Voltage
real(8):: F,f_b !Memrsistor F function
real(8):: current, omega, constant
real(8),parameter:: mu_Over_D2=1
!**********

do i=1,N_in; do j=1,N_bulk

Vm=V_in(i)-V_bulk(j)  !Aqui tiene convención  de signo opuesto a mi programa 
current=Vm/R_InBulk(i,j)
omega=(Rmax_InBulk(i,j)-R_InBulk(i,j))/(Rmax_InBulk(i,j)-Rmin_InBulk(i,j))
Constant= (Rmax_InBulk(i,j)-Rmin_InBulk(i,j))*Rmin_InBulk(i,j)*mu_Over_D2 !mu/D**2=1
 

f_b=0
IF ((omega>0.0d0) .and. (omega<1.d0)) THEN
f_b=bvalue
if ((Vm<vt0) .and. (Vm>-vt1)) then
f_b=avalue
endif
ELSE
if ((Vm>vthr0) .and. (omega .le. 0.0d0)) then !R=Roff
f_b=bvalue
endif

if ((Vm<-vt1) .and. (omega .ge. 1.d0)) then !R=Ron
f_b=bvalue
endif
ENDIF
!f= beta_ij(i,j)*(Vm-0.5*(abs(Vm+VT_ij(i,j))-abs(Vm-VT_ij(i,j)) ))

R_InBulk(i,j)=R_InBulk(i,j)-f_b*current*Constant*dt



if (R_InBulk(i,j)> Rmax_InBulk(i,j)) R_InBulk(i,j)= Rmax_InBulk(i,j)
if (R_InBulk(i,j)< Rmin_InBulk(i,j)) R_InBulk(i,j)= Rmin_InBulk(i,j)
enddo; enddo



do j=1,N_bulk; do k=1,N_out
Vm=V_bulk(j)-V_out(k)
current=Vm/R_BulkOut(j,k)
omega=(Rmax_BulkOut(j,k)-R_BulkOut(j,k))/(Rmax_BulkOut(j,k)-Rmin_BulkOut(j,k))
Constant= (Rmax_BulkOut(j,k)-Rmin_BulkOut(j,k))*Rmin_BulkOut(j,k)*mu_Over_D2 !mu/D**2=1



f_b=0
IF ((omega>0.0d0) .and. (omega<1.0d0)) THEN
f_b=bvalue
if ((Vm<vt0) .and. (Vm>-vt1)) then
f_b=avalue
endif
ELSE
if ((Vm>vthr0) .and. (omega .le. 0.0d0)) then !R=Roff
f_b=bvalue
endif

if ((Vm<-vt1) .and. (omega .ge. 1.0d0)) then !R=Ron
f_b=bvalue
endif
ENDIF
!f= beta_ij(i,j)*(Vm-0.5*(abs(Vm+VT_ij(i,j))-abs(Vm-VT_ij(i,j)) ))

R_BulkOut(j,k)=R_BulkOut(j,k)-f_b*current*Constant*dt


if (R_BulkOut(j,k)> Rmax_BulkOut(j,k)) R_BulkOut(j,k)= Rmax_BulkOut(j,k)
if (R_BulkOut(j,k)< Rmin_BulkOut(j,k)) R_BulkOut(j,k)= Rmin_BulkOut(j,k)
enddo; enddo
CondR_InBulk=1.d0/R_InBulk
CondR_BulkOut=1.d0/R_BulkOut

end subroutine

!***********************************************************************
! 11- Make The matrix of coefficients for the equations
!***********************************************************************

subroutine MatixMake
use variables
implicit none
integer EqNumber !Equation for which I generate coefficients
! The equation system has VoltageVariables=N_in+N_bulk+N_out variables,
!which are the voltages.
! 2 equations are given by the voltages in chosenInNode and TestOutNode
!the other VoltageVariables-2 equations are the no net current on all other
!N_in+N_bulk+N_out-2 nodes. The current on a network edge is expressed 
!as the voltage difference multiplied by the conductivity of the
!corresponding memristor.



MatrixEq=0
SolutionVector=0

!***********************************************************************
!First N_in equations: !No net current on input nodes.
!sum over j of  I_InBulk(i,j) is 0 , except for i=chosenInNode
do i=1,N_in; do j=1,N_bulk
EqNumber=i
MatrixEq(EqNumber,i)=MatrixEq(EqNumber,i)+CondR_InBulk(i,j)
MatrixEq(EqNumber,j+N_in)=-CondR_InBulk(i,j)
enddo; enddo


i=chosenInNode
MatrixEq(i,:)=0
MatrixEq(i,i)=1
SolutionVector(i,1)=V0 

!***********************************************************************
!Next N_bulk equations (N_in+1 to N_in+N_bulk)
 !Conserved charge in Bulk node j 
 !Sum over i of I_inBulk(i,j) plus sum sover  k of I_BulkOut(X,j) equals zero.
! la variable I_inBulk es la VoltageVariables+j+(i-1)*N_bulk
do j=1,N_bulk
EqNumber=j+N_in

do i=1,N_in
MatrixEq(EqNumber,i)=CondR_InBulk(i,j)
MatrixEq(EqNumber,j+N_in)=MatrixEq(EqNumber,j+N_in)-CondR_InBulk(i,j)
enddo

do k=1,N_out
MatrixEq(EqNumber,j+N_in)=MatrixEq(EqNumber,j+N_in)-CondR_BulkOut(j,k)
MatrixEq(EqNumber,k+N_in+N_bulk)=CondR_BulkOut(j,k)
enddo

enddo
!***********************************************************************
!Last N_out Equations: 
do k=1,N_out
EqNumber=N_in+N_bulk+k
do j=1,N_bulk
MatrixEq(EqNumber,j+N_in)=CondR_BulkOut(j,k)
MatrixEq(EqNumber,K+N_in+N_bulk)=MatrixEq(EqNumber,K+N_in+N_bulk)-CondR_BulkOut(j,k)
enddo
enddo

k=TestOutNode
EqNumber=N_in+N_bulk+k
MatrixEq(EqNumber,:)=0
MatrixEq(EqNumber,N_in+N_bulk+k)=1


end subroutine

!***********************************************************************
! 12- Retrieve voltage and current results from the solution to equations
!***********************************************************************
subroutine ShowResults
use variables
implicit none

do i=1,N_in; V_in(i)=SolutionVector(i,1); enddo
do j=1,N_bulk; V_bulk(j)=SolutionVector(N_in+j,1); enddo
do k=1,N_out; V_out(k)=SolutionVector(N_in+N_bulk+k,1); enddo

do i=1,N_in; do j=1,N_bulk
I_inBulk(i,j)= (V_in(i)-V_bulk(j))*condR_InBulk(i,j) 
enddo; enddo

do j=1,N_bulk; do k=1,N_out
I_BulkOut(j,k)= (V_bulk(j)-V_out(k))*condR_BulkOut(j,k) 
enddo; enddo

end subroutine
