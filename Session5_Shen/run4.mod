;; 1. Based on:
;; 2. Description: PD model PK PARAM
;; x1. Author: user
;; 3. Label:

$PROBLEM PK

$INPUT C ID TIME NTIME DV AMT CMT EVID MDV TYPE
       GENO SEX ETHNICITY AGE WT HT CRCL DOSE BMI KTRI CLI VCI

$DATA r4pmx_dataset_simulation1_pkparam.csv IGNORE=C  IGNORE(CMT.EQ.5)

$SUBROUTINES ADVAN13 TOL=12
$MODEL
  COMP =(DEPOT)
  COMP =(TRANS1)
  COMP =(TRANS2)
  COMP =(TRANS3)
  COMP =(CENT)
  COMP =(RESP)

$PK

KTR    = KTRI
CL     = CLI
VC     = VCI  

TVIC50 = THETA(1)
TVKIN  = THETA(2)
TVKOUT = THETA(3)

IC50 = TVIC50 * EXP(ETA(1))
KIN  = TVKIN  * EXP(ETA(2))
KOUT = TVKOUT * EXP(ETA(3))


S6=1

A_0(6) = KIN/KOUT


$DES
CP = A(5)/VC

DADT(1)= - KTR * A(1)
DADT(2)=   KTR * (A(1) - A(2))
DADT(3)=   KTR * (A(2) - A(3))
DADT(4)=   KTR * (A(3) - A(4))

DADT(5)=   KTR * A(4) - CL*CP

INH = CP/(IC50 + CP)

DADT(6)=   KIN*(1-INH) - KOUT*A(6)


$ERROR
RESP   = A(6)
IPRED = RESP
Y   = IPRED + IPRED*ERR(1)

$THETA
(0, 1)  ; IC50
(0, 50)  ; KIN
(0, 5) ; KOUT

$OMEGA
 0.3  ; IIV KOUT
 0.2  ; IIV KIN
 0 FIX; IIV KOUT


$SIGMA
 0.04 ; Prop ERR

$EST METHOD=1 INTER PRINT=5 MAX=9999 NSIG=3 SIGL=9 POSTHOC
$COV SIGL=12 MATRIX=R PRINT=E

$TABLE ID TIME DV MDV EVID IPRED PRED CWRES
KTR CL VC IC50 KIN KOUT ETA1 ETA2 ETA3
GENO SEX ETHNICITY AGE WT HT CRCL DOSE BMI
ONEHEADER NOPRINT FILE=sdtab004