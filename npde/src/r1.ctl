$PROB Example 1
; Create a model to fit the data that you find more suitable according the the EDA. 
$INPUT
        ; In this block we tell NM to label our dataset accordingly

ID      ; Individual identifier number 
TIME    ; Time column (This can be clock or nominal time) 
AMT     ; 
DROP      ; Dependent variable column 
DV       ; Dosage 
MDV     ; Missing DVs 


; After we did add the labels that matchs our dataset in the index we can move to set the 
; directory

$DATA  data/dat1.csv  ignore=@        ; This refers to ignoring all non-numerical values 
                        ; You can replace this with C if thats your column 
                        ; or you can do #. 
; Additional values can be added that is specific to dataset so NMTRAN can know what to do. 
; e.g. ignore(CMT == 1) - I think its self explantory. 

$SUBR   ; Here we add subroutines that are defined and written as functions in
        ; PREDPP. you can check the values for them in nonmem guide or in your terminal
        ;       "nmhelp ADVAN" 

ADVAN1  ; subroutine PK library to implement one-compartment linear model 
TRANS2  ; Translator that reads the PK parameters as CL and V from the PK block


$PK
        TVCL    = THETA(1)      ; Population typical value for CL
        TVV     = THETA(2)      ; Population typical value for V

        ; Individual model 

        CL      =       TVCL * EXP(ETA(1))
        V       =       TVV  * EXP(ETA(2))
; For scaling parameter if we want the prediction and observation to have the same unit since F (prediction) will be devided by this value. 
        S1      =       V  

$ERROR

; In the error block, here we define the statistical model for the observation
; e.g. y = f(x) + error 
IPRED   =       F
Y       =       IPRED * (1 + ERR(1))

; We just assigned the value here since we are interested in the ipred for diagnostic plots. 

; We define the initial values
; try to have best guess. 

$THETA NAMES(CL,V) 
(0, 5) 
(0. 40) 

$OMEGA BLOCK(2) NAMES(CL,V) values(0.1,0.01) 
$SIGMA 0.01 

; In the estimation block; we define the algorithm that will approximate the value of maximum likelihood and estimate the PK parameters. Different options out there, and for this course, we are going to stick with first-order one. 

$EST    MET=1 INTER NOABORT PRINT=1 
; Table block request a dataframe (Tab separated - or actually you can request whatever you want) for post-processing in whatever other languages for data wrangling!! 
$TAB   ID TIME DV AMT MDV IPRED PRED CWRES 
        NPDE NPD ESAMPLE=1000 ONEHEADER NOPRINT NOAPPEND 
        FILE=r1.tab 



