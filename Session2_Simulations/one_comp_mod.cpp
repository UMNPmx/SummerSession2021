[set] delta = 0.01, end = 50 //set options for simulations

[PROB]
// $PROB
- author: Shen Cheng
- date: `r Sys.Date()`
- comments: this is a demonstration
- one compartment model + first order absportion
- random effects added on `CL` and `VC` as exponential distributions
- proportional error model
- model covariate: SEX ~ CL

$PARAM
TVCL = 1, TVVC = 20, TVKA = 1, SEXCL=0.5, SEX=0

$CMT  @annotated 
// similar to $MODEL in NONMEM
// if @annotated is not used, simply declare compartment names with a space

DEPOT  : Depot compartment (mass)
CENT   : Central compartment (mass)

$MAIN
// similar to $PK block in NONMEM

double CL   = TVCL*exp(ETA(1))*pow(SEXCL, SEX); //power function in C++ witten as pow(,)
double VC   = TVVC*exp(ETA(2))                ; // "double" means defining the data structure as double
double KA   = TVKA                            ; // ";" is required at the end of each defined variable


$GLOBAL 
  // declare variables you want to either use further in computation
  // OR derive variables you want calculated from computations
  
  // bool cure = FALSE; // declaring variables needs ";" at the end
  // double CL, VC, KA; // declaring variables so no "double" is need in $MAIN
#define CP (CENT/VC)  // # define does not require ";"
  
  //[PREMABLE]
  // similar to $GLOBAL to some extent
  // However, unlike $MAIN $TABLE and $ODE, $PREMABLE only call once during the simulation
  
  $OMEGA @block
  // @block generates a block matrix with the formulation (1,1), (1,2), (2,2)....
  // @correlation will specify values to be corrleations instead as variances
  
  0.1 0.02 0.3
  
  $SIGMA 
    // sigma variances go here
    
    0.01
  
  $ODE
    //similar to $DES in NONMEM
    //differential equations go here
    
    dxdt_DEPOT = -KA*DEPOT;
  dxdt_CENT  = KA*DEPOT - CL*CP;
  
  $TABLE
    //another block could be used for defining secondary parameters
    
    double DV = CP * (1 + EPS(1));
  
  $CAPTURE
    //This is a block to identify variables that should be captured in the simulated output
    
    CP, DV