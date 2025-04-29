// Basic RBC Model with: 
//  1) EZ preferences
//  2) Adjustment costs on capital
//  3) Growth, Unit root 
//  4) Taylor rule for inflation
// 
// Jesus Fernandez-Villaverde
//
// Bala Cynwyd, August 16, 2010

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables
//----------------------------------------------------------------

var 

// Utility variables
u v ev m

// Allocation variables 
y c l k i  

// Input prices
r w 

// Bond yield variables
R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 R17 R18 R19 R20 yield 

// Inflation
infl 

// Shocks
z;

//----------------------------------------------------------------
// 2. Exogenous variables
//----------------------------------------------------------------

varexo e pphi;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Utility function
cte nnu bbeta ggamma eis ttheta

// Technology 
aalpha ddelta llambda ssigma ttau a1 a2 

// Taylor rule
ggamma_R ggamma_inf ggamma_y targetinf ssigmam y_ss R_ss;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences 
nnu     = 0.5201;
bbeta   = 0.994;
ggamma  = 79.34;
eis     = 1.731;
ttheta  = (1-ggamma)/(1-1/eis);

// Technology
aalpha  = 0.3;
ddelta  = 0.0294;
llambda = 0.0045;
ssigma  = 0.008;
ttau    = 0.4;//0.032;
a1      = (exp(llambda)-1+ddelta)/(1-ttau);
a2      = (exp(llambda)-1+ddelta)^(1/ttau);
  
// Taylor Rule
ggamma_R   = 0.8;
ggamma_inf = 1.5;
ggamma_y   = 0.25;
targetinf  = 1.005;
ssigmam    = 0.005;

//----------------------------------------------------------------
// 5. Steady State
//----------------------------------------------------------------

psi   = exp(llambda)-1+ddelta;
omega = ((1/(aalpha*(exp(llambda)^(1-aalpha))))*...
        (1/(bbeta*(exp(llambda)^(nnu*(1-ggamma)/ttheta-1)))-1+ddelta))^(1/(aalpha-1));
phi   = (nnu/(1-nnu))*(1-aalpha)*(exp(llambda)^(1-aalpha))*(omega^aalpha);
l_ss  = phi/((exp(llambda)^(1-aalpha))*(omega^aalpha)-omega*psi+phi);
k_ss  = omega*l_ss;
i_ss  = psi*k_ss;
y_ss  = (k_ss^aalpha)*((exp(llambda)*l_ss)^(1-aalpha));
r_ss  = aalpha*y_ss/k_ss;
w_ss  = (1-aalpha)*y_ss/l_ss;
c_ss  = y_ss-i_ss;
m_ss  = bbeta*(exp(llambda)^(nnu*(1-ggamma)/ttheta-1));
R_ss  = targetinf/m_ss;
norm  = (1/(1-bbeta*(exp(llambda)^(nnu*(1-ggamma)/ttheta)))^(ttheta/(1-ggamma)))*(c_ss^nnu)*((1-l_ss)^(1-nnu));
cte   = 1/norm; % normalizes v_ss to 1
u_ss  = cte*(c_ss^nnu)*((1-l_ss)^(1-nnu));
v_ss  = (1/(1-bbeta*(exp(llambda)^(nnu*(1-ggamma)/ttheta)))^(ttheta/(1-ggamma)))*u_ss;
ev_ss = v_ss^(1-ggamma);

//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model; 
  
  // 1. Utility
  u = cte*(c^nnu)*((1-l)^(1-nnu));
  
  // 2. Value function
  v = (u^((1-ggamma)/ttheta)
      +bbeta*(z^(nnu*(1-ggamma)/ttheta))*(ev^(1/ttheta)))^(ttheta/(1-ggamma));
  
  // 3. Expected value function
  ev = v(+1)^(1-ggamma);

  // 4. Static leisure-consumption
  ((1-nnu)/nnu)*c/(1-l) = w;
  
  // 5. Pricing Kernel
  m = bbeta*((z(-1))^(nnu*(1-ggamma)/ttheta-1))
           *((u/u(-1))^((1-ggamma)/ttheta))*(c(-1)/c)
           *((v^(1-ggamma)/ev(-1))^(1-1/ttheta));

  // 6. Euler equation for capital
  (i/k(-1))^(1/ttau) = m(+1)*
  (a2*(r(+1))+((i(+1)/k)^(1/ttau))*(1-ddelta+a1+(a2/(ttau-1))*((i(+1)/k)^(1-(1/ttau))))); 
   
  // 7. Euler equations for nominal bonds
  1/R1     = m(+1)/(1+infl(+1));
  1/(R2^2) = (m(+1)*(1/R1(+1)))/(1+infl(+1)); 
  1/(R3^3) = (m(+1)*(1/(R2(+1)^2)))/(1+infl(+1)); 
  1/(R4^4) = (m(+1)*(1/(R3(+1)^3)))/(1+infl(+1)); 
  1/(R5^5) = (m(+1)*(1/(R4(+1)^4)))/(1+infl(+1)); 
  1/(R6^6) = (m(+1)*(1/(R5(+1)^5)))/(1+infl(+1)); 
  1/(R7^7) = (m(+1)*(1/(R6(+1)^6)))/(1+infl(+1)); 
  1/(R8^8) = (m(+1)*(1/(R7(+1)^7)))/(1+infl(+1));
  1/(R9^9) = (m(+1)*(1/(R8(+1)^8)))/(1+infl(+1)); 
  1/(R10^10) = (m(+1)*(1/(R9(+1)^9)))/(1+infl(+1)); 
  1/(R11^11) = (m(+1)*(1/(R10(+1)^10)))/(1+infl(+1)); 
  1/(R12^12) = (m(+1)*(1/(R11(+1)^11)))/(1+infl(+1)); 
  1/(R13^13) = (m(+1)*(1/(R12(+1)^12)))/(1+infl(+1)); 
  1/(R14^14) = (m(+1)*(1/(R13(+1)^13)))/(1+infl(+1)); 
  1/(R15^15) = (m(+1)*(1/(R14(+1)^14)))/(1+infl(+1)); 
  1/(R16^16) = (m(+1)*(1/(R15(+1)^15)))/(1+infl(+1)); 
  1/(R17^17) = (m(+1)*(1/(R16(+1)^16)))/(1+infl(+1)); 
  1/(R18^18) = (m(+1)*(1/(R17(+1)^17)))/(1+infl(+1)); 
  1/(R19^19) = (m(+1)*(1/(R18(+1)^18)))/(1+infl(+1)); 
  1/(R20^20) = (m(+1)*(1/(R19(+1)^19)))/(1+infl(+1)); 
  
  // 8. Yield Curve difference between 20 quarters and 1 quarter
  yield = R20-R4;

  // 9. Production function
  y = (k(-1)^aalpha)*((z*l)^(1-aalpha));
  
  // 10. Rental rate of capital
  r = aalpha*y/(k(-1));

  // 11. Wage
  w = (1-aalpha)*y/l;
  
  // 12. Resource constraint
  c+i = y; 

  // 13. Law of motion for capital
  z*k = (1-ddelta)*(k(-1))+(a1+(a2/(1-1/ttau))*((i/k(-1))^(1-(1/ttau))))*(k(-1));  

  // 14. Law of motion for productivity
  z = exp(llambda+ssigma*e);
  
  // 15. Taylor Rule
  R1/R_ss = ((R1(-1)/R_ss)^ggamma_R)*
            ((((1+infl)/targetinf)^ggamma_inf)*((y/y_ss)^ggamma_y))^(1-ggamma_R)
            *exp(ssigmam*pphi);

end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------

initval;
  l = l_ss;
  k = k_ss;
  y = y_ss;
  c = c_ss;
  i = i_ss;
  r = r_ss;
  w = w_ss;
  m = m_ss;
  u = u_ss;
  v = v_ss;
  ev = ev_ss;
  z = exp(llambda);
  e = 0;
  R1 = R_ss;
  R2 = R_ss;
  R3 = R_ss;
  R4 = R_ss;
  R5 = R_ss;
  R6 = R_ss;
  R7 = R_ss;
  R8 = R_ss;
  R9 = R_ss;
  R10 = R_ss;  
  R11 = R_ss;
  R12 = R_ss;
  R13 = R_ss;
  R14 = R_ss;
  R15 = R_ss;
  R16 = R_ss;
  R17 = R_ss;
  R18 = R_ss;
  R19 = R_ss;
  R20 = R_ss;
  yield = 0;
  infl = targetinf-1;
  pphi = 0;
end;
    
shocks;
  var e = 1;
  var pphi = 1;
end;

steady;  

check;

stoch_simul(hp_filter = 1600, irf = 20, order = 2) infl R1 R4 R20 yield k m y c;