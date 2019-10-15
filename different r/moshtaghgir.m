%  
clc
clear all
 syms am bm cm dm D1 D2 D3 D4 D5 D6 D7 D8 C1 C2 C3 C4 C5 C6 C7 C8 p x c C A m E x y d;

 U1 = ( exp(am*x)*( C1*sin(bm*x)+C2*cos(bm*x) )+ exp(-1*am*x)*( C3*sin(bm*x)+C4*cos(bm*x))...
    + ( C5*sin(cm*x)+C6*cos(cm*x))+ C7*exp(dm*x) + C8*exp(-1*dm*x)  )* m* cos(m*y)  ;
 U2 =( exp(am*x)*(C1*bm*cos(bm*x) - C2*bm*sin(bm*x)) + exp(-am*x)*(C3*bm*cos(bm*x) - C4*bm*sin(bm*x)) +...
      am*exp(am*x)*(C2*cos(bm*x) + C1*sin(bm*x)) - am*exp(-am*x)*(C4*cos(bm*x) + C3*sin(bm*x)) +...
      C5*cm*cos(cm*x) + C7*dm*exp(dm*x) - C8*dm*exp(-dm*x) - C6*cm*sin(cm*x) )* sin(m*y) ;
 
 PHI=( exp(am*x)*( C1*sin(bm*x)+C2*cos(bm*x) )+ exp(-1*am*x)*( C3*sin(bm*x)+C4*cos(bm*x))...
 +( C5*sin(cm*x)+C6*cos(cm*x))+ C7*exp(dm*x) + C8*exp(-1*dm*x) ) * sin(m*y) ;
 %cu2,11-Au2,1111
  AA = C*diff(U2,x,2)-A*diff(U2,x,4) ;
  subs(AA,x,c) ;
  BB = diff(U2,x,1) ;
  subs(BB,x,0) ;
  CC= diff(U2,x,3);
  subs(CC,x,c) ;
 