clc
clear all

%%%%%%%%%%%numeric and symbolic variables

syms  C1 C2 C3 C4 C5 C6 x y;
C=150;%%% C/MIU
A=15; %%%A/MIU

%%%dimensions
d=1;
c=1;

%%%%coefficients of the general equation

PHI=zeros(21,11) ;
U1=0 ;
U2=0 ;

%%%for loop to calculate coefficients corresponding to fourier series for each iteration 
for n=1:2:61           

    %coefficients found from the original pde (all are functions of n and m ) 
    m=(n*pi)/(2*d);
    P= ((C^3)/(27*(A^3))) ;
    Q= ((C)/(6*(A^2))) ;
    M= (P+((m^2)/(2*A))-Q)^2 ;
    N= (-(C^2)/(9*(A^2)))+(1/(3*A)) ;
    O = ((P+((m^2)/(2*A))+sqrt(M+(N^3)))-Q)^(1/3) ;
    R= (P+((m^2)/(2*A))+sqrt((N^3)+M)-Q)^(1/3) ;
    S= 0.5*(((N/O)+R)^(1/3)) ;
    F= (N/(2*O))+(C/(3*A))-(R/2) ;
    am= sqrt((C/(3*A))-(N/O)+R);
    cm= S/sqrt(2*(F+sqrt((S^2)+(F^2)))) ;
    bm= sqrt((F+sqrt((S^2)+(F^2)))/2) ;
    
    %%%it is notable that we have totally 6 coefficients to derive 
                  
    %%%u1,111(c,y)=0    "rate of changes in curvature is zero at the boundaries"            (1) 
 
    eq1=3*bm^2*exp(bm*c)*(C3*cm*cos(c*cm) - C4*cm*sin(c*cm)) - exp(-bm*c)*(C5*cm^3*cos(c*cm) - C6*cm^3*sin(c*cm)) - exp(bm*c)*(C3*cm^3*cos(c*cm) - C4*cm^3*sin(c*cm)) + 3*bm^2*exp(-bm*c)*(C5*cm*cos(c*cm) - C6*cm*sin(c*cm)) + bm^3*exp(bm*c)*(C4*cos(c*cm) + C3*sin(c*cm)) - bm^3*exp(-bm*c)*(C6*cos(c*cm) + C5*sin(c*cm)) - 3*bm*exp(bm*c)*(C4*cm^2*cos(c*cm) + C3*cm^2*sin(c*cm)) + 3*bm*exp(-bm*c)*(C6*cm^2*cos(c*cm) + C5*cm^2*sin(c*cm)) + C1*am^3*exp(am*c) - C2*am^3*exp(-am*c);

    %%%u1,111(-c,y)=-u1,111(c,y)  "rate of changes in curvature is symmetric"               (2)

    eq2=3*bm^2*exp(-bm*c)*(C3*cm*cos(c*cm) + C4*cm*sin(c*cm)) - exp(bm*c)*(C5*cm^3*cos(c*cm) + C6*cm^3*sin(c*cm)) - exp(-bm*c)*(C3*cm^3*cos(c*cm) + C4*cm^3*sin(c*cm)) + 3*bm^2*exp(bm*c)*(C5*cm*cos(c*cm) + C6*cm*sin(c*cm)) + bm^3*exp(-bm*c)*(C4*cos(c*cm) - C3*sin(c*cm)) - bm^3*exp(bm*c)*(C6*cos(c*cm) - C5*sin(c*cm)) - 3*bm*exp(-bm*c)*(C4*cm^2*cos(c*cm) - C3*cm^2*sin(c*cm)) + 3*bm*exp(bm*c)*(C6*cm^2*cos(c*cm) - C5*cm^2*sin(c*cm)) + C1*am^3*exp(-am*c) - C2*am^3*exp(am*c);
 
    %%%m1=fourier series of 5    "bending momentum(m1) at the boundaries is 5"              (3)
    %%%derivatives found from another code
    
    EQ1=bm^2*exp(bm*c)*(C4*cos(c*cm) + C3*sin(c*cm)) - exp(-bm*c)*(C6*cm^2*cos(c*cm) + C5*cm^2*sin(c*cm)) - exp(bm*c)*(C4*cm^2*cos(c*cm) + C3*cm^2*sin(c*cm)) + bm^2*exp(-bm*c)*(C6*cos(c*cm) + C5*sin(c*cm)) + 2*bm*exp(bm*c)*(C3*cm*cos(c*cm) - C4*cm*sin(c*cm)) - 2*bm*exp(-bm*c)*(C5*cm*cos(c*cm) - C6*cm*sin(c*cm)) + C1*am^2*exp(am*c) + C2*am^2*exp(-am*c);
    EQ2=exp(bm*c)*(C4*cm^4*cos(c*cm) + C3*cm^4*sin(c*cm)) + exp(-bm*c)*(C6*cm^4*cos(c*cm) + C5*cm^4*sin(c*cm)) + 4*bm^3*exp(bm*c)*(C3*cm*cos(c*cm) - C4*cm*sin(c*cm)) - 4*bm^3*exp(-bm*c)*(C5*cm*cos(c*cm) - C6*cm*sin(c*cm)) + bm^4*exp(bm*c)*(C4*cos(c*cm) + C3*sin(c*cm)) + bm^4*exp(-bm*c)*(C6*cos(c*cm) + C5*sin(c*cm)) - 4*bm*exp(bm*c)*(C3*cm^3*cos(c*cm) - C4*cm^3*sin(c*cm)) + 4*bm*exp(-bm*c)*(C5*cm^3*cos(c*cm) - C6*cm^3*sin(c*cm)) + C1*am^4*exp(am*c) + C2*am^4*exp(-am*c) - 6*bm^2*exp(bm*c)*(C4*cm^2*cos(c*cm) + C3*cm^2*sin(c*cm)) - 6*bm^2*exp(-bm*c)*(C6*cm^2*cos(c*cm) + C5*cm^2*sin(c*cm));                               
    eq3=((C*EQ1)-((A)*EQ2));
     
    %%%u1(0,+-b)=0    "deformation in x direction is zero at the x=0 boundary"              (4)
   
    eq4 = C1+C2+C4+C6 ;

    %%%m2(c,y)=0    "bending momentum (m2) at the boundaries is 0"                          (5)
 
    EQ3=3*bm^2*exp(bm*c)*(C3*cm*cos(c*cm) - C4*cm*sin(c*cm)) - exp(-bm*c)*(C5*cm^3*cos(c*cm) - C6*cm^3*sin(c*cm)) - exp(bm*c)*(C3*cm^3*cos(c*cm) - C4*cm^3*sin(c*cm)) + 3*bm^2*exp(-bm*c)*(C5*cm*cos(c*cm) - C6*cm*sin(c*cm)) + bm^3*exp(bm*c)*(C4*cos(c*cm) + C3*sin(c*cm)) - bm^3*exp(-bm*c)*(C6*cos(c*cm) + C5*sin(c*cm)) - 3*bm*exp(bm*c)*(C4*cm^2*cos(c*cm) + C3*cm^2*sin(c*cm)) + 3*bm*exp(-bm*c)*(C6*cm^2*cos(c*cm) + C5*cm^2*sin(c*cm)) + C1*am^3*exp(am*c) - C2*am^3*exp(-am*c);
    EQ4=exp(bm*c)*(C3*cm^5*cos(c*cm) - C4*cm^5*sin(c*cm)) + exp(-bm*c)*(C5*cm^5*cos(c*cm) - C6*cm^5*sin(c*cm)) + 5*bm^4*exp(bm*c)*(C3*cm*cos(c*cm) - C4*cm*sin(c*cm)) + 5*bm^4*exp(-bm*c)*(C5*cm*cos(c*cm) - C6*cm*sin(c*cm)) + bm^5*exp(bm*c)*(C4*cos(c*cm) + C3*sin(c*cm)) - bm^5*exp(-bm*c)*(C6*cos(c*cm) + C5*sin(c*cm)) + 5*bm*exp(bm*c)*(C4*cm^4*cos(c*cm) + C3*cm^4*sin(c*cm)) - 5*bm*exp(-bm*c)*(C6*cm^4*cos(c*cm) + C5*cm^4*sin(c*cm)) + C1*am^5*exp(am*c) - C2*am^5*exp(-am*c) - 10*bm^3*exp(bm*c)*(C4*cm^2*cos(c*cm) + C3*cm^2*sin(c*cm)) - 10*bm^2*exp(bm*c)*(C3*cm^3*cos(c*cm) - C4*cm^3*sin(c*cm)) + 10*bm^3*exp(-bm*c)*(C6*cm^2*cos(c*cm) + C5*cm^2*sin(c*cm)) - 10*bm^2*exp(-bm*c)*(C5*cm^3*cos(c*cm) - C6*cm^3*sin(c*cm));                               
    eq5=((C*EQ3)-((A^2)*EQ4));

    %%% U1(e,d)=-U1(e,-d)     "nature of problem is symmetric"                              (6)
 
  eq6=(  C1*exp(am*c) + C2*exp(-1*am*c) + exp(bm*c)*( C3*sin(cm*c)+C4*cos(cm*c) )+ exp(-1*bm*c)*( C5*sin(cm*c)+C6*cos(cm*c) ) )+...
      (C1*exp(-am*c) + C2*exp(am*c) + exp(-bm*c)*(C4*cos(c*cm) - C3*sin(c*cm)) + exp(bm*c)*(C6*cos(c*cm) - C5*sin(c*cm)));
 
  %%%fourier series for a coefficient
 
  T = (20/(n*pi))*((-1)^((n-1)/2));    
 
  %%% solve 6 variables equations and finding the coefficients
                  
  Q = vpasolve([ (eq1)==T/(5*A) , (eq2)==T/(5*A) , ((pi/(2*d))*eq3)==(T/3) , eq4==0 , ((pi/(2*d))*eq5)==-T/3 , eq6==0 ], [C1,C2,C3,C4,C5,C6]);                 
  
  %%%coefficient matrix assigning
  
  D1 = Q.C1 ;
  D2=  Q.C2 ;
  D3 = Q.C3 ;
  D4 = Q.C4 ;
  D5 = Q.C5 ;
  D6=  Q.C6 ;
   
  %%%finding PHI and U1 and U2
    
  PHI(round(10*x+11),round(10*y+6)) = PHI(round(10*x+11),round(10*y+6))+...
  (  D1*exp(am*x) + D2*exp(-1*am*x) + exp(bm*x)*( D3*sin(cm*x)+D4*cos(cm*x) )+...
  exp(-1*bm*x)*( D5*sin(cm*x)+D6*cos(cm*x) ) ) * sin(m*y) ;
  %%%
  
  U1= U1-...
    (  D1*exp(am*x) + D2*exp(-1*am*x) + exp(bm*x)*( D3*sin(cm*x)+D4*cos(cm*x) )+...
     exp(-1*bm*x)*( D5*sin(cm*x)+D6*cos(cm*x) ) ) * m *cos(m*y) ;
 
  %%%
  
  U2 =U2-...
  ( exp(bm*x)*(D3*cm*cos(cm*x) - D4*cm*sin(cm*x)) + exp(-bm*x)*(D5*cm*cos(cm*x)-...
  D6*cm*sin(cm*x)) + bm*exp(bm*x)*(D4*cos(cm*x) + D3*sin(cm*x)) - bm*exp(-bm*x)*(D6*cos(cm*x)+...
  D5*sin(cm*x)) + D1*am*exp(am*x) - D2*am*exp(-am*x) )* sin(m*y) ;

  %%%
  
end


