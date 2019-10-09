clc
clear all
syms C1 C2 C3 C4 C5 C6 C7 C8 ; 
 %%%%%%values
 C=150 ;
 A=150 ;
 E=15 ;
 d=1;
 c=1;
 %%%%%%%%%%%%matrix
 
PHI=zeros(21,11) ;
U1=zeros(21,11)  ;
U2=zeros(21,11)  ;

%%%%%%%%%%%%
%%%coefficients 
%%%%%%%%%%%%%%%%
 
 for x=-1:0.1:1     
       for y=-.5:.1:.5 
            for n=1:2:11          
%---------------------
m=( (n*pi)/(2*d )) ; 
T24 = A*(m^2)  + C ;
T23 = C*(m^2) + 1 ;
T22 = E *(m^2) + 2 *(m^2) ;       
T21 = (3*(T24^4))/(256*(A^4)) ;
T20 =((T24^2)*T23)/(16*(A^3)) ;
T19 = ((T24^3)/(8*(A^3)))+(T22/A)+((T24*T23)/(2*(A^2))) ;
T18 =(3*(T24^2))/(8*(A^2))-(T23/A) ;
T17 = T21-((m^4)/A)+((T24*T22)/(4*(A^2)))-T20 ;
T16 = (2+E)*(m^2) ;
T15 = sqrt( (27*(T19^4))+ (T16*T17*(T18^4))+(256*(T17^3))-(4*(T18^3)*(T19^2))+(128*(T18^2)*(T17^2))-(144*T18*(T19^2)*T17)) ;
T14 =((T24^3)/(8*(A^3))) + (T16/A)-((T24*T23)/(2*(A^2)));
T13 = T21-((m^4)/A)-T20+((T24*T16)/(4*(A^2))) ;
T12 =(9*(T24^4))/(64*(A^4)) ;
T11 = (3*(T24^2)*(T23))/(4*(A^3)) ;
T10 =((T19^2)/2)+((sqrt(3)*T15)/18)-((4*T18*T17)/3)-((T18^3)/27) ;
T9 = ((T14^2)/2)-((4*T18*T13)/3)+ (sqrt(3)/18)*sqrt((128*(T18^2)*(T13^2))+(27*(T14^4))+(16*(T18^4)*T13)+(256*(T13^3))-...
(4*(T18^3)*(T14^2))-(144*T18*(T14^2)*T13))-((T18^3)/27) ;
T8 = sqrt((6*T18*(T10^(1/3)))+(9*(T10^(2/3)))-T12+((12*(m^4))/A)+(T18^2)-((3*T24*T22)/(A^2))+T11) ;
T7=(6*(T9^(1/6)))*(( (6*T18*(T9^(1/3)))+(9*(T9^(2/3)))-T12+((12*(m^4))/A)+(T18^2)+(T11)-((3*T24*T16)/(A^2)))^(1/4)) ;
T6= 3*sqrt(6)*T19*sqrt(27*(T19^2)+3*sqrt(3)*T15 -72*T17*T18-(2*(T18^3)));
T5 = (T8/(6*(T9^(1/6)))) ;
an= sqrt( (T24/(4*A))-T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                   
bn= sqrt( (T24/(4*A))-T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
cn= sqrt( (T24/(4*A))+T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                                                      
dn= sqrt( (T24/(4*A))+T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
am=(an+bn)/2 ;
bm=(bn-an)/(2*1i) ;
cm=cn/1i ;
dm=dn;
%---------------------

T = (20/(n*pi))*((-1)^((n-1)/2)); %%%%havasam bashe

%%%%%%%%%%%%%%%boundary conditions
%---------------------------------------------------------------------------------------------------
%---------------U1(X=0,Y=d)=0 
eq1 = (C2 + C4 + C6 + C7 + C8) ;
%---------------U1(X=c,Y=d)+U1(X=-c,Y=d)=0
eq2 =2*C6*cos(c*cm) + C7*exp(c*dm) + C7*exp(-c*dm) + C8*exp(c*dm) + C8*exp(-c*dm) + exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) ;
%--------------r(X=c,Y=0)=0 
eq3 =(exp(am*c)*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) );
%---------------r(X=0,Y=0)=0 (definite)
eq4 =C2*am^3 - C4*am^3 - C1*bm^3 - C3*bm^3 - C5*cm^3 + C7*dm^3 - C8*dm^3 + 3*C1*am^2*bm - 3*C2*am*bm^2 + 3*C3*am^2*bm + 3*C4*am*bm^2 ;
%---------------p(x=c,y=0)=P
EQ1=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ2=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) ;
eq5 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
%--------------p(x=c,y=0)=-P
EQ1=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ2=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(-am*c)*(C1*bm^3*cos(bm*c) + C2*bm^3*sin(bm*c)) + exp(am*c)*(C3*bm^3*cos(bm*c) + C4*bm^3*sin(bm*c)) - 3*am^2*exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - 3*am^2*exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) - am^3*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am^3*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + 3*am*exp(-am*c)*(C2*bm^2*cos(bm*c) - C1*bm^2*sin(bm*c)) - 3*am*exp(am*c)*(C4*bm^2*cos(bm*c) - C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(-c*dm) + C8*dm^3*exp(c*dm) + C6*cm^3*sin(c*cm)) ;
eq6 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
% eq6=C*m*(C4*am^3 - C2*am^3 + C1*bm^3 + C3*bm^3 + C5*cm^3 - C7*dm^3 + C8*dm^3 - 3*C1*am^2*bm + 3*C2*am*bm^2 - 3*C3*am^2*bm - 3*C4*am*bm^2) - m*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) + m*(E + 1)*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) ;
%---------------dU2/dx = 0 (x=0,y=d)
% eq7 =(C2*am^2 + C4*am^2 - C2*bm^2 - C4*bm^2 - C6*cm^2 + C7*dm^2 + C8*dm^2 + 2*C1*am*bm - 2*C3*am*bm);
 EQ6=am^2*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) - exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) + am^2*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 2*am*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 2*am*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - C6*cm^2*cos(c*cm) + C7*dm^2*exp(c*dm) + C8*dm^2*exp(-c*dm) - C5*cm^2*sin(c*cm) ;
 EQ7=exp(am*c)*(C2*bm^4*cos(bm*c) + C1*bm^4*sin(bm*c)) + exp(-am*c)*(C4*bm^4*cos(bm*c) + C3*bm^4*sin(bm*c)) + 4*am^3*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 4*am^3*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am^4*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^4*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) - 4*am*exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + 4*am*exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) + C6*cm^4*cos(c*cm) + C7*dm^4*exp(c*dm) + C8*dm^4*exp(-c*dm) - 6*am^2*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 6*am^2*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^4*sin(c*cm) ;
 eq7=-1*((C*EQ6)-((A)*EQ7));
%--------------U2(X=c,Y=d)+U2(X=-c,Y=d)=0
eq8 =exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C7*dm*exp(c*dm) - C7*dm*exp(-c*dm) + C8*dm*exp(c*dm) - C8*dm*exp(-c*dm) - 2*C6*cm*sin(c*cm) ;
%---------------
%%%%%%%%%%%solve the equation

 Q = vpasolve([ eq1==0 , eq2==0 , eq3==T , eq4==0 , eq5==4*T , eq6==4*T , eq7==0 ,eq8==0  ], [C1,C2,C3,C4,C5,C6,C7,C8]);  

   
   D1 = Q.C1 ;
   D2=  Q.C2 ;
   D3 = Q.C3 ;
   D4 = Q.C4 ;
   D5 = Q.C5 ;
   D6=  Q.C6 ;
   D7 = Q.C7 ;
   D8=  Q.C8 ;
   
   %%%%%%%%
% 
%     PHI(round(10*x+11),round(10*y+6)) = PHI(round(10*x+11),round(10*y+6))+...
%     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
%     + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x) ) * sin(m*y) ;
% 
     U1(round(10*x+11),round(10*y+6))  = U1(round(10*x+11),round(10*y+6))-...
     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
    + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x)  ) * m* cos(m*y) ;
     
     U2(round(10*x+11),round(10*y+6)) =U2(round(10*x+11),round(10*y+6))+...
      ( exp(am*x)*(D1*bm*cos(bm*x) - D2*bm*sin(bm*x)) + exp(-am*x)*(D3*bm*cos(bm*x) - D4*bm*sin(bm*x)) +...
      am*exp(am*x)*(D2*cos(bm*x) + D1*sin(bm*x)) - am*exp(-am*x)*(D4*cos(bm*x) + D3*sin(bm*x)) +...
      D5*cm*cos(cm*x) + D7*dm*exp(dm*x) - D8*dm*exp(-dm*x) - D6*cm*sin(cm*x) ) * sin(m*y) ;
                   
            end
       end
 end
 %---------------horizental lines1
  x=-1.2:0.12:1.2;
  j=0 ;
  for i=1
  j=j+1 ;    
  P(1)= plot(x,0.344*((U2(:,i))+(6-j)*0.1),'r');
  axis equal
  hold on
  P(1)= plot(x,0.344*(-(U2(:,i))-(6-j)*0.1),'r');
  end
  hold on
  %---------------vertical lines (ommiting distribution)
  hold on
  for i=0
      
  x1=1.2-(i*0.12) ;
  x2=-(1.2-i*0.12) ;
  y1=-0.344*((U2(i+1,1))+0.5) ;
  y2=0.344*((U2(i+1,1))+0.5);
  P(1)= line([x1 x1], [y1 y2],'Color','red'); 
  hold on
  P(1)= line([x2 x2], [y1 y2],'Color','red'); 
  hold on
  end
  axis([-1.3 1.3 -1.3 1.3])

  
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
clc
clear x y C1 C2 C3 C4 C5 C6 C7 C8
syms C1 C2 C3 C4 C5 C6 C7 C8 ; 
 %%%%%%values
 C=150 ;
 A=150 ;
 E=150 ;
 d=1;
 c=1;
 %%%%%%%%%%%%matrix
 
PHI=zeros(21,11) ;
U1=zeros(21,11)  ;
U2=zeros(21,11)  ;

%%%%%%%%%%%%
%%%coefficients 
%%%%%%%%%%%%%%%%
 
 for x=-1:0.1:1     
       for y=-.5:.1:.5 
            for n=1:2:11          
%---------------------
m=( (n*pi)/(2*d )) ; 
T24 = A*(m^2)  + C ;
T23 = C*(m^2) + 1 ;
T22 = E *(m^2) + 2 *(m^2) ;       
T21 = (3*(T24^4))/(256*(A^4)) ;
T20 =((T24^2)*T23)/(16*(A^3)) ;
T19 = ((T24^3)/(8*(A^3)))+(T22/A)+((T24*T23)/(2*(A^2))) ;
T18 =(3*(T24^2))/(8*(A^2))-(T23/A) ;
T17 = T21-((m^4)/A)+((T24*T22)/(4*(A^2)))-T20 ;
T16 = (2+E)*(m^2) ;
T15 = sqrt( (27*(T19^4))+ (T16*T17*(T18^4))+(256*(T17^3))-(4*(T18^3)*(T19^2))+(128*(T18^2)*(T17^2))-(144*T18*(T19^2)*T17)) ;
T14 =((T24^3)/(8*(A^3))) + (T16/A)-((T24*T23)/(2*(A^2)));
T13 = T21-((m^4)/A)-T20+((T24*T16)/(4*(A^2))) ;
T12 =(9*(T24^4))/(64*(A^4)) ;
T11 = (3*(T24^2)*(T23))/(4*(A^3)) ;
T10 =((T19^2)/2)+((sqrt(3)*T15)/18)-((4*T18*T17)/3)-((T18^3)/27) ;
T9 = ((T14^2)/2)-((4*T18*T13)/3)+ (sqrt(3)/18)*sqrt((128*(T18^2)*(T13^2))+(27*(T14^4))+(16*(T18^4)*T13)+(256*(T13^3))-...
(4*(T18^3)*(T14^2))-(144*T18*(T14^2)*T13))-((T18^3)/27) ;
T8 = sqrt((6*T18*(T10^(1/3)))+(9*(T10^(2/3)))-T12+((12*(m^4))/A)+(T18^2)-((3*T24*T22)/(A^2))+T11) ;
T7=(6*(T9^(1/6)))*(( (6*T18*(T9^(1/3)))+(9*(T9^(2/3)))-T12+((12*(m^4))/A)+(T18^2)+(T11)-((3*T24*T16)/(A^2)))^(1/4)) ;
T6= 3*sqrt(6)*T19*sqrt(27*(T19^2)+3*sqrt(3)*T15 -72*T17*T18-(2*(T18^3)));
T5 = (T8/(6*(T9^(1/6)))) ;
an= sqrt( (T24/(4*A))-T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                   
bn= sqrt( (T24/(4*A))-T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
cn= sqrt( (T24/(4*A))+T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                                                      
dn= sqrt( (T24/(4*A))+T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
am=(an+bn)/2 ;
bm=(bn-an)/(2*1i) ;
cm=cn/1i ;
dm=dn;
%---------------------

T = (20/(n*pi))*((-1)^((n-1)/2)); %%%%havasam bashe

%%%%%%%%%%%%%%%boundary conditions
%---------------------------------------------------------------------------------------------------
%---------------U1(X=0,Y=d)=0 
eq1 = (C2 + C4 + C6 + C7 + C8) ;
%---------------U1(X=c,Y=d)+U1(X=-c,Y=d)=0
eq2 =2*C6*cos(c*cm) + C7*exp(c*dm) + C7*exp(-c*dm) + C8*exp(c*dm) + C8*exp(-c*dm) + exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) ;
%--------------r(X=c,Y=0)=0 
eq3 =(exp(am*c)*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) );
%---------------r(X=0,Y=0)=0 (definite)
eq4 =C2*am^3 - C4*am^3 - C1*bm^3 - C3*bm^3 - C5*cm^3 + C7*dm^3 - C8*dm^3 + 3*C1*am^2*bm - 3*C2*am*bm^2 + 3*C3*am^2*bm + 3*C4*am*bm^2 ;
%---------------p(x=c,y=0)=P
EQ1=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ2=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) ;
eq5 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
%--------------p(x=c,y=0)=-P
EQ1=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ2=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(-am*c)*(C1*bm^3*cos(bm*c) + C2*bm^3*sin(bm*c)) + exp(am*c)*(C3*bm^3*cos(bm*c) + C4*bm^3*sin(bm*c)) - 3*am^2*exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - 3*am^2*exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) - am^3*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am^3*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + 3*am*exp(-am*c)*(C2*bm^2*cos(bm*c) - C1*bm^2*sin(bm*c)) - 3*am*exp(am*c)*(C4*bm^2*cos(bm*c) - C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(-c*dm) + C8*dm^3*exp(c*dm) + C6*cm^3*sin(c*cm)) ;
eq6 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
% eq6=C*m*(C4*am^3 - C2*am^3 + C1*bm^3 + C3*bm^3 + C5*cm^3 - C7*dm^3 + C8*dm^3 - 3*C1*am^2*bm + 3*C2*am*bm^2 - 3*C3*am^2*bm - 3*C4*am*bm^2) - m*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) + m*(E + 1)*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) ;
%---------------dU2/dx = 0 (x=0,y=d)
% eq7 =(C2*am^2 + C4*am^2 - C2*bm^2 - C4*bm^2 - C6*cm^2 + C7*dm^2 + C8*dm^2 + 2*C1*am*bm - 2*C3*am*bm);
 EQ6=am^2*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) - exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) + am^2*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 2*am*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 2*am*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - C6*cm^2*cos(c*cm) + C7*dm^2*exp(c*dm) + C8*dm^2*exp(-c*dm) - C5*cm^2*sin(c*cm) ;
 EQ7=exp(am*c)*(C2*bm^4*cos(bm*c) + C1*bm^4*sin(bm*c)) + exp(-am*c)*(C4*bm^4*cos(bm*c) + C3*bm^4*sin(bm*c)) + 4*am^3*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 4*am^3*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am^4*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^4*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) - 4*am*exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + 4*am*exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) + C6*cm^4*cos(c*cm) + C7*dm^4*exp(c*dm) + C8*dm^4*exp(-c*dm) - 6*am^2*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 6*am^2*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^4*sin(c*cm) ;
 eq7=-1*((C*EQ6)-((A)*EQ7));
%--------------U2(X=c,Y=d)+U2(X=-c,Y=d)=0
eq8 =exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C7*dm*exp(c*dm) - C7*dm*exp(-c*dm) + C8*dm*exp(c*dm) - C8*dm*exp(-c*dm) - 2*C6*cm*sin(c*cm) ;
%---------------
%%%%%%%%%%%solve the equation

 Q = vpasolve([ eq1==0 , eq2==0 , eq3==T , eq4==0 , eq5==4*T , eq6==4*T , eq7==0 ,eq8==0  ], [C1,C2,C3,C4,C5,C6,C7,C8]);  

   
   D1 = Q.C1 ;
   D2=  Q.C2 ;
   D3 = Q.C3 ;
   D4 = Q.C4 ;
   D5 = Q.C5 ;
   D6=  Q.C6 ;
   D7 = Q.C7 ;
   D8=  Q.C8 ;
   
   %%%%%%%%
% 
%     PHI(round(10*x+11),round(10*y+6)) = PHI(round(10*x+11),round(10*y+6))+...
%     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
%     + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x) ) * sin(m*y) ;
% 
     U1(round(10*x+11),round(10*y+6))  = U1(round(10*x+11),round(10*y+6))-...
     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
    + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x)  ) * m* cos(m*y) ;
     
     U2(round(10*x+11),round(10*y+6)) =U2(round(10*x+11),round(10*y+6))+...
      ( exp(am*x)*(D1*bm*cos(bm*x) - D2*bm*sin(bm*x)) + exp(-am*x)*(D3*bm*cos(bm*x) - D4*bm*sin(bm*x)) +...
      am*exp(am*x)*(D2*cos(bm*x) + D1*sin(bm*x)) - am*exp(-am*x)*(D4*cos(bm*x) + D3*sin(bm*x)) +...
      D5*cm*cos(cm*x) + D7*dm*exp(dm*x) - D8*dm*exp(-dm*x) - D6*cm*sin(cm*x) ) * sin(m*y) ;
                   
            end
       end
 end
 %---------------horizental lines1
  x=-1.0493:0.10493:1.0493;
  j=0 ;
  for i=1
  j=j+1 ;    
  P(2)=plot(x,(0.9963*(U2(:,i))+(6-j)*0.083),'b');
  axis equal
  hold on
  P(2)=plot(x,0.9963*(-(U2(:,i))-(6-j)*0.083),'b');
  end
  hold on
  %---------------vertical lines (ommiting distribution)
  hold on
  for i=0
      
  x1=1.0493-(i*0.10493) ;
  x2=-(1.0493-i*0.10493) ;
  y1=-0.9963*(((U2(i+1,1))+0.415)) ;
  y2=0.9963*(((U2(i+1,1))+0.415));
  P(2)=line([x1 x1], [y1 y2],'Color','blue'); 
  hold on
  P(2)= line([x2 x2], [y1 y2],'Color','blue'); 
  hold on

  end
  axis([-1.3 1.3 -1.3 1.3])

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
clc
clear x y C1 C2 C3 C4 C5 C6 C7 C8
syms C1 C2 C3 C4 C5 C6 C7 C8 ; 
 %%%%%%values
 C=150 ;
 A=150 ;
 E=500 ;
 d=1;
 c=1;
 %%%%%%%%%%%%matrix
 
PHI=zeros(21,11) ;
U1=zeros(21,11)  ;
U2=zeros(21,11)  ;

%%%%%%%%%%%%
%%%coefficients 
%%%%%%%%%%%%%%%%
 
 for x=-1:0.1:1     
       for y=-.5:.1:.5 
            for n=1:2:11          
%---------------------
m=( (n*pi)/(2*d )) ; 
T24 = A*(m^2)  + C ;
T23 = C*(m^2) + 1 ;
T22 = E *(m^2) + 2 *(m^2) ;       
T21 = (3*(T24^4))/(256*(A^4)) ;
T20 =((T24^2)*T23)/(16*(A^3)) ;
T19 = ((T24^3)/(8*(A^3)))+(T22/A)+((T24*T23)/(2*(A^2))) ;
T18 =(3*(T24^2))/(8*(A^2))-(T23/A) ;
T17 = T21-((m^4)/A)+((T24*T22)/(4*(A^2)))-T20 ;
T16 = (2+E)*(m^2) ;
T15 = sqrt( (27*(T19^4))+ (T16*T17*(T18^4))+(256*(T17^3))-(4*(T18^3)*(T19^2))+(128*(T18^2)*(T17^2))-(144*T18*(T19^2)*T17)) ;
T14 =((T24^3)/(8*(A^3))) + (T16/A)-((T24*T23)/(2*(A^2)));
T13 = T21-((m^4)/A)-T20+((T24*T16)/(4*(A^2))) ;
T12 =(9*(T24^4))/(64*(A^4)) ;
T11 = (3*(T24^2)*(T23))/(4*(A^3)) ;
T10 =((T19^2)/2)+((sqrt(3)*T15)/18)-((4*T18*T17)/3)-((T18^3)/27) ;
T9 = ((T14^2)/2)-((4*T18*T13)/3)+ (sqrt(3)/18)*sqrt((128*(T18^2)*(T13^2))+(27*(T14^4))+(16*(T18^4)*T13)+(256*(T13^3))-...
(4*(T18^3)*(T14^2))-(144*T18*(T14^2)*T13))-((T18^3)/27) ;
T8 = sqrt((6*T18*(T10^(1/3)))+(9*(T10^(2/3)))-T12+((12*(m^4))/A)+(T18^2)-((3*T24*T22)/(A^2))+T11) ;
T7=(6*(T9^(1/6)))*(( (6*T18*(T9^(1/3)))+(9*(T9^(2/3)))-T12+((12*(m^4))/A)+(T18^2)+(T11)-((3*T24*T16)/(A^2)))^(1/4)) ;
T6= 3*sqrt(6)*T19*sqrt(27*(T19^2)+3*sqrt(3)*T15 -72*T17*T18-(2*(T18^3)));
T5 = (T8/(6*(T9^(1/6)))) ;
an= sqrt( (T24/(4*A))-T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                   
bn= sqrt( (T24/(4*A))-T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)-T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
cn= sqrt( (T24/(4*A))+T5-(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;                                                                      
dn= sqrt( (T24/(4*A))+T5+(sqrt(( (-T8*(T18^2))-(9*T8*(T10^(2/3)))+(12*T17*T8)+T6+(12*(T10^(1/3))*T8*T18))/T7))) ;
am=(an+bn)/2 ;
bm=(bn-an)/(2*1i) ;
cm=cn/1i ;
dm=dn;
%---------------------

T = (20/(n*pi))*((-1)^((n-1)/2)); %%%%havasam bashe

%%%%%%%%%%%%%%%boundary conditions
%---------------------------------------------------------------------------------------------------
%---------------U1(X=0,Y=d)=0 
eq1 = (C2 + C4 + C6 + C7 + C8) ;
%---------------U1(X=c,Y=d)+U1(X=-c,Y=d)=0
eq2 =2*C6*cos(c*cm) + C7*exp(c*dm) + C7*exp(-c*dm) + C8*exp(c*dm) + C8*exp(-c*dm) + exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) ;
%--------------r(X=c,Y=0)=0 
eq3 =(exp(am*c)*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) );
%---------------r(X=0,Y=0)=0 (definite)
eq4 =C2*am^3 - C4*am^3 - C1*bm^3 - C3*bm^3 - C5*cm^3 + C7*dm^3 - C8*dm^3 + 3*C1*am^2*bm - 3*C2*am*bm^2 + 3*C3*am^2*bm + 3*C4*am*bm^2 ;
%---------------p(x=c,y=0)=P
EQ1=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ2=m*(exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(c*dm) - C8*dm*exp(-c*dm) - C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) - 3*am^2*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 3*am^2*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - am^3*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^3*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 3*am*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 3*am*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(c*dm) + C8*dm^3*exp(-c*dm) - C6*cm^3*sin(c*cm)) ;
eq5 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
%--------------p(x=c,y=0)=-P
EQ1=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ2=m*(exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) + exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) - am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + C5*cm*cos(c*cm) + C7*dm*exp(-c*dm) - C8*dm*exp(c*dm) + C6*cm*sin(c*cm)) ;
EQ3=-m*(exp(-am*c)*(C1*bm^3*cos(bm*c) + C2*bm^3*sin(bm*c)) + exp(am*c)*(C3*bm^3*cos(bm*c) + C4*bm^3*sin(bm*c)) - 3*am^2*exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - 3*am^2*exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) - am^3*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am^3*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) + 3*am*exp(-am*c)*(C2*bm^2*cos(bm*c) - C1*bm^2*sin(bm*c)) - 3*am*exp(am*c)*(C4*bm^2*cos(bm*c) - C3*bm^2*sin(bm*c)) + C5*cm^3*cos(c*cm) - C7*dm^3*exp(-c*dm) + C8*dm^3*exp(c*dm) + C6*cm^3*sin(c*cm)) ;
eq6 =((1+E)*EQ1-1*EQ2-C*EQ3)   ;
% eq6=C*m*(C4*am^3 - C2*am^3 + C1*bm^3 + C3*bm^3 + C5*cm^3 - C7*dm^3 + C8*dm^3 - 3*C1*am^2*bm + 3*C2*am*bm^2 - 3*C3*am^2*bm - 3*C4*am*bm^2) - m*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) + m*(E + 1)*(C2*am - C4*am + C1*bm + C3*bm + C5*cm + C7*dm - C8*dm) ;
%---------------dU2/dx = 0 (x=0,y=d)
% eq7 =(C2*am^2 + C4*am^2 - C2*bm^2 - C4*bm^2 - C6*cm^2 + C7*dm^2 + C8*dm^2 + 2*C1*am*bm - 2*C3*am*bm);
 EQ6=am^2*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) - exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) + am^2*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + 2*am*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 2*am*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) - C6*cm^2*cos(c*cm) + C7*dm^2*exp(c*dm) + C8*dm^2*exp(-c*dm) - C5*cm^2*sin(c*cm) ;
 EQ7=exp(am*c)*(C2*bm^4*cos(bm*c) + C1*bm^4*sin(bm*c)) + exp(-am*c)*(C4*bm^4*cos(bm*c) + C3*bm^4*sin(bm*c)) + 4*am^3*exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - 4*am^3*exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am^4*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) + am^4*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) - 4*am*exp(am*c)*(C1*bm^3*cos(bm*c) - C2*bm^3*sin(bm*c)) + 4*am*exp(-am*c)*(C3*bm^3*cos(bm*c) - C4*bm^3*sin(bm*c)) + C6*cm^4*cos(c*cm) + C7*dm^4*exp(c*dm) + C8*dm^4*exp(-c*dm) - 6*am^2*exp(am*c)*(C2*bm^2*cos(bm*c) + C1*bm^2*sin(bm*c)) - 6*am^2*exp(-am*c)*(C4*bm^2*cos(bm*c) + C3*bm^2*sin(bm*c)) + C5*cm^4*sin(c*cm) ;
 eq7=-1*((C*EQ6)-((A)*EQ7));
%--------------U2(X=c,Y=d)+U2(X=-c,Y=d)=0
eq8 =exp(am*c)*(C1*bm*cos(bm*c) - C2*bm*sin(bm*c)) - exp(-am*c)*(C1*bm*cos(bm*c) + C2*bm*sin(bm*c)) - exp(am*c)*(C3*bm*cos(bm*c) + C4*bm*sin(bm*c)) + exp(-am*c)*(C3*bm*cos(bm*c) - C4*bm*sin(bm*c)) + am*exp(am*c)*(C2*cos(bm*c) + C1*sin(bm*c)) - am*exp(-am*c)*(C2*cos(bm*c) - C1*sin(bm*c)) + am*exp(am*c)*(C4*cos(bm*c) - C3*sin(bm*c)) - am*exp(-am*c)*(C4*cos(bm*c) + C3*sin(bm*c)) + C7*dm*exp(c*dm) - C7*dm*exp(-c*dm) + C8*dm*exp(c*dm) - C8*dm*exp(-c*dm) - 2*C6*cm*sin(c*cm) ;
%---------------
%%%%%%%%%%%solve the equation

 Q = vpasolve([ eq1==0 , eq2==0 , eq3==T , eq4==0 , eq5==4*T , eq6==4*T , eq7==0 ,eq8==0  ], [C1,C2,C3,C4,C5,C6,C7,C8]);  

   
   D1 = Q.C1 ;
   D2=  Q.C2 ;
   D3 = Q.C3 ;
   D4 = Q.C4 ;
   D5 = Q.C5 ;
   D6=  Q.C6 ;
   D7 = Q.C7 ;
   D8=  Q.C8 ;
   
   %%%%%%%%
% 
%     PHI(round(10*x+11),round(10*y+6)) = PHI(round(10*x+11),round(10*y+6))+...
%     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
%     + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x) ) * sin(m*y) ;
% 
     U1(round(10*x+11),round(10*y+6))  = U1(round(10*x+11),round(10*y+6))-...
     ( exp(am*x)*( D1*sin(bm*x)+D2*cos(bm*x) )+ exp(-1*am*x)*( D3*sin(bm*x)+D4*cos(bm*x))...
    + ( D5*sin(cm*x)+D6*cos(cm*x))+ D7*exp(dm*x) + D8*exp(-1*dm*x)  ) * m* cos(m*y) ;
     
     U2(round(10*x+11),round(10*y+6)) =U2(round(10*x+11),round(10*y+6))+...
      ( exp(am*x)*(D1*bm*cos(bm*x) - D2*bm*sin(bm*x)) + exp(-am*x)*(D3*bm*cos(bm*x) - D4*bm*sin(bm*x)) +...
      am*exp(am*x)*(D2*cos(bm*x) + D1*sin(bm*x)) - am*exp(-am*x)*(D4*cos(bm*x) + D3*sin(bm*x)) +...
      D5*cm*cos(cm*x) + D7*dm*exp(dm*x) - D8*dm*exp(-dm*x) - D6*cm*sin(cm*x) ) * sin(m*y) ;
                   
            end
       end
 end
 %---------------horizental lines1
  x=-1.0293:0.10293:1.0293;
  j=0 ;
  for i=1
  j=j+1 ;    
  P(3)=plot(x,(0.50*(U2(:,i))+(6-j)*0.097),'k');
  axis equal
  hold on
  P(3)=plot(x,(0.50*(-U2(:,i))-(6-j)*0.097),'k');
  end
  hold on
  %---------------vertical lines (ommiting distribution)
  hold on
  for i=0
      
  x1=1.0293-(i*0.10293) ;
  x2=-(1.0293-i*0.10293) ;
  y1=-0.50*(((U2(i+1,1))+0.485)) ;
  y2=0.50*(((U2(i+1,1))+0.485));
  P(3)=line([x1 x1], [-.5 0.5],'Color','black'); 
  hold on
  P(3)= line([x2 x2], [-.5 0.5],'Color','black'); 
   hold on
  end
  axis([-1.3 1.3 -1.3 1.3])
  legend(P,'E=15','E=150','E=1500')
