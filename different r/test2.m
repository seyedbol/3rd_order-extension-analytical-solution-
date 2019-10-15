 clc
 clear all
for x =-0.1:0.002:0.1
         for y =-0.5:.01:0.5
             p(round(500*x+51),round(100*y+51))=(-abs(alter(x,y)));

         end
end
%-------------------------------------------------------------------
 for i=1:1:101
     for j=1:1:101
      x(i,j)=2.1*(j-1)*0.01 ; 
     end         
 end
 
 for i=1:1:101
     for j=1:1:101
         u=2.1*((j-1)/100); 
         if(i<51)
       y(i,j)=((3-((3/50)*i))*(- 0.022.*(u-1.05).^4 + 1.3e-18.*(u-1.05).^3 - 0.0038.*(u-1.05).^2 - 2.5e-19.*(u-1.05))+(i-1)/100)*(5/6)-.42 ;
         else
       y(i,j)=(-((3/50)*(i-50))*(- 0.022.*(u-1.05).^4 + 1.3e-18.*(u-1.05).^3 - 0.0038.*(u-1.05).^2 - 2.5e-19.*(u-1.05))+(i-1)/100)*(5/6)-.42;
         end
      end
 end
p=p+2.9 ;

colormap 'jet' ;
q = smoothdata(p);
surf(x,y,q) ;

% caxis([min(min(p))+0.2  max(max(p))])
axis equal
view(0,90);