 clc
 clear all
 
%%%assigning a value of strain to every point in xy_plane

for x =0:0.01:1
         for y = 0:1/110:1-(1/110)
             p(round(100*x+1),round(110*y+1))=(-(strain_countour_at_xy_3(x,y)));
         end
end

%%%fiber meshing and boundaries

for i=1:1:101
    for j=1:1:110  
        x(i,j)=(j-1)*0.01-((i-1)*j/(68000)) ;
        if(x(i,j)<0)
           x(i,j)=0 ;
        end    
    end
end
 
for i=1:1:101
    for j=1:1:110
         u=(j-1)*0.01 ;         
         y(i,j)=0.14*u.^3+.0021*u.^2+0.02*u+((i-1)/100) ;
    end
end
%%%combing the shape of fiber with the corresponding contour

colormap 'jet'
surf(x,y,p) ;
view(0,90);
