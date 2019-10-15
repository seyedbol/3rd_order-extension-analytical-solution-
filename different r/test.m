clc
clear all
ss=10;
gg=0.71;
for x =0.2/gg:-.004/gg:0
         for y =-0.5:.01:0.5
             p(round(gg*(10/4)*100*x+51),round(100*y+51))=ss*(-abs(alter1(0.55/gg-x,y)))+0.2;

         end
end
for x =-.004/gg:-.004/gg:-0.2/gg
         for y =-0.5:.01:0.5
             p(round(gg*(10/4)*100*x+51),round(100*y+51))=ss*(-abs(alter1(0.55/gg+x,y)))+0.2;

         end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
pp=130;
gg=1;
for x =0.2/gg:-.004/gg:0
         for y =-0.5:.01:0.5
             if(p(round(gg*(250)*x+51),round((100)*y+51))<.174)
             p(round(gg*(250)*x+51),round((100)*y+51))=p(round(gg*(250)*x+51),round((100)*y+51))+ -pp*(abs(alter2(0.5/gg-x,y)))+0.01;
             end    
         end
end
for x =-.004/gg:-.004/gg:-0.2/gg
         for y =-0.5:.01:0.5
            if(p(round(gg*(250)*x+51),round((100)*y+51))<.174)
             p(round(gg*(250)*x+51),round((100)*y+51))=p(round(gg*(250)*x+51),round((100)*y+51))+ -pp*(abs(alter2(0.5/gg+x,y)))+0.01;
            end
         end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=p+0.25 ;


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
colormap 'jet' ;
q = smoothdata(p);

surf(x,y,q) ;

% caxis( [( min(min(p))+max(max(p)) )/2.8 max(max(p))-.003])
axis equal
view(0,90);
