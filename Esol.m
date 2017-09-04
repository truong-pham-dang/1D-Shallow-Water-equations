function f=Esol(x,t,g,ul,ur,hl,hr)
%case 1 ur=-ul<0, hl=hr)
if and(and(ur==-ul,ur<0),hl==hr)
    h=max(hstar(g,ul,hl));
    if (x/t)<=(ul-h*sqrt(g*(h+hl)/(2*h*hl)))
        f(1)=hl;
        f(2)=hl*ul;
    elseif (x/t)>=hr*sqrt(g*(h+hr)/(2*h*hr))
        f(1)=hr;
        f(2)=hr*ur;
    else f(1)=h;
        f(2)=0;
    end
end
%case 2 ur=-ul>0, hl=hr)
if and(and(ur==-ul,ur>0),hl==hr)
   h=(ul+2*sqrt(g*hl))^2/(4*g);
   u=0;
   aML=ul-sqrt(g*hl);
   aMS=-sqrt(g*h);
   aPS=sqrt(g*h);
   aPR=ur+sqrt(g*hr);
   if (x/t)<=aML
       f(1)=hl;
       f(2)=hl*ul;
   elseif and((x/t)>aML,(x/t)<=aMS)
       f(1)=1/g*(-1/3*x/t+1/3*ul+2/3*sqrt(g*hl))^2;
       f(2)=f(1)*(2/3*x/t+1/3*ul+2/3*sqrt(g*hl));
   elseif and(aMS<(x/t),(x/t)<=aPS)
       f(1)=h;
       f(2)=h*u;
   elseif and((x/t)>aPS,(x/t)<=aPR)
       f(1)=1/g*(1/3*x/t+1/3*u+2/3*sqrt(g*h))^2;
       f(2)=f(1)*(2/3*x/t+1/3*u-2/3*sqrt(g*h));
   else
       f(1)=hr;
       f(2)=hr*ur;
   end
end
   

%case 3 ur=ul=0, hl>hr)
if and(and(ur==ul,ur==0),hl>hr)
    h=hstar3(hl,hr);
    u=ul+2*sqrt(g)*(sqrt(hl)-sqrt(h));
    aML=ul-sqrt(g*hl);
    aMS=u-sqrt(g*h);
    sigma=u+hr*sqrt(g/2*(1/h+1/hr));
    if (x/t)<=aML
        f(1)=hl;
        f(2)=hl*ul;
    elseif and((x/t)>aML,(x/t)<=aMS)
        f(1)=(-x/t+ul+2*sqrt(g*hl))^2/(9*g);
        f(2)=f(1)*(2*x/t+ul+2*sqrt(g*hl))/(3);
    elseif and((x/t)>aMS,(x/t)<=sigma)
        f(1)=h;
        f(2)=h*u;
    elseif (x/t)>sigma
        f(1)=hr;
        f(2)=hr*ur;
    end    
end
%case 4 ur=ul=0, hl<hr)
if and(and(ur==ul,ur==0),hl<hr)
    h=hstar4(hl,hr);
    u=ur-2*sqrt(g)*(sqrt(hr)-sqrt(h));
    aPR=ur+sqrt(g*hr);
    aPS=u+sqrt(g*h);
    sigma=ul-h*sqrt(g/2*(1/h+1/hl));
    if (x/t)<=sigma
        f(1)=hl;
        f(2)=hl*ul;
    elseif and((x/t)>sigma,(x/t)<=aPS)
        f(1)=h;
        f(2)=h*u;
    elseif and((x/t)>aPS,(x/t)<=aPR)
        f(1)=(x/t-u+2*sqrt(g*h))^2/(9*g);
        f(2)=f(1)*(2*x/t+u-2*sqrt(g*h))/(3);
    elseif (x/t)>aPR
        f(1)=hr;
        f(2)=hr*ur;
    end    
end
end


    