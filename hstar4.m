function h=hstar4(hl,hr)
p=[1,0,-9*hl,16*hl*sqrt(hr),-hl^2-8*hl*hr,0,hl^3];
t=roots(p);
for i=1:length(t)
    if imag(t(i))==0
        if and(t(i)<sqrt(hr),t(i)>sqrt(hl))
            h=t(i)^2;
        end
    end
end
%u=ur-2*sqrt(g)*(sqrt(hr)-sqrt(h));
end

