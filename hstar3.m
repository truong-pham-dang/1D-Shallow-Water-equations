function h=hstar3(hl,hr)
p=[1,0,-9*hr,16*hr*sqrt(hl),-hr^2-8*hr*hl,0,hr^3];
t=roots(p);
for i=1:length(t)
    if imag(t(i))==0
        if and(t(i)<sqrt(hl),t(i)>sqrt(hr))
            h=t(i)^2;
        end
    end
end
%u=ul+2*sqrt(g)*(sqrt(hl)-sqrt(h));
end

