function f=hstar(g,u,h0)
p=[g,-g*h0,-(g*h0^2+2*h0*u^2),g*h0^3];
f=roots(p);
end