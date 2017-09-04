function sigma = HancockSlopes(a,b,c,delta)
n = length(a);
%% Three slopes
x = (b-a)/delta;
y = (c-a)/(2*delta);
z = (c-b)/delta;
%% Minmod function
sigma = zeros(n,1);
for i = 1:n
    if min(min(x(i),y(i)),z(i))>0
        sigma(i) = min(min(x(i),y(i)),z(i));
    elseif max(max(x(i),y(i)),z(i))<0
        sigma(i) = max(max(x(i),y(i)),z(i));
    else
        sigma(i) = 0;
    end
end
end