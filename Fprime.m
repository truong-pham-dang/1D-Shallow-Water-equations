function [v] = Fprime(U)  %calcule le max des valeurs propres
  global g;
  if ( U(1)==0. )
    v = 0.;
  else
    v = max(abs(U(2)/U(1)+sqrt(g*U(1))),abs(U(2)/U(1)-sqrt(g*U(1))));
  end
end