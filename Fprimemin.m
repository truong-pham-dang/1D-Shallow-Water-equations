function [v] = Fprimemin(U)  %calcule le max des valeurs propres
  global g;
  if ( U(1)==0. )
    v = 0.;
  else
    v = min((U(2)/U(1)+sqrt(g*U(1))),(U(2)/U(1)-sqrt(g*U(1))));
  end
end