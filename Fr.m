%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Froude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [froude] = Fr(U)
  global g;
  if ( U(1) == 0 )
    froude = 0.;
  else
    froude = U(2)/(U(1)*sqrt(g*U(1)));
  end
end 