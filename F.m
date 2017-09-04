function [Flux] = F(U)
  global g;
  Flux(1)=U(2);
  if ( U(1)==0. )
    Flux(2) = 0.;
  else
    Flux(2)=U(2)*U(2)/U(1)+g*U(1)*U(1)/2;
  end
end