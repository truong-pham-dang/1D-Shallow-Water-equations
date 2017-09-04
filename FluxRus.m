%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% une sélection de flux numériques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rusanov
function [flux] = FluxRus(u,v)	%Flux(U_R,U_L)
  c = max(Fprime(u),Fprime(v));
  flux = (F(u)+F(v))/2 - c*(u-v)/2;
end
