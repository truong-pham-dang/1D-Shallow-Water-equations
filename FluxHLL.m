%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% une sélection de flux numériques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HLL
function [flux] = FluxHLL(u,v)	%Flux(U_R,U_L)
  cR = max(Fprimemax(u),Fprimemax(v));
  cL = min(Fprimemin(u),Fprimemin(v));
  if 0<cL
      flux = F(v);
  elseif cR<0
      flux = F(u);
  else
      flux = (cR*F(v)-cL*F(u) + cL*cR*(u-v))/(cR-cL);
  end
  %flux = (F(u)+F(v))/2 - c*(u-v)/2;
end