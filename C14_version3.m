function [GalModeled] = C14_version3(a,x)

vEUC = x(:,1);
rDeepClim = x(:,2);
AtmC14 = x(:,3);
EUCc = x(:,4);

DeepC = Delta14ToConcentration(-100,2);

Gal0 = 2.270086329970849e-04;
GalModeled = Gal0;
for i=1:length(vEUC)-1
    upC = rDeepClim(i)*DeepC+EUCc(i)*(1-rDeepClim(i));
    AtmIn(i) = (AtmC14(i) - GalModeled(i))*a;
    GalModeled(i+1) = (GalModeled(i)+AtmIn(i))*(1-vEUC(i))+vEUC(i)*upC;
end

end

