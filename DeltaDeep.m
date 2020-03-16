function [ DeltaCombo ] = DeltaDeep(a,Deltas )
%Calculates the average upwelling C14 from the input N/S concentration and
%the deep water concentration, where a(1) is the deep water component, and
%a(2) is the ratio of deep water upwelled

DeltaCombo = a(2)*a(1)*ones(1,length(Deltas))+(1-a(2))*Deltas;
end

