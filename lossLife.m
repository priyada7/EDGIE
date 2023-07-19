function [lossUnoptT,lossOptT]=lossLife(HST1,HST2)

for i=1:length(HST2)
lossUnoptT(i,1) = exp(15000/383-15000/(HST1(i,1)+273));
end


fprintf('Loss of life UNOPT: %.5g h.\n',sum(lossUnoptT))

for i=1:length(HST2)
lossOptT(i,1) = exp(15000/383-15000/(HST2(i,1)+273));
end
fprintf('Loss of life OPT: %.5g h.\n',sum(lossOptT))
end