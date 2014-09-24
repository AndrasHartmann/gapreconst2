function [stable] = checkstable(A);
% Simple stability test of IIR filter coeffitients in A
stable = 1;
for i = 1:size(A,2)
    stable = stable && checkstable([1; A(:,i)]);
end;



