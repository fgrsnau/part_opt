function B = cumwin2(A,W)
% B = cumwin2(A,W)
% sums input A over all windows of size W x W
% values outside A are assumed 0

B = cumwin(cumwin(A,W),W);

end