function IntegerVector = DIG(N,M)
%  Different Integers Generator
%  Generate K different integers in the range of [1:N]

LoveVector = [1:N];
IntegerContainer = zeros(1,M);
for k = 1:M
    index = randint(1,1,[1,length(LoveVector)]);
    IntegerContainer(k) = LoveVector(index);
    LoveVector(index) = [];
end
IntegerVector = sort(IntegerContainer);