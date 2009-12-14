function r = r_pop_IF(TK,ln,N)

ln2 = cumsum([0,ln]);
r = zeros(sum(ln),1);

for i = 1:N
    r(ln2(i)+1:ln2(i+1)) = diff(TK(1:ln(i)+1,i).^2)/2;
end