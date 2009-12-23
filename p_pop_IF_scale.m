function p = p_pop_IF_scale(TK,ln,w,N)

ln2 = cumsum([0,ln]);
p = zeros(sum(ln),1);

for i = 1:N
    p(ln2(i)+1:ln2(i+1)) = w(i)*diff(TK(1:ln(i)+1,i));
end