function G = G_pop_IF_scale(TK,ln,w,N)

G = zeros(sum(ln),sum(ln));
ln2 = cumsum([0,ln]);

for i = 1:N
    for j = 1:N
        Gb = Gblock_IF(TK(1:ln(i)+1,i)',TK(1:ln(j)+1,j)');
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = w(i)*w(j)*Gb;
    end
end