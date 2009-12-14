function q = q_IF_ONOFF(tk1,tk2, b, d, C, tauf, scale, dt);

q1 = C*d - b*diff(tk1);
fb1 = 0*q1;

for i=1:length(tk1)-1
    ln=length(find(tk2<tk1(i+1)));
    if ln>0
        for j = 1:ln
            fb1(i) = fb1(i) + dt*trapz(cross_fb(max(tk1(i)-tk2(j),0):dt:(tk1(i+1)-tk2(j)),5,7,tauf,scale));
        end
    else
        fb1(i) = 0;
    end
end

q1 = q1 - fb1;

q2 = -C*d + b*diff(tk2);
fb2 = 0*q2;

for i=1:length(tk2)-1
    ln=length(find(tk1<tk2(i+1)));
    if ln>0
        for j = 1:ln
            fb2(i) = fb2(i) + dt*trapz(cross_fb(max(tk2(i)-tk1(j),0):dt:(tk2(i+1)-tk1(j)),5,7,tauf,scale));
        end
    else
        fb2(i) = 0;
    end
end

q2 = q2 + fb2;

q = [q1, q2]';