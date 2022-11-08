%%THIS FILE IS FOR PART E

function [iL, Vload] = buckedited(L, C, R, delta_t, z, s, tend)

A_inverse = [1 (delta_t/(L)); (-delta_t/(C)) (1+(delta_t/(R*C)))]^-1;

k = 0;
i = 1;

iL(1) = z;
Vload(1) = s;

while(k <= tend)
    sw(i) = switching(1/10, k);
    C = [delta_t*sw(i)/L; 0]*(480 - 3);
    p = A_inverse * ([iL(i); Vload(i)] + C);
    iL(i+1) = p(1);
    Vload(i+1) = p(2);
    k = k + (delta_t);
    i = i + 1;
end

iL = iL .* (iL >= 0);


end