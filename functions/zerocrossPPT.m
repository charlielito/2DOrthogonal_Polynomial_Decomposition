function time = zerocrossPPT(x)


N = length(x);
time = N;
curr = sign(x(1));
for i = 2:N
    if (curr - sign(x(i)) == 2) %from positive to negative (1 - (-1))
        time = i;
        break;
    end
    curr = sign(x(i));
end
