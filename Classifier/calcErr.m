function err = calcErr(T, Y)
    err = mean(sign(T)~=sign(Y));
end