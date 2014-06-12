function y = regularizedKernLSTest(w, Xtr, kernel, sigma, Xts)
    Ktest = KernelMatrix(Xts, Xtr, kernel, sigma);
    y = Ktest*w;
end
