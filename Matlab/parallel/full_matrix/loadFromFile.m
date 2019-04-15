function [A,b,init,nRLC] = loadFromFile()
    A = load('Koko_exact_A_1925.mat');
    b = load('Koko_exact_b_1925.mat');
    init = load('Koko_exact_ic_1925.mat');
    
    A = A.A;
    b = b.b_const;
    init = init.y0;
    
    nRLC = size(A,1);
end