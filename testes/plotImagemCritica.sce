clear
clc

file1 = "imagem1finalcompleta.txt";
IC1 = fscanfMat(file1, "%lg");
file2 = "imagem2finalcompleta.txt";
IC2 = fscanfMat(file2, "%lg");

// Use "o" para plotar os pontos discretos ou "-" para a curva
plot(IC1(:, 1), IC1(:, 2), "ro"); //Vermelho (r)
set(gca(), "auto_clear", "off");
plot(IC2(:, 1), IC2(:, 2), "ko"); //Preto (k)
