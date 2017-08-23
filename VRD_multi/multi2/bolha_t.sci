mode(-1)
lines(0)
warning off

function [f,penal] = Fi(teta,N,bbeta,flag)

c = (length(teta)-1)/2;
x = teta(1:c);
y = teta((c+1:2*c));
P = teta(2*c+1);
    

//% ----------     Dados de entrada para  Methane + n-butane      -------- %

R   = 8.31434e-3;                                //% Constante universal dos gases em kPa*m3/mol*K
Tc = [190.56 ; 425.12];
Pc = [4599   ; 3796];
w  = [0.011  ; 0.200];
T   = 189.06; //% Temperatura em Kelvin (observe-se que T<Tc do metano)


Tr = T./Tc;
k  = 0.37464 + (1.5422*w) - ((0.26992)*(w.^2));
alfa  = (1+k.*(1-(Tr.^(1/2)))).^2;
ai  = (0.45724*(R^2)*(Tc.^2)./Pc).*alfa;
bi  = 0.07780*R*Tc./Pc;
k  = zeros(c,c);
M  = sqrt(ai*ai').*(1-k);

if flag == 0
    ccomp = x;
else
    ccomp = y;
end

//% ---------                Regras de mistura                --------- %

a = ccomp'*M*ccomp; 
b = ccomp'*bi;



A = (a*P)/(R^2*T^2);
B = (b*P)/(R*T);

//% ---------                      Derivadas                   ----------- %

//%%  i=1, i=2 
dA = bi;
dB = 2*sqrt(ai*ai')*ccomp

//%----------                        Raizes                       ----------%

Z = poly([-(A*B-B^3-B^2) A-2*B-3*B^2 -(1-B) 1],"z","coeff") // Eq. de Peng-Robinson na forma cúbica

z = roots(Z);                  // fator de compressibilidade
aux = (z == real(z));          // verdadeiro para raízes reais e falso para complexas
aux = aux.*z;                  // multiplicação elemento-a-elemento de modo a eliminar raízes complexas
aux2 = find(aux~=0);           // localiza posições onde há raízes reais
z = real(z(aux2));             // seleciona raízes reais
                 
if length(z)==1
    penal = 1;
else
    penal = 0;
end

if flag==0 // seleciona raiz para a fase líquida
    z = min(real(z));
elseif flag==1 // seleciona raiz para a fase vapor
    z = max(real(z));
end

Vol = z*R*T/P;

//% ---------              Coeficiente de fugacidade              ----------

f=exp((-log ((P*Vol)/(R*T) - (P*b)/(R*T))) + (1./b.*dA*(((P*Vol)/(R*T))-1))+...
      ((1/(2*sqrt(2)))*(a./(b*R*T))*((1./a).*dB-(1./b).*dA))*log((Vol+b*(1-sqrt(2)))./(Vol+b*(1+sqrt(2)))));

endfunction

function g = ef(teta,N,bbeta)
c = (length(teta)-1)/2;
g1 = (N - (1-bbeta)*teta(1:c) - bbeta*teta((c+1):(2*c)));
[Fiv,penalL] = Fi(teta,N,bbeta,1);
[Fil,penalV] = Fi(teta,N,bbeta,0);
g2 = Fiv.*teta((c+1):(2*c)) - Fil.*teta(1:c);
g3 = sum(teta((c+1):(2*c))) - 1;
g  = [g1 ; g2 ; g3];

endfunction

function jj = jacobiana(teta,N,bbeta);
c = (length(teta)-1)/2;
jac  = [];
for k = 1:(2*c+1)
base = ef(teta,N,bbeta)
teta(k) = teta(k) + 1e-5;
avan = ef(teta,N,bbeta);
teta(k) = teta(k) - 1e-5;
der = (avan - base)/1e-5;
jac = [jac der];
end
jj = jac;
endfunction

// programa principal


vpalto = [];
vy1alto = [];
bbeta = 0; // calcula ponto de orvalho

   

vpbaixo = [];
vy1baixo = [];
i=1; // contador de iterações


for N1 = 0.2:0.01:0.939
N = [N1 ; 1- N1];
if i==1
teta0 = [N1;1-N1;.01;.99;100] // trecho de pressão baixa
end
erro = 100;
cont = 0;
while (erro > 1e-8)
   valor = ef(teta0,N,bbeta);
   jacob = jacobiana(teta0,N,bbeta);
   testaj  = max(isnan(jacob)*ones(5,5)); // testa valores não numéricos na matriz jacobiana
   testaj2 = max(isinf(jacob)*ones(5,5)); // testa valores infinitos na matriz jacobiana
   if (testaj == 0)&(testaj2 == 0)
     if rank(jacob) == 5
       flanewton = 0;
       passo = 1;
       while flanewton == 0 
         novoteta = teta0 - passo*real(inv(jacob)*valor)
         if (novoteta(1)>0)&(novoteta(2)>0)&(novoteta(3)>0)&(novoteta(4)>0)&(novoteta(5)>0)&(novoteta(1)<1)&(novoteta(2)<1)&(novoteta(3)<1)&(novoteta(4)<1)
           flanewton = 1;
           //disp('ponto aceito')
         else
           passo = passo/2;
           disp('Redução de passo no Newton-Raphson')
         end
       end
       erro = norm(novoteta - teta0)
       teta0 = novoteta;
       cont = cont + 1;
     else
       disp('Zebra.Problema singular.'),erro = -1;
     end
   else
     disp('Zebra. Valor não numérico na matriz jacobiana.'), erro = -1;
   end
   if cont>= 200
     erro = -1; // ultrapassou o número máximo de iterações do MNR
   end
end


vpbaixo = [vpbaixo ; novoteta(5)]
vy1baixo = [vy1baixo ; novoteta(3)];


end

plot(vy1baixo,vpbaixo,'r-')

vpbaixo = [];
vy1baixo = [];


