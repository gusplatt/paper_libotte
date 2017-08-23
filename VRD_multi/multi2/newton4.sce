mode(-1)
lines(0)
warning off

function [f,penal] = Fi(teta,N,bbeta,flag)

x1 = teta(1);
x = [x1 ; 1-x1];
y = [teta(3);teta(4)];
P = teta(5);      

//% ----------     Dados de entrada para  Methane + n-butane      -------- %

R   = 8.31434e-3;                                //% Constante universal dos gases em kPa*m3/mol*K
Tc1 = 190.56;                                     //% Temperatura critica do Metano em Kelvin
Tc2 = 425.12;                                     //% Temperatura critica do Butano em Kelvin
Pc1 = 4599;                                      // Pressão crítica do metano em kPa
Pc2 = 3796;                                  //% Pressao critica do Butano em kPa
w1  = 0.011;                                     //% Fator acentrico para o Metano
w2  = 0.200;                                     //% fator acentrico para o Butano
T   = 189.06; //% Temperatura em Kelvin (observe-se que T<Tc do metano)


Tr1 = T/Tc1;
Tr2 = T/Tc2;


k1 = 0.37464 + (1.5422*w1) - ((0.26992)*(w1^2));
k2 = 0.37464 + (1.5422*w2) - ((0.26992)*(w2^2));

alfa1 = (1+k1*(1-(Tr1^(1/2))))^2;
alfa2 = (1+k2*(1-(Tr2^(1/2))))^2;


a1 = (0.45724*(R^2)*(Tc1^2)/Pc1) * alfa1;
a2 = (0.45724*(R^2)*(Tc2^2)/Pc2) * alfa2;

k11 = 0;
k22 = 0;
k12 = -0.04186;
k21 = -0.04186;

// modificado em 21/12/06
k12 = 0;
k21 = 0;


a11 = sqrt(a1*a1)*(1-k11);
a12 = sqrt(a1*a2)*(1-k12);
a21 = sqrt(a2*a1)*(1-k21);
a22 = sqrt(a2*a2)*(1-k22);

b1 = 0.07780*R*Tc1/Pc1;
b2 = 0.07780*R*Tc2/Pc2;


if flag == 0
    ccomp = x;
else
    ccomp = y;
end

//% ---------                  Regra de Lorentz                --------- %

//%   i=1  
b11 = ((((b1^(1/3))+(b1^(1/3)))/2)^3);
b12 = ((((b1^(1/3))+(b2^(1/3)))/2)^3) * (1-0.07484);


//%   i=2  
b21 = ((((b2^(1/3))+(b1^(1/3))))/2)^3* (1-0.07484);
b22 = ((((b2^(1/3))+(b2^(1/3))))/2)^3;

// regra de combinação clássica

b11 = b1;
b12 = (b1+ b2)/2;
b21 = b12;
b22 = b2;


//% ---------                Regras de mistura                --------- %

a = ccomp(1)^2*a11 + 2*ccomp(1)*ccomp(2)*a12 + ccomp(2)^2*a22; 
b = ccomp(1)^2*b11 + 2*ccomp(1)*ccomp(2)*b12 + (ccomp(2)*ccomp(1)*b21) + ccomp(2)^2*b22;


A = (a*P)/(R^2*T^2);
B = (b*P)/(R*T);

//% ---------                      Derivadas                   ----------- %

//%%  i=1, i=2 
dA1 = ccomp(1)*b11 + ccomp(2)*b12;
dA2 = ccomp(1)*b21 + ccomp(2)*b22;

//%%  i=1, i=2   
dB1 = 2*ccomp(1)*a11 + 2*ccomp(2)*a12;
dB2 = 2*ccomp(1)*a21 + 2*ccomp(2)*a22;


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


//% ---------              Coeficiente de fugacidade              ----------%

//% Fi do componente 1
f1=exp((-log ((P*Vol)/(R*T) - (P*b)/(R*T)))+((P*b)/(R*T)) + (1/b*dA1*(((P*Vol)/(R*T))-1))+...
      ((1/(2*sqrt(2)))*(a/(b*R*T))*((1/a)*dB1-(1/b)*dA1))*log((Vol+b*(1-sqrt(2)))/(Vol+b*(1+sqrt(2)))));

//% Fi do componente 2
f2=exp((-log ((P*Vol)/(R*T) - (P*b)/(R*T)))+((P*b)/(R*T)) + (1/b*dA2*(((P*Vol)/(R*T))-1))+...
      ((1/(2*sqrt(2)))*(a/(b*R*T))*((1/a)*dB2-(1/b)*dA2))*log((Vol+b*(1-sqrt(2)))/(Vol+b*(1+sqrt(2)))));

f = [f1;f2];

endfunction

function g = ef(teta,N,bbeta)

g1 = (N - (1-bbeta)*teta(1:2) - bbeta*teta(3:4));
[Fiv,penalL] = Fi(teta,N,bbeta,1);
[Fil,penalV] = Fi(teta,N,bbeta,0);
g2 = Fiv(1:2).*teta(3:4) - Fil(1:2).*teta(1:2);
g3 = teta(3) + teta(4) - 1;
g  = [g1 ; g2 ; g3];

endfunction

function jj = jacobiana(teta,N,bbeta);

jac  = [];

for k = 1:5
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
vx = [];


for N1 = 0.0001:0.0001:0.025
N = [N1 ; 1- N1];
teta0 = [N1;1-N1;.01;.99;100] // trecho de pressão baixa
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
vx = [vx ; N1];

end

plot(vy1baixo,vpbaixo,'r-')
plot(vx,vpbaixo,'b-')

