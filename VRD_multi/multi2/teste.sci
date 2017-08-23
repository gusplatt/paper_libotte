vn1 = [];
vn2 = [];
matn1 = [];
matn2 = [];
N1i = 0.90;
passo1 = 0.01;
N1f    = 0.98;
N2i    = 0.001;
valo   = (N1f-N1i)/passo1;

for N1 = N1i:passo1:N1f
  passo2 = (1-N1-N2i)/valo;
  for N2 = N2i:passo2:(1-N1+1e-6)
     vn1 = [vn1 ; N1];
     vn2 = [vn2 ; N2];
end
//size(vn1),pause

matn1 = [matn1 vn1];
matn2 = [matn2 vn2];
vn1 = [];
vn2 = [];
end

