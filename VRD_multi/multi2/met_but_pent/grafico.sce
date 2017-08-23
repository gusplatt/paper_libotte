[nlin,ncol]=size(maty); 

linha = input('entre o valor da linha:');

xset('window',0)
plot(maty(linha,:),matp(linha,:),'k-')
xset('window',1)
plot(maty2(linha,:),matp(linha,:),'k-')
xset('window',2)
plot(1-maty(linha,:)-maty2(linha,:),matp(linha,:),'k-')

passo_mx = (matx(linha,ncol) - matx(linha,ncol-1))
passo_mx2 = (matx2(linha,ncol) - matx2(linha,ncol-1))

