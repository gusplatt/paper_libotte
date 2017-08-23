[p,q]=size(matp)

mataux = matp;

for linha = 1:p
	for coluna = 1:q
			elemento = matp(linha,coluna);
			if elemento == -1
				// ponto não convergido
				mataux(linha,coluna)=1;
				elseif abs(elemento - 1583.65)<1
				// raiz "1"
				mataux(linha,coluna)=2;
				elseif abs(elemento - 3075.16)<1
				// raiz "2"
				mataux(linha,coluna)=3;
				elseif abs(elemento - 4282.97)<1
				// raiz "3"
				mataux(linha,coluna)=24;
				elseif elemento == -2
				//solução trivial
				mataux(linha,coluna)=7;
			end
	end
end

Matplot(mataux)
