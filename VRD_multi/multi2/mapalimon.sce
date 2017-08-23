// constrói o mapa de raízes
// composição y1 = 0.998966 (fração molar de etano na fase vapor)

[p,q]=size(matp)

mataux = matp;

for linha = 1:p
	for coluna = 1:q
			elemento = matp(linha,coluna);
			if elemento == -1
				// ponto não convergido
				mataux(linha,coluna)=1;
				elseif abs(elemento - 4859.5)<1
				// raiz "1"
				mataux(linha,coluna)=2;
				elseif abs(elemento - 4931.2)<1
				// raiz "2"
				mataux(linha,coluna)=3;
				elseif abs(elemento - 619.24)<1
				// raiz "3"
				mataux(linha,coluna)=24;
				elseif abs(elemento - 5007.3)<1
				// raiz "4"
				mataux(linha,coluna)=5;
				elseif elemento == -2
				//solução trivial
				mataux(linha,coluna)=7;
			end
	end
end

Matplot(mataux)
