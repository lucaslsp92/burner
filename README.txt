Criado por Eric Henrique Favero

Caso deseje gerar GRIDs com mais do que 20 milhões de partículas, altere o #define NUMMAXPARTICLES 20000000 para o valor desejado.

O limite rodando em Windows é de 120 milhões, no linux o limite é bem superior, passando dos 500 milhões.

A fim de conhecimento os arquivos que devem ser editados para gerar a geometria desejado são:
	- ioside.cpp --> define os volumes (fluidos ou vazios) e gera uma casca sólida (parede e dummies) no cortorno do volume definido
	- ioflow.cpp --> realiza operações apenas sobre a casca, returne id caso não deseje fazer nenhuma operação, id + X, para transformar o id da parede e retorne -1 para retirar a casca na região desejada.
			- Obs. todo o código em ioflow não faz nenhuma alteração sobre as regiões fluídas e vazias declaradas no ioside
			
O arquivo de entrada padrão é o data/burner.txt e sua estrutura é:

# modo, int --> dimensão do grid, 2 para casos 2D e 3 para casos 3D, caso o caso seja 2D será considerado apenas o plano xy (z = 0) das regiões definidas nos arquivos ioside e ioflow
# dp, double --> distância entre partículas

# numero de parametros geometria, int --> define o tamanho da lista acessivel dentro dos arquivos ioside e ioflow sobre o nome geometry[]
# lista com n parametros, double

# numero de materiais fluidos, int --> define a quantidade de materiais fluidos
# lista com n ids dos materiais fluidos, int

# numero de materias sólidos, int --> define a quantidade de materias sólidos (sólido, sólido livre e sólido forçado)
# lista dos n ids dos materias sólidos, int --> entre apenas com o valor do id da parede, a dummy equivalente será o id + 1

# ponto a partir do qual o burner começará a procurar as partículas, double double double --> deve coincidir com um ponto onde haja partículas, caso contrário o programa acaba sem gerar nenhuma partícula.