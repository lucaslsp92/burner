Criado por Eric Henrique Favero

Caso deseje gerar GRIDs com mais do que 20 milh�es de part�culas, altere o #define NUMMAXPARTICLES 20000000 para o valor desejado.

O limite rodando em Windows � de 120 milh�es, no linux o limite � bem superior, passando dos 500 milh�es.

A fim de conhecimento os arquivos que devem ser editados para gerar a geometria desejado s�o:
	- ioside.cpp --> define os volumes (fluidos ou vazios) e gera uma casca s�lida (parede e dummies) no cortorno do volume definido
	- ioflow.cpp --> realiza opera��es apenas sobre a casca, returne id caso n�o deseje fazer nenhuma opera��o, id + X, para transformar o id da parede e retorne -1 para retirar a casca na regi�o desejada.
			- Obs. todo o c�digo em ioflow n�o faz nenhuma altera��o sobre as regi�es flu�das e vazias declaradas no ioside
			
O arquivo de entrada padr�o � o data/burner.txt e sua estrutura �:

# modo, int --> dimens�o do grid, 2 para casos 2D e 3 para casos 3D, caso o caso seja 2D ser� considerado apenas o plano xy (z = 0) das regi�es definidas nos arquivos ioside e ioflow
# dp, double --> dist�ncia entre part�culas

# numero de parametros geometria, int --> define o tamanho da lista acessivel dentro dos arquivos ioside e ioflow sobre o nome geometry[]
# lista com n parametros, double

# numero de materiais fluidos, int --> define a quantidade de materiais fluidos
# lista com n ids dos materiais fluidos, int

# numero de materias s�lidos, int --> define a quantidade de materias s�lidos (s�lido, s�lido livre e s�lido for�ado)
# lista dos n ids dos materias s�lidos, int --> entre apenas com o valor do id da parede, a dummy equivalente ser� o id + 1

# ponto a partir do qual o burner come�ar� a procurar as part�culas, double double double --> deve coincidir com um ponto onde haja part�culas, caso contr�rio o programa acaba sem gerar nenhuma part�cula.