/*
	Programa: st_seq.c  (Versão 2) com base na última versão CUDA
	Descrição: Implementa o Algoritmo para árvore geradora.
	Programadora: Jucele Vasconcellos
	Data: 18/08/2017

	Compilacao:	gcc -o st_seq.exe st_seq.c
	Execucao:	./st_seq.exe in/grafo/grafo1000a seq.out
	
	Entrada de dados: Este programa lê os dados de um grafo no formato
	8
	16
	4 5 0.35
	4 7 0.37
	5 7 0.28
	0 7 0.16
	1 5 0.32
	0 4 0.38
	2 3 0.17
	1 7 0.19
	0 2 0.26
	1 2 0.36
	1 3 0.29
	2 7 0.34
	6 2 0.40
	3 6 0.52
	6 0 0.58
	6 4 0.93

	sendo a primeira linha o número de vértices, a segunda linha o número de arestas
	e as linhas subsequentes as arestas v1 v2 custo
		
	Saída de Dados: Este programa produz um arquivo de saída as arestas que compõem a árvore geradora 
*/

#include <stdio.h> // printf
#include<stdbool.h> // true, false
#include <stdlib.h> //malloc
#include <time.h> //clock

// Grafo Original
typedef struct { 
	unsigned int v, u;
	unsigned int grau_st;
} aresta;


typedef struct { 
	int n, m;
	aresta *arestas;
} grafo;

typedef struct {
	int v1, v2;
} aresta_E;


// Funções e Procedimentos

#define CHECK_ERROR(err) (HandleError(err, __FILE__, __LINE__))

grafo LeGrafo(char *);
void Inicializa_arestasE_C(aresta*, int*, int, int*, aresta_E*, unsigned int*);
void AtualizaC_1(aresta_E*, int*, int*, int); 
void DefineNovosVU(int*, int, char *, unsigned int *); 
void AtualizaC_3(int*, int, char *); 
void MarcarArestas(aresta *, int*, int, int);



// Função Principal
int main (int argc, char** argv){
	grafo G;
	double tempoTotal, tempo1, tempo2;
	double tempo1p, tempo2p;
	int *SolutionEdgeSet;
	int SolutionSize, i, it;
	int num_zerodiff;
	FILE *Arq;
		
	int dimBloco, dimGrid;
	int *menorAresta;
	aresta *arestas_trabalho;
	int n, x, y;
	
	
	// Passo 1: Verificação de parâmetros
	// Passo 2: Leitura dos dados do grafo 
	// Passo 3: Criação do grafo bipartido correspondente às arestas recebidas
	// Passo 4: Encontra a solução
		// Passo 4.1: Escolher arestas que comporão a strut
		// Passo 4.2: Calcular o num_zerodiff e computar novas componenetes conexas
		// Passo 4.3: Compactar o grafo
	
	
	// ==============================================================================
	// Passo 1: Verificação de parâmetros
	// ==============================================================================
	
	//Verificando os parametros
	if(argc < 3 ){
	   printf( "\nParametros incorretos\n Uso: ./cms_seq.exe <ArqEntrada> <ArqSaida> <dimBloco> <S/N> onde:\n" );
	   printf( "\t <ArqEntrada> (obrigatorio) - Nome do arquivo com as informações do grafo (número de vértices, número de arestas e arestas.\n" );
		printf( "\t <ArqSaida> (obrigatorio) - Nome do arquivo de saida.\n" );
		printf( "\t <S ou N> - Mostrar ou não as arestas da MST.\n" );
		return 0;
	} 	
	//Define a dimensão do bloco fixada em 32 threads
	dimBloco = 32;
	
	// ==============================================================================
	// Passo 2: Leitura dos dados do Grafo G
	// ==============================================================================
	tempo1p = (double) clock( ) / CLOCKS_PER_SEC;
	G = LeGrafo(argv[1]);
	n = G.n;
// 	printf("****************************\n");
// 	printf("****** GRAFO ORIGINAL ******\n");
// 	MostraGrafoOriginal(GO);
	
// 	printf("Grafo de entrada lido\n");
	
	//Alocação de variável para armazenar solução
	SolutionEdgeSet = (int *) malloc((G.n-1)*sizeof(int)); 
	SolutionSize = 0;

	tempo2p = (double) clock( ) / CLOCKS_PER_SEC;
// 	printf("Tempo Passo 2: %lf\n", tempo2p - tempo1p);
	
	
	// ==============================================================================
	// Passo 4: Encontra solução
	// ==============================================================================
	tempo1 = (double) clock( ) / CLOCKS_PER_SEC;

	menorAresta = (int *) malloc((G.n)*sizeof(int)); 
	arestas_trabalho = (aresta *) malloc(G.m*sizeof(aresta));
	
	for(i = 0; i < G.m; i++)
		arestas_trabalho[i] = G.arestas[i];
	
	it = 0;
	num_zerodiff = 0;
	while(num_zerodiff != 1)
	{
		// ==============================================================================
		// Passo 4.1: Escolher arestas que comporão a strut
		// ==============================================================================
		tempo1p = (double) clock( ) / CLOCKS_PER_SEC;

		//Inicializa menorAresta com -1
		for(i = 0; i < n; i++)
			menorAresta[i] = -1;
		
		for(i = 0; i < G.m; i++)
		{
			if(arestas_trabalho[i].v != G.n)
			{
				if((menorAresta[arestas_trabalho[i].v] == -1) || (menorAresta[arestas_trabalho[i].v] > i))
					menorAresta[arestas_trabalho[i].v] = i;
				if((menorAresta[arestas_trabalho[i].u] == -1) || (menorAresta[arestas_trabalho[i].u] > i))
					menorAresta[arestas_trabalho[i].u] = i;
			}
		}
		
		for(i = 0; i < n; i++)
			arestas_trabalho[menorAresta[i]].grau_st++;
		
		tempo2p = (double) clock( ) / CLOCKS_PER_SEC;
// 		printf("Tempo Passo 4.1: %lf\n", tempo2p - tempo1p);
		
		// ==============================================================================
		// Passo 4.2: Calcular o num_zerodiff
		// ==============================================================================
		tempo1p = (double) clock( ) / CLOCKS_PER_SEC;

		num_zerodiff = 0;
		
		for(i = 0; i < G.m; i++)
		{
			if(arestas_trabalho[i].grau_st > 0)
			{
				SolutionEdgeSet[SolutionSize] = i;
				SolutionSize++;
				if(arestas_trabalho[i].grau_st == 2)
					num_zerodiff++;
			}
		}
		tempo2p = (double) clock( ) / CLOCKS_PER_SEC;
//  		printf("Tempo Passo 4.2: %lf     it = %d     num_zerodiff = %d     SolutionSize = %d\n", tempo2p - tempo1p, it, num_zerodiff, SolutionSize);
		
		// ==============================================================================
		// Passo 4.3: Compactar o grafo
		// ==============================================================================
		if(num_zerodiff != 1)
		{
			// ==============================================================================
			// Passo 4.3.1: Computar componenetes conexas
			// ==============================================================================
			tempo1p = (double) clock( ) / CLOCKS_PER_SEC;

			//Declara variável tam_arestasE
			unsigned int tam_arestasE;
			tam_arestasE = 0;
			aresta_E *arestasE;
			
			arestasE = (aresta_E *)malloc(n * sizeof(aresta_E));

			//Declara variável C
			int *C;
			
			C = (int *) malloc(n * sizeof(int));

			//Inicializa arestasE e C
			for(i = 0;  i < n; i ++)
			{
				if((arestas_trabalho[menorAresta[i]].grau_st == 1) || (((arestas_trabalho[menorAresta[i]].grau_st == 2)) && (i == arestas_trabalho[menorAresta[i]].v)))
				{
					arestasE[tam_arestasE].v1 = arestas_trabalho[menorAresta[i]].v;
					arestasE[tam_arestasE].v2 = arestas_trabalho[menorAresta[i]].u;
					tam_arestasE++;
				}
				C[i] = i;
			}
			
			int fim, c1, c2;
			do
			{
				fim = 0;
				for(i = 0; i < tam_arestasE; i++)
				{
					x = arestasE[i].v1;
					y = arestasE[i].v2;
					c1 = C[x];
					c2 = C[y];
					if(c1 < c2)
					{
						if(C[y] > c1)
							C[y] = c1;
						fim = 1;
					}
					else if(c2 < c1)
					{
						if(C[x] > c2)
							C[x] = c2;
						fim = 1;
					}
				}
			}while (fim == 1);

			int *marcador;
			marcador = (int *) malloc (n * sizeof(int));
			x = 0;
			for(i = 0; i < n; i++)
			{
				if(C[i] == i)
				{					
					C[i] = x;
					x++;
					marcador[i] = 1;
				}
				else
					marcador[i] = 0;
			}	
			
			for(i = 0; i < n; i++)
				if(marcador[i] == 0)
					C[i] = C[C[i]];

			tempo2p = (double) clock( ) / CLOCKS_PER_SEC;
// 			printf("Tempo Passo 4.3.1: %lf\n", tempo2p - tempo1p);
			
			// ==============================================================================
			// Passo 4.3.2: Marcar arestas
			// ==============================================================================
			tempo1p = (double) clock( ) / CLOCKS_PER_SEC;

			for(i = 0; i < G.m; i++)
			{
				if(arestas_trabalho[i].v != G.n)
				{
					x = C[arestas_trabalho[i].v];
					y = C[arestas_trabalho[i].u];
				
					if( x !=  y )
					{
						arestas_trabalho[i].v = x;
						arestas_trabalho[i].u = y;
					}
					else
					{
						arestas_trabalho[i].v = G.n;
						arestas_trabalho[i].u = G.n;
					}
					arestas_trabalho[i].grau_st = 0;
				}

			}
			tempo2p = (double) clock( ) / CLOCKS_PER_SEC;
// 			printf("Tempo Passo 4.3.2: %lf\n", tempo2p - tempo1p);		
		}

		it++;
		n = num_zerodiff;
	} // fim while(num_zerodiff != 1)
	
	
	tempo2 = (double) clock( ) / CLOCKS_PER_SEC;
	tempoTotal = tempo2 - tempo1;

	printf("%lf\n", tempoTotal);

	Arq = fopen(argv[2], "a");
	fprintf(Arq, "\n*** Arquivo de entrada: %s\n", argv[1]);
	fprintf(Arq, "Tempo Total: %lf\n", tempoTotal); 
	fprintf(Arq, "Número de iterações: %d\n", it);
	fprintf(Arq, "SolutionSize: %d\n", SolutionSize);

  	if((argc == 4) && (argv[3][0] == 'S' || argv[3][0] == 's'))
	{
  		fprintf(Arq, "*** ST formada pelas %d arestas\n", SolutionSize);
  		for(i = 0; i < SolutionSize; i++)
  			fprintf(Arq, "Aresta %d - %d\n", G.arestas[SolutionEdgeSet[i]].v, G.arestas[SolutionEdgeSet[i]].u);
  	}
  	fclose(Arq);

	
	// Liberando variávais alocadas no device
	
	// Liberando variávais alocadas no host
	free(G.arestas);
	free(arestas_trabalho);
	free(menorAresta);
	free(SolutionEdgeSet);
	
	return 0;

}


// ==============================================================================
// Função LeGrafo:  Lê as informações do Grafo de um arquivo e armazena em uma 
//                  estrutura
// ==============================================================================
grafo LeGrafo(char *Arquivo){
	int i, aux;
	grafo G;
	float custo;
   FILE *Arq;
    
   Arq = fopen(Arquivo, "r");

   i = 0;
	fscanf(Arq,"%d",&i);
	G.n = i;
	
	fscanf(Arq,"%d",&i);
	G.m = i;
	
	G.arestas = (aresta *) malloc(G.m*sizeof(aresta)); 
	
	for(i = 0; i < G.m; i++){
		fscanf(Arq,"%hu",&G.arestas[i].u);
		fscanf(Arq,"%hu",&G.arestas[i].v);
		if(G.arestas[i].v > G.arestas[i].u)
		{
			aux = G.arestas[i].v;
			G.arestas[i].v = G.arestas[i].u;
			G.arestas[i].u = aux;
		}
		fscanf(Arq,"%f",&custo);
		G.arestas[i].grau_st = 0;
	}
	
	fclose(Arq);
   return G;
}



