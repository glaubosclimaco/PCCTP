/*
 * grasp.cpp
 *
 * Módulo que implementará o GRASP, e os subalgoritmos usados em cada fase.
 *
 *  Created on: 20/05/2012
 *      Author: Rogerio
 */

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <ctime>
//#include <fstream>
//#include <iostream>
//#include <string>
#include "mineracao.h"
#include "random_provider.h"
#include "Utils.h"
#include "PR.h"
#include "./model/Graph.h"

using namespace std;

#define MAX 501

#define debug 0

int valor_otimo = 0;
int qtd_de_padroes = 0;

int comp_vizinhos(const void* p1, const void* p2);
void preenche_matriz(int m[MAX][MAX], int* rota);
void print_matriz(int m[MAX][MAX], int n);
string intToString(int);
string solutionToString(int*);

//Tipo Estruturado para vértice
//typedef struct vertice {
//	int id; //Precisei..então botei., guarda o número do vertices (0, 1, 2, ... N)
//	int x;
//	int y;
//	char tipo[4];
//	int premio;
//	int *cobre; //Conjunto de vertices W, que este vértice V cobre
//	int *coberto_por; //Conjunto de vertices V, que cobre o este W
//	int qtd_cobre;
//	int qtd_coberto_por;
//	int reduzido;
//} Vertice;
//
////Tipo usado para gerar NP
//typedef struct vizinho {
//	int vertice;
//	float distancia;
//} Vizinho;

//Headers
int in_int_vetor(int a, int *v, int tam);
int copy_int_vetor(int *vetor_from, int *vetor_to, int size);
int antecessor_de(int *s, int vertice);
int vertices_da_solucao(int*& vertices_solucao, int*& solucao);
float custo(Matriz* C, int *solucao);
int get_premio_vertice(int vertice, Vertice *vertices_tsp);
int qtd_vertices_em_s(int *s, int qtd_n);

//GENIUS
void GENI_US(vector<Vertex>, Vertice *vertices, Matriz* C, float alfa,
		int *&rota, int prize, float d, int qtd_n, int qtd_t, int minerou,
		DMPadrao *padroes, int qtd_de_padroes, int*& subProbSize);
int remover_genius(Vertice *vertices, int *&rota, int vertice_v, Matriz *C,
		int qtd_n, int tamanho_da_rota, int *vertices_solucao);
void GENI(Vertice* vertices, int *&rota, Matriz *C, int vertice_v, int p,
		int qtd_n);
int vertices_de_LC(int *LC, int *&vertices_lc, int qtd_n);
void criar_LC_4_genius(Vertice *vertices, int *rota, int *LC, int qtd_n,
		int somente_t);
int criar_LRC_4_genius(Vertice *vertices, int *vertices_lc, int qtd_lc,
		int *rota, int qtd_n, Matriz *C, float alfa);

//Enum para tipo de Movimento
enum movimento {

	SHIFT = 0, //Mudar um vertice aleatorio de posicao dentro da solucao

	SWAP = 1, //Trocar Posicao entre dois vertices da solucao

	OR_OPT = 2, //Transferir N vertices consecutivos para outra posicao dentro da solucao

	TWO_OPT = 3, //Remove duas arestas nao adjacentes e inserir duas novas

	THREE_OPT = 4, //Remove tres areas nao adjacentes e inserir tres novas

	DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA = 5, //Procura um V/T candidato que subtitua 2 já presente da rota

	REMOVER_RE_INSERCAO_MAIS_BARATA = 6, //Remove um vertice V e reinsere usando insercao mais barata

	REMOVER_INSERCAO_MAIS_BARATA = 7, //Remove um V/T e insere um novo vertice de LC em S;

	REMOVER_RE_INSERIR_GENIUS = 8, //Remove um V/T e o re-insere GENI

	REMOVER_INSERIR_GENIUS = 9, //Remove um V/T e insere um novo vertice de LC em S via GENI

	REMOVER_MAIS_BARATA_RE_INSERIR_GENIUS = 10, //Remover o vertice de maior custo e reinserí-lo via GENI

	REMOVER_GENIUS_RE_INSERIR_GENIUS = 11,

	REMOVER_GENIUS_INSERIR_GENIUS = 12,

	REMOVER_GENIUS_RE_INSERCAO_MAIS_BARATA = 13, // Remover Genius Reinserir IB

	REMOVER_GENIUS_INSERCAO_MAIS_BARATA = 14,

};

enum algoritmo_construtivo {

	INSERCAO_MAIS_BARATA = 0,

	GENIUS = 1,

};

std::string nome_movimento(int movimentoEnum) {

	switch (movimentoEnum) {
	case SHIFT:
		return "SHIFT";
		break;
	case SWAP:
		return "SWAP";
		break;
	case OR_OPT:
		return "OR_OPT";
		break;
	case TWO_OPT:
		return "TWO_OPT";
		break;
	case THREE_OPT:
		return "SHIFT";
		break;
	case DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA:
		return "Dupla Remoção Simples - Inserir IB";
		break;
	case REMOVER_RE_INSERCAO_MAIS_BARATA:
		return "Remover Simples - Reinserir IB";
		break;
	case REMOVER_INSERCAO_MAIS_BARATA:
		return "Remover Simples - Inserir IB";
		break;
	case REMOVER_RE_INSERIR_GENIUS:
		return "Remover Simples - Reinserir GENIUS";
		break;
	case REMOVER_INSERIR_GENIUS:
		return "Remover Simples - Inserir GENIUS";
		break;
	case REMOVER_MAIS_BARATA_RE_INSERIR_GENIUS:
		return "Remover IB - Reinserir GENIUS";
		break;
	case REMOVER_GENIUS_RE_INSERIR_GENIUS:
		return "Remover e Reinserir GENIUS";
		break;
	case REMOVER_GENIUS_INSERIR_GENIUS:
		return "Remover e Inserir GENIUS";
		break;
	case REMOVER_GENIUS_RE_INSERCAO_MAIS_BARATA:
		return "Remover GENIUS - Reinserir IB";
		break;
	case REMOVER_GENIUS_INSERCAO_MAIS_BARATA:
		return "Remover GENIUS - Inserir IB";
		break;
	}

	return "";
}

int vertices_w_descobertos(Vertice *vertices, int *rota, int qtd_n,
		int *&vertices_w) {

	int qtd_w_descoberto = 0;
	int coberto = 0;

	for (int i = 0; i < qtd_n; ++i) {
		if (strcmp(vertices[i].tipo, "w") == 0 and vertices[i].reduzido == 0) { //se for W
		//Vertices que cobrem este W
			for (int j = 0; j < vertices[i].qtd_coberto_por; ++j) {
				//if ( in_int_vetor(vertices[i].coberto_por[j], rota, qtd_n) == 1) { //esta coberto
				if (rota[vertices[i].coberto_por[j]] != -1) { //Ha na rota pelo menos um j que cobre i
					coberto = 1;
					break;
				}
			}

			//Se não encontrou nada Rota nenhum dos o cobre
			if (coberto == 0) {
				//printf("Vertice W = %d esta descoberto!\n", vertices[i].id);

				if (qtd_w_descoberto == 0) {
					vertices_w = (int *) malloc(sizeof(int));
				} else {
					vertices_w = (int *) realloc(vertices_w,
							(qtd_w_descoberto + 1) * sizeof(int));
				}
				vertices_w[qtd_w_descoberto] = vertices[i].id;

				qtd_w_descoberto++;
			} else {
				coberto = 0;
			}

		}
	}

	return qtd_w_descoberto;

}

int cobre_algum_w_descoberto(Vertice *vertices, Matriz* C, float d,
		int *s_inicial, int qtd_n, int vertice) {
	/*
	 * Retorna 1 Se o vertice recebido cobre algum W ainda descoberto(ou seja, fora da solucao)
	 */
	if (existe_w_descoberto(vertices, C, d, s_inicial, qtd_n) == 1) {
		return 1;
	} else {
		return 0;
	}
}

void imprimir_solucao_em_arquivo(int *solucao, Matriz *C, Vertice *vertices) {
	int p = 0;

//	FILE* out = fopen("solucao_x.txt", "wt");
//
//	while (solucao[p] == -1) {
//		p++;
//	}
//
//	//p e primeiro vertice pertencente a solucao
//	int antecessor = p;
//	int sucessor = solucao[antecessor];
//	fprintf(out, "%d;%d;%.0f;%d,%s\n", antecessor, sucessor,
//			get(C, antecessor, sucessor),
//			get_premio_vertice(antecessor, vertices),
//			vertices[antecessor].tipo);
//	////printf("%d;%d;%.0f\n",antecessor,sucessor, get(C, antecessor, sucessor));
//	while (sucessor != p) {
//		antecessor = sucessor;
//		sucessor = solucao[sucessor];
//		fprintf(out, "%d;%d;%.0f;%d,%s\n", antecessor, sucessor,
//				get(C, antecessor, sucessor),
//				get_premio_vertice(antecessor, vertices),
//				vertices[antecessor].tipo);
//	}
//
//	fclose(out);

}

void imprimir_solucao(int *solucao, Vertice* vertices) {
	//printf("\n*** Print da Solucao ***\n");
	int p = 0;
	while (solucao[p] == -1) {
		p++;
	}

	//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = solucao[antecessor];
	printf("\n%d  --> %d", antecessor, sucessor);
	int count = 1;
	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = solucao[sucessor];
		printf(" --> %d", sucessor);

		count++;
//		if (count > 100)
//			scanf("%d", &count);

	}

	printf("\n%d vertices\n", count);

}

void printSolVector(vector<edgeSol> sol) {
	printf("\n");
	int python = 1;
	unsigned narestas = 0;
	double totalCost = 0.0;
	for (int i = 0; i < sol.size(); ++i) {
		printf("%d -> %d [%d]: %0.f -> \t", sol[i].o, sol[i].d,
				sol[i].frequencia, sol[i].custo);
		totalCost += sol[i].custo;
		narestas++;
	}
	printf("\n");
	printf("arestas = %d \n", narestas);
	printf("custo: %0.f \n", totalCost);

	if (python) {
		for (int i = 0; i < sol.size(); ++i) {
			printf("(\"%d\", \"%d\", \"1\"),\t", sol[i].o, sol[i].d);

		}
		printf("\n");
	}

}

void atualizar_LC(Vertice *vertices, int *rota, int *LC, int k, int qtd_n) {

	/*
	 * 1) Atualiza a lista de candidatos, retirando os ja inseridos da solucao;
	 * 2) Enquanto houve T fora da solucao, somente eles sao candidatos
	 * */

	//Atualizar a Lista de Candidatos
	LC[k] = 0;

	if (existe_t_fora_da_solucao(vertices, rota, qtd_n) == 1) {

		//LC[k] = 0;

	} else {
		for (int i = 0; i < qtd_n; ++i) {

			//Se for v/t e nao estiver em s_inicial
			//if ( (strcmp(vertices[i].tipo,"v/t") == 0) && ( in_int_vetor(i, rota, qtd_n) == 0 ) ){
			if ((strcmp(vertices[i].tipo, "v/t") == 0) && (rota[i] != -1)) {
				LC[i] = 1;
			}

		}
	}
}

void criar_LC(Vertice *vertices, int a, int b, int *LC, int qtd_n) {
	/* Vetor Binário, onde o indice identifica o vertice e valor é binário,
	 * onde 1 indica vertice candidato
	 * a e b sao os dois vertices iniciais de s_inicial, escolhidos aleatoriamente entre os T
	 */

	/*
	 *
	 * Criar_LC: Criar a lista inicial de candidatos incluindo somente com os Obrigatórios(T).
	 *
	 */

	for (int i = 0; i < qtd_n; ++i) {
		//Se for W nao é candidato
		if ((strcmp(vertices[i].tipo, "t") == 0)) {
			LC[i] = 1;
		} else {
			LC[i] = 0;
		}
	}
	LC[a] = 0;
	LC[b] = 0;

}

float custo_inserir_k(Vertice *vertices, Matriz* C, int v_ant, int k,
		int v_pos) {

	float custo = get(C, v_ant, k) + get(C, k, v_pos) - get(C, v_ant, v_pos);
	/*
	 if ( strcmp(vertices[k].tipo, "v/t") == 0 ){
	 custo *= 1.75;
	 }
	 */

	return custo;
}

void definir_custos_dos_candidatos(Vertice *vertices, Matriz* C, int *rota,
		int *LC, Matriz *C_insert, int *pmin, int *pmax, int qtd_n) {
	/*
	 * C_insert tem dimensao 2 x |LC|, ou seja, para cada vertice candidato N pertencente LC,
	 * é guardado o menor custo, bem como a posicao que gera esse menor custo.
	 */
	//Preencher C_insert com -1, em todas as posicoes.
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < qtd_n; j++) {
			setMatriz(C_insert, i, j, -1);
		}
	}

	//Preencher C_insert
	for (int k = 0; k < qtd_n; k++) {

		if (LC[k] == 1) { //Se candidato

			int p = 0; //auxiliar de controle de posicao P.

			//procurar primeiro vertice constante na solucao
			while (rota[p] == -1) {
				p++;
			}

			//primeira aresta v_ant --> v_pos
			////printf("Valor do P = %d\n",p);
			////printf("S_inicial[p] = %d\n",s_inicial[p]);
			int v_ant = p;
			int v_pos = rota[v_ant];
			//Custo de inserir k, entre v_ant e v_pos
			////printf("Valores %d\t%d\t%d\n", v_ant, v_pos, k);
			////printf("Valores %.2f\t%.2f\t%.2f\n", C[v_ant][k], C[k][v_pos], C[v_ant][v_pos]);
			float k_custo_min = custo_inserir_k(vertices, C, v_ant, k, v_pos);
			//Posicao do menor custo de inserir k na solucao atual(vertice da esquerda)
			int k_pos_custo_min = v_ant;

			////printf("\nk = %d\n", k);
			////printf("Primeira aresta v_ant: %d, v_pos: %d, custo de insercao: %.2f\n", v_ant, v_pos, k_custo_min);
			//Fazer o mesmo para todas as arestas na solucao (s_inicial)
			while (v_pos != p) { // ou seja, quando fechar o ciclo
				v_ant = v_pos;
				v_pos = rota[v_ant];
				//custo na aresta atual da solucao, de v_ant até v_pos
				float k_custo_x = custo_inserir_k(vertices, C, v_ant, k, v_pos);
				if (k_custo_x < k_custo_min) {
					k_custo_min = k_custo_x;
					k_pos_custo_min = v_ant;
				}
			}
			//Guarda melhor posicao de inserir k na solucao, e o custo;
			//Ponderar o custo_min, de tal forma ao vertices T, serem privilegiados.
			setMatriz(C_insert, 0, k, k_custo_min);
			setMatriz(C_insert, 1, k, k_pos_custo_min);
			/*
			 //printf("Menor custo de k:%d == %.2f\n", k, C_insert[0][k]);
			 //printf("Melhor Posicao de k:%d == %d\n", k, (int)C_insert[1][k]);
			 //printf("k = %d finalizado\n",k);
			 */

		}
	}

	//Identificar menor e maior custo dentro os coletado em C_insert
	int c = 0;
	for (int i = 0; i < qtd_n; ++i) {
		if (LC[i] == 1) {
			c = i;
			break;
		}
	}

	float min_geral = get(C_insert, 0, c);
	float max_geral = get(C_insert, 0, c);
	*pmin = c;
	*pmax = c;
	for (int i = 1; i < qtd_n; ++i) {
		if (LC[i] == 1) {
			if (get(C_insert, 0, i) < min_geral) {
				min_geral = get(C_insert, 0, i);
				*pmin = i;
			} else if (get(C_insert, 0, i) > max_geral) {
				max_geral = get(C_insert, 0, i);
				*pmax = i;
			}
		}
	}
}

int definir_LRC_e_get_k_aleatoriamente(int *LC, Matriz *C_insert, int min,
		int max, float alfa, int qtd_n) {
	/*
	 * Definir quais vertices entram em LCR e efetua o sorteio.
	 * Recebe: Matriz de menores custos de insercao, indices dos vertices de maior e menor custo, entre os menores,
	 * fator de aleatorieadade alfa e LC, que indica quem sao os candidatos
	 */

	//Tamanho de lista definir por valor(e não por cardinalidade)
	//C_min + a * (S_max - S_min)
	////printf("INICIO LRC\n");
	float ref = (float) get(C_insert, 0, min)
			+ alfa * ((float) (get(C_insert, 0, max) - get(C_insert, 0, min)));
	////printf("Valor de Ref. para LRC = %.2f\n", ref);

	float custo_k = 0;

	//o tamanho de LRC é variavel, logo tem alocacao dinamica incremental
	int * LRC;
	int size_LRC = 0;

	//TO-REFACTOR
	for (int k = 0; k < qtd_n; k++) {
		if (LC[k] == 1) { //Se candidatos
			custo_k = get(C_insert, 0, k);
			if (custo_k <= ref) {
				if (size_LRC == 0) {
					LRC = (int *) malloc(sizeof(int));
					LRC[size_LRC] = k;
				} else {
					LRC = (int *) realloc(LRC, (size_LRC + 1) * sizeof(int));
					LRC[size_LRC] = k;
				}
				size_LRC++;
				////printf("Vertice: %d entra na LRC\n", k);
			}
		}
	}

	//Sortear k em LRC
	////printf("Tamanho da LRC: %d\n", size_LRC);
	//srand(seed);
	int sorteado = genrand_int32() % (size_LRC);
	////printf("Sorteado posicao LRC: %d\n", sorteado);
	int v_sorteado = LRC[sorteado];

	free(LRC);
	return v_sorteado;
}

void inserir_k_na_solucao(int *&solucao, int k, int pos) {

	int sucessor_pos = solucao[pos];

	solucao[pos] = k;
	solucao[k] = sucessor_pos;
}

void remover_k_da_solucao(int *&rota, int k) {
	int antecessor = antecessor_de(rota, k);
	rota[antecessor] = rota[k];
	rota[k] = -1;
}

void insercao_mais_barata(int qtd_t, int *&rota, Vertice *vertices, int qtd_n,
		Matriz* C, float d, int prize, float alfa) {

	int vertices_t[qtd_t];
	int LC[qtd_n];
	int count_t = 0;
	int count_w = 0;
	//int count_v = 0;
	int count_vt = 0;

	//Matriz com cursos de insercao, Duas dimensoes: 0 - Custo; 1 - Posicao
	Matriz* C_insert = create_matriz(2, qtd_n);

	//Separando vertices T para selecionar 2 aleatoriamente
	for (int i = 0; i < qtd_n; ++i) {
		rota[i] = -1; //prepar s_inicial
		if (strcmp(vertices[i].tipo, "t") == 0) {
			vertices_t[count_t] = i;
			count_t++;
		}
		if (strcmp(vertices[i].tipo, "w") == 0 and vertices[i].reduzido == 0) {
			count_w++;
		}
		if (strcmp(vertices[i].tipo, "v/t") == 0) {
			count_vt++;
		}
	}

	//count_v = count_t + count_vt;
	////printf("\n--> |V| = %d (|T| = %d; |V/T| = %d); |W| = %d\n\n", count_v, count_t, count_vt, count_w);

	//Passo 1: Pegar dois vertices aleatorios de T
	int a = genrand_int32() % qtd_t;
	int b = genrand_int32() % qtd_t;
	while (a == b)
		a = genrand_int32() % qtd_t;
	//sorteio duas posicoes em vertices_t, o valor em cada posicao é o indice do vertice em vertices
	////printf("Posicao a sorteado: %d, Vertice: %d \n", a, vertices_t[a]);
	////printf("Posicao b sorteado: %d, Vertice: %d \n", b, vertices_t[b]);
	//vertices
	//int i = vertices[vertices_t[a]].id;
	int i = vertices_t[a];
	//int j = vertices[vertices_t[b]].id;
	int j = vertices_t[b];
	//Adicionar na solucao i e j na solucao
	rota[i] = j;
	rota[j] = i;
	////printf("Vertices da Solucao Inicial: \n");
	////printf("%d (%s), %d (%s)", i, vertices[i].tipo, j, vertices[j].tipo);
	//Passo 2: Definir LC, excluindo i e j
	//Somente com vertices obrigatorios
	criar_LC(vertices, i, j, LC, qtd_n);
	//**** DEBUG ****
	////printf("SOLUCAO: %d --> %d",i, j);
	/*
	 //printf("2 Vertices T iniciais da solucao= %d, %d\n",i,j);
	 //printf("Demais Vertices de V, ou seja, Candidatos:\n");
	 for (int i = 0; i < qtd_nos; ++i) {
	 if (LC[i] == 1) {
	 //printf("-> %d\n",i);
	 }
	 }

	 //printf("\n\nSolucao inicial 2 vertices:\n");
	 int p = 0;
	 int ant, pos;
	 while (s_inicial[p] == -1)
	 p++;
	 ant = p;
	 pos = s_inicial[ant];
	 //printf("aresta: %d --> %d\n",ant,pos);
	 ant = pos;
	 pos = s_inicial[pos];
	 //printf("aresta: %d --> %d\n",ant,pos);
	 */
	//**** DEBUG ****
	int cont = 0;
	int prize_obtido = prize_acumulado(vertices, rota, qtd_n);
	int k_custo_min = -1;
	int k_custo_max = -1;
	while (existe_w_descoberto(vertices, C, d, rota, qtd_n)
			|| existe_t_fora_da_solucao(vertices, rota, qtd_n)
			|| (prize_obtido < prize)) {
		/* Calcular o custo de inserir todos os vertices
		 * em LC em qualquer posicao entre vertices adjacentes
		 * Função de Custo: Distancia e Premio
		 */
		//C_insert guarda a melhor posicao e o custo de inserir cada k pertencente a LC.
		//Indices do vertices com menor e maior custo em LC
		definir_custos_dos_candidatos(vertices, C, rota, LC, C_insert,
				&k_custo_min, &k_custo_max, qtd_n);
		// **** DEBUG ****
		/*
		 //printf("\nn*** Tabela de menor custo de k pertence a LC:****\n");
		 for (int i = 0; i < qtd_nos; i++) {
		 if (LC[i] == 1)
		 //printf("Vertice: %d ( %s ), Depois do Vertice: %d, Custo:%.2f\n", i, vertices[i].tipo,(int)C_insert[1][i], (float)C_insert[0][i]);
		 }

		 //printf("\n--------MENOR E MAIOR--------\n");
		 //printf("Menor Custo de Insercao--> Vertice: %d, Depois do Vertice: %d, Custo:%.2f\n", k_min, (int)C_insert[1][k_min],(float)C_insert[0][k_min]);
		 //printf("Maior Custo de Insercao--> Vertice: %d, Depois do Vertice: %d, Custo:%.2f\n", k_max, (int)C_insert[1][k_max],(float)C_insert[0][k_max]);

		 //imprimir_solucao(s_inicial);
		 */
		//Construir LCR apartir de LC;
		////printf("***Definindo LRC e Sorteando k ***\n");
		int k = definir_LRC_e_get_k_aleatoriamente(LC, C_insert, k_custo_min,
				k_custo_max, alfa, qtd_n);
		////printf("Vertice k sorteado: %d, posicao = %d\n", k, (int)C_insert[1][k]);
		//Inserir k na s_inicial, na melhor posicao
		inserir_k_na_solucao(rota, k, (int) (get(C_insert, 1, k)));
		////printf(", %d(%s)", k, vertices[k].tipo);
		//atualizar LC
		atualizar_LC(vertices, rota, LC, k, qtd_n);
		cont++;

		//Acumular premio do novo vertice
		prize_obtido += vertices[k].premio;

		//free_matriz(&C_insert, 2, qtd_n);
		/*
		 //printf("\nPRIZE Solicitado: %d\n", prize);
		 //printf("\nPremio acumulado: %d\n", prize_obtido);
		 //printf("Existe W descoberto? %d \n", existe_w_descoberto(vertices, C, d, s_inicial, qtd_n));
		 //printf("Existe T fora? %d \n", existe_t_fora_da_solucao(vertices, s_inicial, qtd_n));
		 //printf("\nPremio acumulado: %d\n", prize_obtido);
		 */
	}

	free_matriz(C_insert);

}

int cobre_algum_vertice_exclusivamente(int *rota, Vertice *vertices, Matriz *C,
		float d, int v, int qtd_n) {

	/*
	 * Dada um solucao S e um vertice V, verificar se V cobre algum vertice do conjunto W, sozinho na solucao;
	 * Estrategia, pegar os Ws que V cobre e verificar se de todos os vertices que o cobre cada W
	 * somente V esta na solucao;
	 * Status = OK - Precisa Refatorar
	 */
	////printf("\nVertice V/T: %d cobre %d vertices de W: ", v, vertices[v].qtd_cobre);
	/*
	 * Se é V/T e nao cobre niguem e o seu premio é menor
	 * que premio excedente entao retorna logo 0 pra remover o vertice da solucao
	 */
	if (vertices[v].qtd_cobre == 0) {
		return 0;
	}

	int qtd_q_cobre_um_dado_w = 0;

	for (int i = 0; i < vertices[v].qtd_cobre; ++i) { //Ver se cada vertce W que V cobre só coberto por ele
		////printf("%d ((", vertices[v].cobre[i]);

		//Todos os vertices que cobrem um dos W que V cobre e que esteja na solucao corrente;
		for (int j = 0; j < vertices[vertices[v].cobre[i]].qtd_coberto_por;
				++j) {
			if ((in_int_vetor(vertices[vertices[v].cobre[i]].coberto_por[j],
					rota, qtd_n) == 1)
					and (vertices[vertices[v].cobre[i]].coberto_por[j] != v)) {
				qtd_q_cobre_um_dado_w++;
				////printf("%d, ", vertices[vertices[v].cobre[i]].coberto_por[j]);
			}
		}
		//printf(")%d), ", qtd_q_cobre_um_dado_w);

		if (qtd_q_cobre_um_dado_w == 0) { //ou seja, ha um w que so e coberto por V
			return 1;
		} else {
			qtd_q_cobre_um_dado_w = 0; //Zera para contar pro proximo W que V cobre
		}

	}

	return 0;

}

void limpar_solucao(int *&rota, Vertice *vertices, int qtd_n, Matriz* C,
		float d, int prize) {

	/*
	 *LIMPAR_SOLUCAO:
	 * Procurar na solucao vertices V/T cujo premio seja menor, ou igual, ao premio_excedente
	 * e que nao cubra (exclusivamente) nenhum W
	 * Remover esses vertices da solucao;
	 * Status = Tá funcionando, precisa refatorar
	 */
	int p = 0;
	while (rota[p] == -1) {
		p++;
	}

	int premio_excedente = (int) prize_acumulado(vertices, rota, qtd_n) - prize;

	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

	int vertice = p;
	int sucessor = rota[vertice];
	//printf("Vertices que podem ser removidos: ");
	if (strcmp(vertices[vertice].tipo, "v/t") == 0
			and vertices[vertice].premio < premio_excedente) {
		if (cobre_algum_vertice_exclusivamente(rota, vertices, C, d, vertice,
				qtd_n) == 0) {

			float custo_v = get(C, antecessor_de(rota, vertice), vertice)
					+ get(C, vertice, rota[vertice])
					- get(C, antecessor_de(rota, vertice), rota[vertice]);

			//remover o vertice
			printf("\nREMOVER Vertice %d da solucao (Economia = %.0f)!!\n",
					vertice, custo_v);

			remover_k_da_solucao(rota, vertice);

			p = sucessor;
			vertice = p;
			sucessor = rota[vertice];
			tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

		}
	}

	int count = 1;
	premio_excedente = (int) prize_acumulado(vertices, rota, qtd_n) - prize;

	while (sucessor != p) {

		if (strcmp(vertices[vertice].tipo, "v/t") == 0
				and vertices[vertice].premio < premio_excedente) {
			if (cobre_algum_vertice_exclusivamente(rota, vertices, C, d,
					vertice, tamanho_da_rota) == 0) {

				float custo_v = get(C, antecessor_de(rota, vertice), vertice)
						+ get(C, vertice, rota[vertice])
						- get(C, antecessor_de(rota, vertice), rota[vertice]);

				//remover o vertice
				printf("\nREMOVER Vertice %d da solucao (Economia = %.0f)!!\n",
						vertice, custo_v);

				remover_k_da_solucao(rota, vertice);

				if (vertice == p)
					p = sucessor;

			}
		}
		vertice = sucessor;
		sucessor = rota[sucessor];
		count++;

		premio_excedente = (int) prize_acumulado(vertices, rota, qtd_n) - prize;
		tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

	}

	////printf("Premio Excedente: %d (%.1f porcento) \n", premio_excedente, (((float)premio_excedente/prize) * 100));
	free(vertices_solucao);

}

void construir_solucao_inicial(vector<Vertex> toFix, Vertice *vertices,
		Matriz* C, float alfa, int *rota, int qtd_t, int prize, float d,
		int qtd_n, int algoritmo, int minerou, DMPadrao *padroes,
		int*& subProbSize) {
	//Algoritmo de Construcao Inicial: Insercao Mais Barata (Matheus..)

	switch (algoritmo) {
	case INSERCAO_MAIS_BARATA:
		insercao_mais_barata(qtd_t, rota, vertices, qtd_n, C, d, prize, alfa);
		break;
	case GENIUS:
		GENI_US(toFix, vertices, C, alfa, rota, prize, d, qtd_n, qtd_t, minerou,
				padroes, qtd_de_padroes, subProbSize);
		break;
	default:
		break;
	}

	//Fazer um clean na solucao recem criada
	//imprimir_solucao(rota, vertices);
	//printf("\nPRIZE Solicitado: %d\n", prize);

	//printf("\nPRIZE Coletado: %d\n", (int)prize_acumulado(vertices, rota, qtd_n));

	//printf("\nCusto da Solucao INICIAL: %.1f\n", custo(C, rota));

	//limpar_solucao(rota,vertices, qtd_n, C, d, prize);

	//imprimir_solucao(rota, vertices);

	//printf("\nPRIZE Solicitado: %d\n", prize);

	//printf("\nPRIZE Coletado: %d\n", (int)prize_acumulado(vertices, rota, qtd_n));

	//printf("\nCusto da Solucao INICIAL: %.3f\n", custo(C,rota));

}

void inverter_sentido_rota(int*& rota, int vertice_b, int vertice_c) {
	//Inverter sentido da rota entre B e C
	int v = vertice_b;
	int suc = rota[v]; //sucessor de v
	int suc_suc = rota[suc]; //sucessor do sucessor de v
	rota[suc] = v;
	//Verbose //printf("%d --> %d \n", suc, v);
	while (suc != vertice_c) {
		v = suc;
		suc = suc_suc;
		suc_suc = rota[suc];
		rota[suc] = v;
		////printf("%d --> %d \n", suc, v);
	}
}

int obter_lc_dado_uma_solucao_viavel(Vertice *vertices, int *vertices_solucao,
		int tamanho_rota, int LC[], int qtd_n) {
	int tamanho_lc = 0;
	for (int i = 0; i < qtd_n; ++i) {
		if (strcmp(vertices[i].tipo, "v/t") == 0) {
			LC[i] = 1;
			tamanho_lc++;
			//printf("*%d, ", i);
		} else {
			LC[i] = 0;
		}
	}

	for (int i = 0; i < tamanho_rota; ++i) {
		if (strcmp(vertices[vertices_solucao[i]].tipo, "v/t") == 0) {
			LC[vertices_solucao[i]] = 0;
			//printf("#%d, ", vertices_solucao[i]);
			tamanho_lc--;
		}
	}

	return tamanho_lc;
}

int obter_vertices_vt_habilitados(Vertice* vertices, Matriz* C, int *rota,
		int tamanho_rota, int qtd_n, int *vertices_vt, int *vertices_solucao,
		int qtd_vt, int prize, int *&vertices_habilitados) {
	/* Dado um conjunto de V/T e a rota ao qual pertencem, verificar em LC quais estao
	 * habilitados a substitui-lo: Critério Cobertura e Premio
	 */

	//1) Montar LC; 2) Pegar premio requerido; e 3) Pegar W descobertos com saida dos V/T
	//Criar uma LC com todos os Opcionais que esteja fora da rota
	int tamanho_lc = 0;
	int qtd_habilitados = 0;
	int LC[qtd_n];
	int LC_Habilitados[qtd_n]; //Vertices (V/T) que realmente sao candidatos

	tamanho_lc = obter_lc_dado_uma_solucao_viavel(vertices, vertices_solucao,
			tamanho_rota, LC, qtd_n);

	if (tamanho_lc < qtd_vt)
		return 0;

	qtd_habilitados = tamanho_lc;
	copy_int_vetor(LC, LC_Habilitados, qtd_n);

	//printf("Premio Antes de Remover: %d\n", prize_acumulado(vertices, rota, qtd_n));

	//Remover todos os VT recebidos: para saber quais W ficaram descobertos e PRIZE ainda requerido
	for (int i = 0; i < qtd_vt; ++i) {
		remover_k_da_solucao(rota, vertices_vt[i]);
	}

	//PRIZE requerido com saida dos VT
	int premio_requerido = prize - prize_acumulado(vertices, rota, qtd_n);
	//printf("Premio requerido: %d\n", premio_requerido);

	//Pegar W descobertos
	int *vertices_w;
	int qtd_wd = vertices_w_descobertos(vertices, rota, qtd_n, vertices_w);

	for (int i = 0; i < qtd_n; ++i) {

		if (LC[i] == 0)
			continue;
		//printf("Vertices %d e candidato.", i);

		if (vertices[i].premio < premio_requerido) { //Verificar se o premio  é suficiente
			//printf(" Porém seu premio é insuficiente.");
			LC_Habilitados[i] = 0;
			qtd_habilitados--;
		} else if (qtd_wd > 0) { //Verificar se I cobre os W que ficaram descobertos com a saida de V
			for (int j = 0; j < qtd_wd; ++j) { //todos os W descobertos devem pertencer ao conjunto dos vertices que i cobre
				if (in_int_vetor(vertices_w[j], vertices[i].cobre,
						vertices[i].qtd_cobre) == 0) {
					//printf(" Mas ele nao cobre o vertice w = %d que esta descoberto com saida dos V/T que sairam", vertices_w[j]);
					LC_Habilitados[i] = 0;
					qtd_habilitados--;
					break;
				}
			}
		}
		//printf("\n");

		//if (LC_Habilitados[i] == 1)
		//printf("Vertice Habilitado = %d (P = %d)\n", i, vertices[i].premio);

	}

	if (qtd_habilitados > 0)
		vertices_habilitados = (int *) malloc(qtd_habilitados * sizeof(int));

	int index = 0;

	//Adicionar os habilitados no conjunto de retorno
	//if (qtd_habilitados >= qtd_vt) // Se encontrou a qtd necessário de substitutos
	for (int i = 0; i < qtd_n; ++i) {
		if (LC_Habilitados[i] == 1) {
			vertices_habilitados[index] = i;
			index++;
		}
	}

	if (qtd_wd > 0)
		free(vertices_w);

	//Verbose
	//printf("Vertices Retirados (qtd=%d): ", qtd_vt);
	//for (int i = 0; i < qtd_vt; ++i) {
	//printf(" %d, ", vertices_vt[i]);
	//}
	//printf("\n");

	//printf("Vertices Habilitados (qtd=%d): ", qtd_habilitados);
	//for (int i = 0; i < qtd_habilitados; ++i) {
	//printf(" %d, ", vertices_habilitados[i]);
	//}
	//printf("\n");

	return qtd_habilitados;
}

void movimento_shift(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &s_vizinho_encontrado, Matriz* C, Vertice* vertices) {

	/*SHIFT:
	 * 1) Sorteia um vertice na solucao e o reinsere na melhor posicao, se houver.
	 */

	//TO_FIX: Shift já tá retornando um vizinho melhor..e nao um vizinho aleatorio
	int melhor_vizinho_encontrado = 0;

	//Sortear um vertice pertecente na solucao e reinseri-lo em um posicao de que acrescente menor custo
	int vertice_v = -1;
	int count = 0;

	do {
		count++;

		vertice_v = vertices_solucao[genrand_int32() % tamanho_rota];

		int antecessor_de_v = antecessor_de(s_vizinho_encontrado, vertice_v);
		int sucessor_de_v = s_vizinho_encontrado[vertice_v];

		float custo_v_na_solucao = get(C, antecessor_de_v, vertice_v)
				+ get(C, vertice_v, sucessor_de_v)
				- get(C, antecessor_de_v, sucessor_de_v);
		//Nova posicao candidata de V
		int vertice_anterior = sucessor_de_v;
		int vertice_posterior = s_vizinho_encontrado[vertice_anterior];

		do {

			//Custo k de inserir v, entre vertice_anterior e vertice_posterior
			float k_custo = get(C, vertice_anterior, vertice_v)
					+ get(C, vertice_v, vertice_posterior)
					- get(C, vertice_anterior, vertice_posterior);

			if (k_custo < custo_v_na_solucao) {

				//verbose
				/*//imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo da Solucao: %.0f\n\n",custo(vertices, C, s_vizinho_encontrado));
				 */

				//Mudar V de posicao
				s_vizinho_encontrado[antecessor_de_v] = sucessor_de_v;
				s_vizinho_encontrado[vertice_anterior] = vertice_v;
				s_vizinho_encontrado[vertice_v] = vertice_posterior;
				melhor_vizinho_encontrado = 1;
				//Fim-Tirar V da solucao

				//Verbose
				/*//printf("SHIFT do Vertice V: %d\n", vertice_v);
				 //printf("Posicao atual entre os vertices: %d, %d\n", antecessor_de_v, sucessor_de_v);
				 //printf("Posicao Melhor, após: %d\n", vertice_anterior);
				 //printf("Custo atual V na solucao: %.1f \n", custo_v_na_solucao);
				 //printf("Novo de Custo de V na solucao: %.1f \n", k_custo);
				 //imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo da Solucao: %.0f\n\n",custo(vertices, C, s_vizinho_encontrado));
				 */

			} else { //senao avanca na solucao

				vertice_anterior = vertice_posterior;
				vertice_posterior = s_vizinho_encontrado[vertice_posterior];

			}

		} while (melhor_vizinho_encontrado == 0
				and vertice_anterior != vertice_v);

	} while (count < tamanho_rota and melhor_vizinho_encontrado == 0);

}

void movimento_shift_v2(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &s_vizinho_encontrado, Matriz* C, Vertice* vertices) {

	/* SHIFT: Versão na qual tentar o SHIFT apenas em um percetual do vertices mais próximos e não em todos.
	 *      1) Se tamanho_rota <= 20 usar versão atual do shift
	 */

	if (tamanho_rota <= 20) {
		movimento_shift(vertices_solucao, tamanho_rota, qtd_n,
				s_vizinho_encontrado, C, vertices);
	} else {
		//Usar p vértices mais próximos
		//printf("\nSolução mais de 20 vértices usar V2 do SHIFT! \n");

		int p = 20;
		int vertice_v;
		int vertice_k;
		Vizinho vizinhos[tamanho_rota];

		for (int i = 0; i < tamanho_rota; ++i) {

			vertice_v = vertices_solucao[i];

			int antecessor_de_v = antecessor_de(s_vizinho_encontrado,
					vertice_v);
			int sucessor_de_v = s_vizinho_encontrado[vertice_v];

			float custo_v_na_solucao = get(C, antecessor_de_v, vertice_v)
					+ get(C, vertice_v, sucessor_de_v)
					- get(C, antecessor_de_v, sucessor_de_v);

			for (int i = 0; i < tamanho_rota; i++) {
				vizinhos[i].vertice = vertices_solucao[i];
				vizinhos[i].distancia = get(C, vertice_v, vertices_solucao[i]);
			}

			//Ordenar as distancias
			qsort(vizinhos, tamanho_rota, sizeof(Vizinho), comp_vizinhos);

			for (int k = 1; k <= p; k++) {

				vertice_k = vizinhos[k].vertice;
				int sucessor_k = s_vizinho_encontrado[vertice_k];

				float k_custo = get(C, vertice_k, vertice_v)
						+ get(C, vertice_v, sucessor_k)
						- get(C, vertice_k, sucessor_k);

				if (k_custo < custo_v_na_solucao) {

					//Mudar V de posicao
					s_vizinho_encontrado[antecessor_de_v] = sucessor_de_v;
					s_vizinho_encontrado[vertice_k] = vertice_v;
					s_vizinho_encontrado[vertice_v] = sucessor_k;

					return;
				}

			}

		}

	}
}

void movimento_swap(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &s_vizinho_encontrado, Matriz* C, Vertice* vertices) {

	/* SWAP: Sortear vertice da solucao e percorre a solucao inteiro a procura de vertice
	 * para trocar de posicao.
	 * TODO: Mudar para sortear so o primeiro, e sair varrendo a solucao a procura de uma troca boa
	 */

	int vertice_a = -1;
	int vertice_b = -1;
	int antecessor_a = -1;
	int antecessor_b = -1;
	int sucessor_a = -1;
	int sucessor_b = -1;

	vertice_a = vertices_solucao[genrand_int32() % tamanho_rota];
	vertice_b = s_vizinho_encontrado[vertice_a];

	int count = 0;
	int melhor_vizinho_encontrado = 0; //indica se jah encontrou um vizinho melhor.

	while (count < tamanho_rota and melhor_vizinho_encontrado == 0) {
		count++;
		antecessor_a = antecessor_de(s_vizinho_encontrado, vertice_a);
		sucessor_a = s_vizinho_encontrado[vertice_a];
		while (vertice_b != vertice_a) {

			sucessor_b = s_vizinho_encontrado[vertice_b];
			antecessor_b = antecessor_de(s_vizinho_encontrado, vertice_b);

			int custo_atual = get(C, antecessor_a, vertice_a)
					+ get(C, vertice_a, sucessor_a)
					+ get(C, antecessor_b, vertice_b)
					+ get(C, vertice_b, sucessor_b);
			int novo_custo = get(C, antecessor_a, vertice_b)
					+ get(C, vertice_b, sucessor_a)
					+ get(C, antecessor_b, vertice_a)
					+ get(C, vertice_a, sucessor_b);

			if (novo_custo < custo_atual) {
				//FAZER: SWAP
				//Verbose
				/*//printf("********* SWAP *******\n");
				 //imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo antes SWAP: %f\n",custo(vertices, C, s_vizinho_encontrado));
				 */
				//Vertices adjancentes
				if (sucessor_a == vertice_b) { //b é sucessor de a
					s_vizinho_encontrado[antecessor_a] = vertice_b;
					s_vizinho_encontrado[vertice_b] = vertice_a;
					s_vizinho_encontrado[vertice_a] = sucessor_b;
				} else if (sucessor_b == vertice_a) { //a é sucessor de b
					s_vizinho_encontrado[antecessor_b] = vertice_a;
					s_vizinho_encontrado[vertice_a] = vertice_b;
					s_vizinho_encontrado[vertice_b] = sucessor_a;
				} else { //Vertices nao ajdacentes
					//Trocar Posicao dos vertices
					s_vizinho_encontrado[antecessor_a] = vertice_b;
					s_vizinho_encontrado[vertice_b] = sucessor_a;
					s_vizinho_encontrado[antecessor_b] = vertice_a;
					s_vizinho_encontrado[vertice_a] = sucessor_b;
				}

				//Verbose
				/*
				 //printf("Vertice A = %d, antecessor: %d e sucessor: %d\n",vertice_a, antecessor_a, sucessor_a);
				 //printf("Vertice B = %d, antecessor: %d e sucessor: %d\n",vertice_b, antecessor_b, sucessor_b);

				 //imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo apos SWAP: %f\n",custo(vertices, C, s_vizinho_encontrado));
				 */
				//encerrar
				melhor_vizinho_encontrado = 1;
				break;

			} else {
				vertice_b = s_vizinho_encontrado[vertice_b];
			}

		}

		//pula pro proximo vertice
		vertice_a = s_vizinho_encontrado[vertice_a];

	}

}

void movimento_or_opt(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &s_vizinho_encontrado, Matriz* C, Vertice* vertices) {

	//Realocar um sequecia de vértices na rota..

	int qtd_vertices_a_mover = 2;

	/* Vertices que poderao ser sorteados como nova posicao inicial da sequencia
	 * Sao os vertices a serem movidos mais o vertice anterior
	 */
	int vertices_excluidos[qtd_vertices_a_mover + 1];

	//sortear um vertice a partir do qual X vertices serao realocados: Primeiro Vertice
	int vertice_inicial = vertices_solucao[genrand_int32() % tamanho_rota];

	vertices_excluidos[0] = vertice_inicial;

	//Pegar ultimo vertices da sequencia a ser movida
	int vertice_final = vertice_inicial;

	for (int i = 0; i < qtd_vertices_a_mover - 1; ++i) {
		vertice_final = s_vizinho_encontrado[vertice_final];
		vertices_excluidos[i + 1] = vertice_final;
	}

	int vertice_anterior = antecessor_de(s_vizinho_encontrado, vertice_inicial);
	int vertice_posterior = s_vizinho_encontrado[vertice_final];

	int novo_vertice_anterior = -1;
	int count = 1;

	while (count == tamanho_rota - 4) {

		//sortear nova posicao para sequencia de vertices
		novo_vertice_anterior =
				vertices_solucao[genrand_int32() % tamanho_rota];
		while (in_int_vetor(novo_vertice_anterior, vertices_excluidos,
				(qtd_vertices_a_mover + 1))) {
			novo_vertice_anterior = vertices_solucao[genrand_int32()
					% tamanho_rota];
		}
		int novo_vertice_posterior = s_vizinho_encontrado[novo_vertice_anterior];

		float custo_liberado = get(C, vertice_anterior, vertice_inicial)
				+ get(C, vertice_final, vertice_posterior)
				+ get(C, novo_vertice_anterior, novo_vertice_posterior);
		float custo_adicionado = get(C, vertice_anterior, vertice_posterior)
				+ get(C, novo_vertice_anterior, vertice_inicial)
				+ get(C, vertice_final, novo_vertice_posterior);
		//se o custo da movimentacao for menor faz o OR-OPT
		if (custo_adicionado < custo_liberado) {

			/*Verbose
			 //imprimir_solucao(s_vizinho_encontrado);
			 //printf("Custo: %f\n", custo(vertices, C, s_vizinho_encontrado));
			 */

			//Retinar o vertices que serao movidos da solucao, ou seja, reconectar a solucao sem eles
			s_vizinho_encontrado[vertice_anterior] = vertice_posterior;

			//Reinserir os vertices na nova posicao
			s_vizinho_encontrado[novo_vertice_anterior] = vertice_inicial;
			s_vizinho_encontrado[vertice_final] = novo_vertice_posterior;

			/*verbose
			 //printf("Primeiro Vertice: %d\n", vertice_inicial );
			 //printf("Ultimo vertice: %d \n", vertice_final);
			 //printf("Vertice Anterior: %d\n", vertice_anterior);
			 //printf("Vertice Posterior: %d\n", vertice_posterior);
			 //printf("Novo Vertice Anterior: %d\n", novo_vertice_anterior);
			 //imprimir_solucao(s_vizinho_encontrado);
			 //printf("Custo: %f\n", custo(vertices, C, s_vizinho_encontrado));
			 */
			//encerra-exploraracao
			break;

		}

		count++;

	}

}

void movimento_two_opt(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &rota, Matriz* C, Vertice* vertices) {
	//Remover duas arestas nao-adjacentes e inserir duas novas
	/* Pegar dois vertices (A e C), nao-adjancentes e pegar seus posteriores (B e D), respec.
	 * Depois ligar A ao D e B ao C;
	 */
	int vertice_a = -1;
	int vertice_b;
	int vertice_c;
	int vertice_d;
	int antecessor_do_vertice_a;
	int encontrou_vizinho_melhor = 0;
	//int vertice_sorteado = -1;

	//vertice_a = vertices_solucao[genrand_int32() % tamanho_rota];

	//Pára o FOR ao encontrar um vizinho melhor
	for (int cont = 0; cont < tamanho_rota; ++cont) {

		/*vertice_sorteado = vertices_solucao[genrand_int32() % tamanho_rota];
		 while (vertice_sorteado == vertice_a)
		 vertice_sorteado = vertices_solucao[genrand_int32() % tamanho_rota];
		 vertice_a = vertice_sorteado;*/
		vertice_a = vertices_solucao[genrand_int32() % tamanho_rota];

		vertice_b = rota[vertice_a];

		//Pegar vertice anterior ao vertice A
		antecessor_do_vertice_a = antecessor_de(rota, vertice_a);

		//Procurar o primeiro melhor vizinho 2_opt para vertice_a
		/* Sair varrendo sequencialmente a solucao realizando movimentos 2_opt
		 * a partir do vertice posterior a B até chegar no anterior a A
		 * Para ao achar uma melhora*/

		vertice_c = rota[vertice_b];
		vertice_d = rota[vertice_c];

		//Ate achar vizinho melhor ou esgotar as possibilidades para o vertice A sorteado
		while (vertice_c != antecessor_do_vertice_a
				and encontrou_vizinho_melhor == 0) {

			if (get(C, vertice_a, vertice_b) + get(C, vertice_c, vertice_d)
					> get(C, vertice_a, vertice_c)
							+ get(C, vertice_b, vertice_d)) {

				//Verbose
				/*//imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo: %f\n", custo(vertices, C, s_vizinho_encontrado));

				 //printf("Vertice a, b, c e d: = %d, %d, %d, %d\n",vertice_a, vertice_b, vertice_c, vertice_d );
				 */
				//Inverter sentido da rota entre B e C
				inverter_sentido_rota(rota, vertice_b, vertice_c);

				//Ligar novas arestas
				rota[vertice_a] = vertice_c;
				rota[vertice_b] = vertice_d;

				//Marcar para parar o processo
				encontrou_vizinho_melhor = 1;

				/*Verbose
				 //imprimir_solucao(s_vizinho_encontrado);
				 //printf("Custo: %f\n", custo(vertices, C, s_vizinho_encontrado));
				 */

			} else {
				//e avanca para proxima aresta
				vertice_c = rota[vertice_c];
				vertice_d = rota[vertice_c];
			}

		}

		if (encontrou_vizinho_melhor == 1) {
			break;
		}

		//vertice_a = rota[vertice_a];

	}

}

int movimento_three_opt(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &rota, Matriz* C, Vertice* vertices) {

	//Remover tres arestas nao-adjacentes e inserir tres novas (ha 4 possibilidades)

	//imprimir_solucao(rota, vertices);

	int vertice_a = vertices_solucao[genrand_int32() % tamanho_rota];
	int vertice_b = rota[vertice_a];
	int vertice_c = -1;
	int vertice_d = -1;
	int vertice_e = -1;
	int vertice_f = -1;
	int sucessor_vertice_b = -1;
	int sucessor_vertice_d = -1;
	float custo_atual_abcdef = -1;
	float novo_custo = -1;

	//Para o FOR ao encontrar um vizinho melhor
	for (int i = 0; i < tamanho_rota - 5; ++i) {

		//vertice_b = rota[vertice_a];
		sucessor_vertice_b = rota[vertice_b];

		//Primeira possibilidade de aresta (C, D)
		vertice_c = sucessor_vertice_b;
		vertice_d = rota[vertice_c];

		//Qtd viavel de arestas (C, D) , tamanho-rota - 5
		int qtd_arestas_cd = tamanho_rota - (5 + i);
		////printf("Qtd arestas CD = %d\n", qtd_arestas_cd);
		for (int j = 0; j < qtd_arestas_cd; ++j) {
			////printf("Aresta CD n. %d \n", j);
			sucessor_vertice_d = rota[vertice_d];
			//Primeira possibilidade de aresta (E, F)
			vertice_e = sucessor_vertice_d;
			vertice_f = rota[vertice_e];

			//Qtd viavel de areasta (E, F)
			int qtd_arestas_ef = qtd_arestas_cd - j;
			////printf("qtd arestas EF = %d\n", qtd_arestas_ef);
			for (int k = 0; k < qtd_arestas_ef; ++k) {
				////printf("Aresta EF n. %d \n", k);
				//Custo de (A, B) + (C, D) + (E, F);
				custo_atual_abcdef = get(C, vertice_a, vertice_b)
						+ get(C, vertice_c, vertice_d)
						+ get(C, vertice_e, vertice_f);

				//printf("3 - OPT \n");
				//printf("Aresta A, B = (%d --> %d) \n", vertice_a, vertice_b);
				//printf("Aresta C, D = (%d --> %d) \n", vertice_c, vertice_d);
				//printf("Aresta E, F = (%d --> %d) \n", vertice_e, vertice_f);
				//printf("Custo Atual: %1.f \n", custo_atual_abcdef);

				//Verificar as 4 formas diferentes de movimentos 3-opt
				for (int opcao = 0; opcao < 4; ++opcao) {
					switch (opcao) {
					case 0: //add (A, C), (B, E) e (D, F) e inverter (B ... C) e (D ... C)

						novo_custo = get(C, vertice_a, vertice_c)
								+ get(C, vertice_b, vertice_e)
								+ get(C, vertice_d, vertice_f);

						if (novo_custo < custo_atual_abcdef) {

							//mudar rota;
							//printf( "3-Opt custo melhor ! - Tipo 1\n");
							//printf("Novo Custo = %.1f", novo_custo);

							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							inverter_sentido_rota(rota, vertice_b, vertice_c);
							inverter_sentido_rota(rota, vertice_d, vertice_e);
							rota[vertice_a] = vertice_c;
							rota[vertice_b] = vertice_e;
							rota[vertice_d] = vertice_f;

							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							return 0;

						}

						break;

					case 1: //add (A, D), (E, B) e (C, F)

						novo_custo = get(C, vertice_a, vertice_d)
								+ get(C, vertice_e, vertice_b)
								+ get(C, vertice_c, vertice_f);

						if (novo_custo < custo_atual_abcdef) {

							//mudar rota;
							//printf( "3-Opt custo melhor ! - Tipo 2\n");
							//printf("Novo Custo = %.1f", novo_custo);
							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							rota[vertice_a] = vertice_d;
							rota[vertice_e] = vertice_b;
							rota[vertice_c] = vertice_f;
							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							return 0;

						}

						break;

					case 2: // add (A, E), (D, B) e (C, F) e inverter (D ... E)

						novo_custo = get(C, vertice_a, vertice_e)
								+ get(C, vertice_d, vertice_b)
								+ get(C, vertice_c, vertice_f);

						if (novo_custo < custo_atual_abcdef) {

							//mudar rota;
							//printf( "3-Opt custo melhor ! - Tipo 3\n");
							//printf("Novo Custo = %.1f", novo_custo);
							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							inverter_sentido_rota(rota, vertice_d, vertice_e);
							rota[vertice_a] = vertice_e;
							rota[vertice_d] = vertice_b;
							rota[vertice_c] = vertice_f;

							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							return 0;

						}

						break;

					case 3: // add (B, F), (A, D) e (E, C), inverter (B ... C)

						novo_custo = get(C, vertice_b, vertice_f)
								+ get(C, vertice_a, vertice_d)
								+ get(C, vertice_e, vertice_c);

						if (novo_custo < custo_atual_abcdef) {

							//mudar rota;
							//printf( "3-Opt custo melhor ! - Tipo 4\n");
							//printf("Novo Custo = %.1f", novo_custo);

							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							inverter_sentido_rota(rota, vertice_b, vertice_c);
							rota[vertice_b] = vertice_f;
							rota[vertice_a] = vertice_d;
							rota[vertice_e] = vertice_c;

							//imprimir_solucao(rota, vertices);
							//printf("Custo Rota: %.0f", custo(C, rota));

							return 0;

						}

						break;

					}
				}
				//fim opcoes de troca de arestas

				//Avanca uma aresta (E, F)
				vertice_e = rota[vertice_e];
				vertice_f = rota[vertice_e];

			} //for (E, F)

			//avanca uma aresta (C, D)
			vertice_c = rota[vertice_c];
			vertice_d = rota[vertice_c];

		} //for (C, D)

		//Avanca proxima areasta (A, B)
		vertice_a = rota[vertice_a];
		vertice_b = rota[vertice_a];

	} //for (A, B)

	return 0;

}

int movimento_double_brigde(int* vertices_solucao, int tamanho_rota, int qtd_n,
		int* &rota, Matriz* C, Vertice* vertices, int prize) {

	//Pegar 4 arestas, nao-adjacentes, e realizar a double-brigde

	//imprimir_solucao(rota, vertices);

	int vertice_a = vertices_solucao[genrand_int32() % tamanho_rota];
	int vertice_b = rota[vertice_a];
	int vertice_c = -1;
	int vertice_d = -1;
	int vertice_e = -1;
	int vertice_f = -1;
	int vertice_g = -1;
	int vertice_h = -1;
	int sucessor_vertice_b = -1;
	int sucessor_vertice_d = -1;
	int sucessor_vertice_f = -1;
	float custo_atual_abcdefgh = -1;
	float novo_custo = -1;
	int count = 0;
	//int *vertices_habilitados;
	//int qtd_habilitados = 0;

	//Para o FOR ao encontrar um vizinho melhor
	for (int i = 0; i < tamanho_rota - 7; ++i) {

		//vertice_b = rota[vertice_a];
		sucessor_vertice_b = rota[vertice_b];

		//Primeira possibilidade de aresta (C, D)
		vertice_c = sucessor_vertice_b;
		vertice_d = rota[vertice_c];

		//Qtd viavel de arestas (C, D) , tamanho-rota - 7
		int qtd_arestas_cd = tamanho_rota - (7 + i);
		////printf("Qtd arestas CD = %d\n", qtd_arestas_cd);
		for (int j = 0; j < qtd_arestas_cd; ++j) {
			////printf("Aresta CD n. %d \n", j);
			sucessor_vertice_d = rota[vertice_d];
			//Primeira possibilidade de aresta (E, F)
			vertice_e = sucessor_vertice_d;
			vertice_f = rota[vertice_e];

			//Qtd viavel de areasta (E, F)
			int qtd_arestas_ef = qtd_arestas_cd - j;
			////printf("qtd arestas EF = %d\n", qtd_arestas_ef);
			for (int k = 0; k < qtd_arestas_ef; ++k) {
				////printf("Aresta EF n. %d \n", k);
				sucessor_vertice_f = rota[vertice_f];

				//Primeira possibilidade de aresta (G, H)
				vertice_g = sucessor_vertice_f;
				vertice_h = rota[vertice_g];

				int qtd_aretas_gh = qtd_arestas_ef - k;
				for (int l = 0; l < qtd_aretas_gh; ++l) {

					////printf("A(%d) - B(%d) --> C(%d) - D(%d) --> E(%d) - F(%d) --> G(%d) - H(%d) \n", vertice_a, vertice_b, vertice_c, vertice_d, vertice_e, vertice_f, vertice_g, vertice_h);
					count++;
					/* Double-Brigde Interna: Faz a reorganização da solução somente com os vertices
					 * que já estão na solução.
					 */
					custo_atual_abcdefgh = get(C, vertice_a, vertice_b)
							+ get(C, vertice_c, vertice_d)
							+ get(C, vertice_e, vertice_f)
							+ get(C, vertice_g, vertice_h);
					novo_custo = get(C, vertice_a, vertice_f)
							+ get(C, vertice_g, vertice_d)
							+ get(C, vertice_e, vertice_b)
							+ get(C, vertice_c, vertice_h);

					/* Double-Brigde Externa: Faz a reorganização da solução substituindo 2 ou 4 (...) vertices V/T
					 * 1) Deve-se obter (AB), (CD), (EF) ou (FG) vertices de V/T
					 * 2) Criar uma funcao que passado 1 conjunto de vertices V/T, retornar todos os candidatos
					 *     capazes de substituí-los na solucao, seja pelo Premio seja pela Cobertura;
					 */
					if ((strcmp(vertices[vertice_a].tipo, "v/t") == 0)
							and (strcmp(vertices[vertice_b].tipo, "v/t"))) { //AB é V/T
						//Pegar V/T habilitados a substituir AB
						//Tentar double-brige externa
						//int vertices_vt[2];
						//vertices_vt[0] = vertice_a;
						//vertices_vt[1] = vertice_b;

						//imprimir_solucao(rota, vertices);
						//qtd_habilitados = obter_vertices_vt_habilitados(vertices,C, rota, tamanho_rota, qtd_n,vertices_vt, vertices_solucao, 2, prize, vertices_habilitados);
						//if (qtd_habilitados > 0)
						//printf("TODO!");

						//scanf("%d", &qtd_n);

					}
					if ((strcmp(vertices[vertice_c].tipo, "v/t") == 0)
							and (strcmp(vertices[vertice_d].tipo, "v/t"))) { //CD é V/T
						//Pegar V/T habilitados a substituir CD
						//Tentar double-brige externa
					}
					if (((strcmp(vertices[vertice_e].tipo, "v/t") == 0)
							and (strcmp(vertices[vertice_f].tipo, "v/t")))) { //EF é V/T
						//Pegar V/T habilitados a substituir EF
						//Tentar double-brige externa
					}
					if (((strcmp(vertices[vertice_g].tipo, "v/t") == 0)
							and (strcmp(vertices[vertice_g].tipo, "v/t")))) {//GH é V/T
						//Pegar V/T habilitados a substituir GH
						//Tentar double-brige externa
					}

					if (novo_custo < custo_atual_abcdefgh) {
						//scanf("%d", &qtd_n);
						//printf("A(%d) - B(%d) --> C(%d) - D(%d) --> E(%d) - F(%d) --> G(%d) - H(%d) \n", vertice_a, vertice_b, vertice_c, vertice_d, vertice_e, vertice_f, vertice_g, vertice_h);
						rota[vertice_a] = vertice_f;
						rota[vertice_g] = vertice_d;
						rota[vertice_e] = vertice_b;
						rota[vertice_c] = vertice_h;
						return 0;
					}

					//Aplicar Double Brigde

					//Avanca uma areasta ( G, H)
					vertice_g = rota[vertice_g];
					vertice_h = rota[vertice_g];
				}

				//Avanca uma aresta (E, F)
				vertice_e = rota[vertice_e];
				vertice_f = rota[vertice_e];

			} //for (E, F)

			//avanca uma aresta (C, D)
			vertice_c = rota[vertice_c];
			vertice_d = rota[vertice_c];

		} //for (C, D)

		//Avanca proxima areasta (A, B)
		vertice_a = rota[vertice_a];
		vertice_b = rota[vertice_a];

	} //for (A, B)

	////printf("Qtd de Arestas = %d (vertices = %d)", count, tamanho_rota);

	//scanf("%d", &qtd_n);
	return 0;

}

int movimento_remover_re_insercao_mais_barata(int* vertices_solucao,
		int tamanho_rota, int qtd_n, int* &rota, Matriz* C, Vertice* vertices) {

	/*
	 * Retira e tentar realocar todos os vertices V da solucao
	 *  */
	//printf("Remover V e  Insercao V Mais Barata \n");
	int vertice_v = -1;
	int antecessor_vertice_v = -1;
	int sucessor_vertice_v = -1;
	float custo_atual_vertice_v = -1;
	float novo_custo_vertice_v = -1;
	//compatibilidade com insercao mais barata;
	int custo_min = -1;
	int custo_max = -1;
	//Matriz com cursos de insercao, Duas dimensoes: 0 - Posicao; 1 - Menor Custo
	Matriz *C_insert = create_matriz(2, qtd_n);

	int LC[qtd_n];
	for (int i = 0; i < qtd_n; ++i) {
		LC[i] = 0;
	}

	for (int i = 0; i < tamanho_rota; ++i) {

		vertice_v = vertices_solucao[i];
		antecessor_vertice_v = antecessor_de(rota, vertice_v);
		sucessor_vertice_v = rota[vertice_v];
		custo_atual_vertice_v = get(C, antecessor_vertice_v, vertice_v)
				+ get(C, vertice_v, sucessor_vertice_v)
				- get(C, antecessor_vertice_v, sucessor_vertice_v);

		//imprimir_solucao(rota, vertices);
		//printf("Custo Rota: %.0f \n", custo(C, rota));

		//remover v da solucao;
		remover_k_da_solucao(rota, vertice_v);
		LC[vertice_v] = 1;

		//Calcular o custo de inserir vertice_v
		definir_custos_dos_candidatos(vertices, C, rota, LC, C_insert,
				&custo_min, &custo_max, qtd_n);

		novo_custo_vertice_v = get(C_insert, 0, vertice_v);

		if (novo_custo_vertice_v < custo_atual_vertice_v) {

			//printf("Vertice v = %d\n", vertice_v);
			//printf("Custo Atual = %.0f\n", custo_atual_vertice_v);
			//printf("Novo Custo = %.0f\n", novo_custo_vertice_v);
			//printf("Posicao Novo Custo = %.0f\n", get(&C_insert, 1, vertice_v));
			inserir_k_na_solucao(rota, vertice_v, get(C_insert, 1, vertice_v));

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo(C, rota));
			free_matriz(C_insert);
			return 0;

			//scanf("Parou ! %d", &custo_max);

		} else {
			LC[vertice_v] = 0;
			inserir_k_na_solucao(rota, vertice_v, antecessor_vertice_v);

		}
	}
	free_matriz(C_insert);
	return 0;
}

int movimento_remover_insercao_mais_barata(int* vertices_solucao,
		int tamanho_rota, int d, int prize, int qtd_n, int* &rota, Matriz* C,
		Vertice* vertices) {

	/*
	 * Retira um vertice V/T e inserir um novo vertice V/T na solucao
	 *  */
	//printf("Remover V e Inserir novo K Mais Barato \n");
	int vertice_v = -1;
	int antecessor_vertice_v = -1;
	int sucessor_vertice_v = -1;
	float custo_atual_vertice_v = -1;
	float novo_custo_vertice_k = -1;

	//Matriz com custos de insercao, Duas dimensoes: 0 - Posicao; 1 - Menor Custo
	Matriz* C_insert = create_matriz(2, qtd_n);

	//compatibilidade com insercao mais barata;
	int k_custo_min = -1; //indice da melhor opcao
	int custo_max = -1;
	//int * vertices_w;
	int *vertices_habilitados;
	int qtd_habilitados = 0;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));

	//Criar uma LC com todos os Opcionais que esteja fora da rota
	int tamanho_lc = 0;
	int LC[qtd_n];
	int LC_Habilitados[qtd_n]; //Vertices (V/T) que realmente sao candidatos
	int vertices_vt[1];

	/* Gerar LC Old
	 for (int i = 0; i < qtd_n; ++i) {
	 if (strcmp(vertices[i].tipo,"v/t") == 0){
	 LC[i] = 1;
	 tamanho_lc++;
	 //printf("*%d, ", i);
	 }else {
	 LC[i] = 0;
	 }
	 }

	 for (int i = 0; i < tamanho_rota; ++i) {
	 if (strcmp(vertices[vertices_solucao[i]].tipo,"v/t") == 0){
	 LC[vertices_solucao[i]] = 0;
	 //printf("#%d, ", vertices_solucao[i]);
	 tamanho_lc--;
	 }
	 }
	 */

	tamanho_lc = obter_lc_dado_uma_solucao_viavel(vertices, vertices_solucao,
			tamanho_rota, LC, qtd_n);
	copy_int_vetor(rota, rota_original, qtd_n);

	if (tamanho_lc > 0)
		for (int i = 0; i < tamanho_rota; ++i) {

			vertice_v = vertices_solucao[i];
			//printf("Vertice V = %d\n", vertice_v);

			if (strcmp(vertices[vertice_v].tipo, "v/t") != 0)
				continue;

			/* 1) Retirar vertices V da solucao
			 * 2) Tornar candidato somente os vertice K cujo premio seja suficente para deixar S viável
			 * 3) E, que cubra os vertice que ficaram descobertos com saida de V
			 */

			antecessor_vertice_v = antecessor_de(rota, vertice_v);
			sucessor_vertice_v = rota[vertice_v];
			custo_atual_vertice_v = get(C, antecessor_vertice_v, vertice_v)
					+ get(C, vertice_v, sucessor_vertice_v)
					- get(C, antecessor_vertice_v, sucessor_vertice_v);
			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo(C, rota));

			//2) Ver qual premio que ainda é requerido
			//printf("Vertice V = %d \n", vertice_v);
			//printf("Prize Solicitado = %d \n", prize);
			//printf("Prize atual sem o vertice V = %d\n", prize_acumulado(vertices, rota, qtd_n));
			//int premio_requerido = prize - (prize_acumulado(vertices, rota, qtd_n));
			//printf("Prize Requerido = %d\n",  premio_requerido);

			//3) Obter vertices aptos a substituir vertice_v

			//copy_int_vetor(LC, LC_Habilitados, qtd_n);

			vertices_vt[0] = vertice_v;

			qtd_habilitados = obter_vertices_vt_habilitados(vertices, C,
					rota_original, tamanho_rota, qtd_n, vertices_vt,
					vertices_solucao, 1, prize, vertices_habilitados);

			if (qtd_habilitados == 0) {
				continue;
			}

			// Remover V da solucao;
			remover_k_da_solucao(rota, vertice_v);

			for (int i = 0; i < qtd_n; ++i) {
				LC_Habilitados[i] = 0;
			}
			for (int h = 0; h < qtd_habilitados; ++h) {
				LC_Habilitados[vertices_habilitados[h]] = 1;
				//printf("Vertice Habilitado = %d\n", vertices_habilitados[h]);
			}
			free(vertices_habilitados);

			//printf("\nQtd Habilitados/Candidato = %d/%d \n", qtd_habilitados, tamanho_lc);

			definir_custos_dos_candidatos(vertices, C, rota, LC_Habilitados,
					C_insert, &k_custo_min, &custo_max, qtd_n);

			novo_custo_vertice_k = get(C_insert, 0, k_custo_min);

			if (novo_custo_vertice_k < custo_atual_vertice_v) {
				//imprimir_solucao(rota, vertices);
				//printf("Vertice v = %d\n", vertice_v);
				//printf("Vertice k = %d\n", k_custo_min);
				//printf("Custo v = %.0f\n", custo_atual_vertice_v);
				//printf("Custo k= %.0f\n", novo_custo_vertice_k);
				//printf("Posicao Novo Custo = %d\n", (int)get(C_insert, 1, k_custo_min));
				inserir_k_na_solucao(rota, k_custo_min,
						(int) get(C_insert, 1, k_custo_min));

				//imprimir_solucao(rota, vertices);
				//printf("Custo Rota: %.0f \n", custo(C, rota));
				//imprimir_solucao_em_arquivo(rota, C, vertices);
				//LC[vertice_v] = 1;
				//LC[k_custo_min] = 0;
				free(rota_original);
				free_matriz(C_insert);
				//scanf("%d", &custo_max);
				//printf("Melhor Mov 7\n");

				return 0;

			} else { //Se nao há um vertices melhor, retorna v pra solucao
				//printf("Sem vertice habilitado que apresente melhora!\n");
				inserir_k_na_solucao(rota, vertice_v, antecessor_vertice_v);
			}

		}

	free(rota_original);
	free_matriz(C_insert);
	return 0;

}

int movimento_remover_re_inserir_genius(int* vertices_solucao, int tamanho_rota,
		int qtd_n, int* &rota, Matriz* C, Vertice* vertices) {

	/*
	 * Retira e tentar realocar todos os vertices V da solucao, via GENI
	 *  */
	//printf("Remover V e  RE-Insercao V via GENIUS \n");
	int vertice_v = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	copy_int_vetor(rota, rota_original, qtd_n);

	for (int i = 0; i < tamanho_rota; ++i) {

		vertice_v = vertices_solucao[i];
		custo_atual_rota = custo(C, rota);
		//imprimir_solucao(rota, vertices);
		//printf("Custo Rota: %.0f \n", custo_atual_rota);

		//remover v da solucao;
		//printf("Remover Vertice = %d\n", vertice_v);
		remover_k_da_solucao(rota, vertice_v);
		//imprimir_solucao(rota, vertices);

		//Re-inserir v via GENI
		//printf("Inserir Vertice = %d via GENI\n", vertice_v);
		GENI(vertices, rota, C, vertice_v, 5, qtd_n);
		//imprimir_solucao(rota, vertices);

		novo_custo_rota = custo(C, rota);

		if (novo_custo_rota < custo_atual_rota) {
			//imprimir_solucao(rota_original, vertices);
			/*printf("Vertice V = %d\n", vertice_v);
			 printf("Custo Atual = %.0f\n", custo_atual_rota);
			 printf("Novo Custo = %.0f\n", novo_custo_rota);

			 imprimir_solucao(rota, vertices);
			 printf("Custo Rota: %.0f \n", custo(C, rota));

			 //scanf("Parou ! %d", &qtd_n);*/
			free(rota_original);
			return 0;
		} else {

			copy_int_vetor(rota_original, rota, qtd_n);

		}
	}
	free(rota_original);
	return 0;
}

int movimento_remover_inserir_genius(int* vertices_solucao, int tamanho_rota,
		int d, int prize, int qtd_n, int* &rota, Matriz* C, Vertice* vertices) {

	/*
	 * Retira um vertice V/T e inserir um novo vertice V/T na solucao
	 *  */
	//printf("Remover V e Inserir novo K via GENIUS \n");
	int vertice_v = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	int *copia_rota = (int *) malloc(qtd_n * sizeof(int));
	int *rota_pos_remocao = (int *) malloc(qtd_n * sizeof(int));
	int vertices_vt[1];

	copy_int_vetor(rota, rota_original, qtd_n);

	//Criar uma LC com todos os Opcionais que esteja fora da rota
	int tamanho_lc = 0;
	int *vertices_habilitados;
	int qtd_habilitados = 0;
	int LC[qtd_n];

	tamanho_lc = obter_lc_dado_uma_solucao_viavel(vertices, vertices_solucao,
			tamanho_rota, LC, qtd_n);

	custo_atual_rota = custo(C, rota);

	if (tamanho_lc > 0)
		for (int i = 0; i < tamanho_rota; ++i) {

			vertice_v = vertices_solucao[i];
			//printf("Vertice V = %d\n", vertice_v);

			if (strcmp(vertices[vertice_v].tipo, "t") == 0)
				continue;
			/* 1) Retirar vertices V/T da solucao
			 * 2) Tornar candidato somente os vertice K cujo premio seja suficente para deixar S viável
			 * 3) E, que cubra os vertice que ficaram descobertos com saida de V
			 */

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo_atual_rota);
			//printf("Remover Vertice  = %d", vertice_v);
			//2) Ver qual premio que ainda é requerido
			//printf("Vertice V = %d \n", vertice_v);
			//printf("Prize Solicitado = %d \n", prize);
			//printf("Prize atual sem o vertice V = %d\n", prize_acumulado(vertices, rota, qtd_n));
			//int premio_requerido = prize - (prize_acumulado(vertices, rota, qtd_n));
			//printf("Prize Requerido = %d\n",  premio_requerido);
			copy_int_vetor(rota, copia_rota, qtd_n);

			//Pegar os Vertices Opcionais da Rota
			vertices_vt[0] = vertice_v;

			qtd_habilitados = obter_vertices_vt_habilitados(vertices, C,
					copia_rota, tamanho_rota, qtd_n, vertices_vt,
					vertices_solucao, 1, prize, vertices_habilitados);

			if (qtd_habilitados == 0) {
				continue;
			}

			//1) Remover V da solucao;
			remover_k_da_solucao(rota, vertice_v);

			copy_int_vetor(rota, rota_pos_remocao, qtd_n);

			//printf("\nQtd Habilitados/Candidato = %d/%d \n", qtd_habilitados, tamanho_lc);

			for (int h = 0; h < qtd_habilitados; ++h) {

				GENI(vertices, rota, C, vertices_habilitados[h], 5, qtd_n);

				novo_custo_rota = custo(C, rota);

				if (novo_custo_rota < custo_atual_rota) {

					//imprimir_solucao(rota_original, vertices);

					//printf("Vertice v = %d\n", vertice_v);
					//printf("Custo Rota = %.0f\n", custo_atual_rota);
					//printf("Custo Novo= %.0f\n", novo_custo_rota);

					//imprimir_solucao(rota, vertices);
					//printf("Custo Rota: %.0f \n", custo(C, rota));
					//imprimir_solucao_em_arquivo(rota, C, vertices);

					free(copia_rota);
					free(rota_original);
					free(rota_pos_remocao);
					free(vertices_habilitados);
					return 0;

				} else { //Se nao há um vertices melhor, retorna v pra solucao
					//printf("Vertice Habilitado (%d) nao melhorou a rota via GENI!\n", vertice_h);
					copy_int_vetor(rota_pos_remocao, rota, qtd_n);
				}
			}

			//Se nao há habilitados, retorna v pra solucao,
			//printf("Sem vertice habilitado!, para substituir o vertice = %d\n", vertice_v);
			copy_int_vetor(rota_original, rota, qtd_n);
			free(vertices_habilitados);

		}

	free(copia_rota);
	free(rota_original);
	free(rota_pos_remocao);

	return 0;

}

int movimento_remover_mais_barata_re_inserir_genius(int* vertices_solucao,
		int tamanho_rota, int d, int prize, int qtd_n, int* &rota, Matriz* C,
		Vertice* vertices) {

	//Pegar o vertice de maior custo e reinseri-lo via GENI
	int custo_vertice_v = -1;
	int maior_custo = 0;
	int vertice_v = -1;
	int vertice_maior_custo = -1;
	int antecessor_vertice_v = -1;
	int sucessor_vertice_v = -1;

	//Pegar vertice de maior custo
	for (int i = 0; i < tamanho_rota; ++i) {
		vertice_v = vertices_solucao[i];
		antecessor_vertice_v = antecessor_de(rota, vertice_v);
		sucessor_vertice_v = rota[vertice_v];
		custo_vertice_v = get(C, antecessor_vertice_v, vertice_v)
				+ get(C, vertice_v, sucessor_vertice_v)
				- get(C, antecessor_vertice_v, sucessor_vertice_v);
		if (custo_vertice_v > maior_custo) {
			maior_custo = custo_vertice_v;
			vertice_maior_custo = vertice_v;
		}
		//printf("Vertice V = %d, Custo = %d (Premio=%d) \n", vertice_v, custo_vertice_v, vertices[vertice_v].premio);

	}

	//printf("Maior Custo: Vertice V = %d, Custo = %d", vertice_maior_custo, maior_custo);
	int custo_antes = custo(C, rota);
	remover_k_da_solucao(rota, vertice_maior_custo);

	GENI(vertices, rota, C, vertice_maior_custo, 5, qtd_n);
	int custo_depois = custo(C, rota);

	if (custo_depois < custo_antes) {

	}
	//printf("TODO!");
	//scanf("%d", &qtd_n);

	return 0;

}

int movimento_remover_genius_re_inserir_genius(int* vertices_solucao,
		int tamanho_rota, int qtd_n, int* &rota, Matriz* C, Vertice* vertices) {

	/*
	 * Retira e tentar realocar todos os vertices V da solucao
	 *  */
	//printf("Remover GENIUS  e RE-Inserir GENIUS \n");
	int vertice_v = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	copy_int_vetor(rota, rota_original, qtd_n);

	for (int i = 0; i < tamanho_rota; ++i) {

		vertice_v = vertices_solucao[i];

		//printf("Vertice V = %d\n", vertice_v);

		custo_atual_rota = custo(C, rota);
		//imprimir_solucao(rota, vertices);

		//printf("Custo Rota: %.0f \n", custo_atual_rota);

		//remover v da solucao com um dos tipo de GENIUS;

		int retorno = remover_genius(vertices, rota, vertice_v, C, qtd_n,
				tamanho_rota, vertices_solucao);

		if (retorno == 1) {
			//printf("GENIUS nao conseguiu remover o vertice\n");
			//free(rota_original);
			continue;
			//return 0;
			//scanf("%d", &qtd_n);
		}

		//imprimir_solucao(rota, vertices);

		//Re-inserir v via GENI
		GENI(vertices, rota, C, vertice_v, 5, qtd_n);

		novo_custo_rota = custo(C, rota);

		if (novo_custo_rota < custo_atual_rota) {
			//imprimir_solucao(rota_original, vertices);
			//printf("Vertice V = %d\n", vertice_v);
			//printf("Custo Atual = %.0f\n", custo_atual_rota);
			//printf("Novo Custo = %.0f\n", novo_custo_rota);

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo(C, rota));

			//scanf("Parou ! %d", &qtd_n);
			free(rota_original);
			return 0;
		} else {

			copy_int_vetor(rota_original, rota, qtd_n);

		}

	}

	free(rota_original);
	return 0;

}

int movimento_remover_genius_inserir_genius(int* vertices_solucao,
		int tamanho_rota, int d, int prize, int qtd_n, int* &rota, Matriz* C,
		Vertice* vertices) {

	/*
	 * Retira um vertice V/T via US e inserir um novo vertice V/T na solucao via GENI,
	 *  */
	//printf("Remover GENIUS V e Inserir novo K via GENIUS \n");
	//imprimir_solucao(rota, vertices);
	int vertices_vt[1];
	int vertice_v = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	int *rota_pos_remocao = (int *) malloc(qtd_n * sizeof(int));
	int *copia_rota = (int *) malloc(qtd_n * sizeof(int));
	int *vertices_habilitados;

	copy_int_vetor(rota, rota_original, qtd_n);

	//Criar uma LC com todos os Opcionais que esteja fora da rota
	int tamanho_lc = 0;
	int qtd_habilitados = 0;
	int LC[qtd_n];

	tamanho_lc = obter_lc_dado_uma_solucao_viavel(vertices, vertices_solucao,
			tamanho_rota, LC, qtd_n);

	custo_atual_rota = custo(C, rota);

	if (tamanho_lc > 0)
		for (int i = 0; i < tamanho_rota; ++i) {

			vertice_v = vertices_solucao[i];
			//printf("Vertice V = %d\n", vertice_v);

			if (strcmp(vertices[vertice_v].tipo, "t") == 0)
				continue;
			/* 1) Retirar vertices V/T da solucao
			 * 2) Tornar candidato somente os vertice K cujo premio seja suficente para deixar S viável
			 * 3) E, que cubra os vertice que ficaram descobertos com saida de V
			 */

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo_atual_rota);
			//printf("Remover Vertice  = %d\n", vertice_v);
			//1) Remover V da solucao;
			int retorno = remover_genius(vertices, rota, vertice_v, C, qtd_n,
					tamanho_rota, vertices_solucao);

			if (retorno == 1) {
				//printf("GENIUS nao conseguiu remover o vertice\n");
				continue;
				//free(rota_original);
				//return 0;
				//scanf("%d", &qtd_n);
			}
			//printf("Vertice Removido  = %d\n", vertice_v);

			copy_int_vetor(rota, rota_pos_remocao, qtd_n);

			vertices_vt[0] = vertice_v;

			//printf("V/T1 e V/T2: %d ** %d  \n", vertice_a, vertice_b);

			//Fazer uma copia para rota para ser usada na obtencao V/T substitutos
			copy_int_vetor(rota_original, copia_rota, qtd_n);
			//Pegar os Vertices Opcionais da Rota

			qtd_habilitados = obter_vertices_vt_habilitados(vertices, C,
					copia_rota, tamanho_rota, qtd_n, vertices_vt,
					vertices_solucao, 1, prize, vertices_habilitados);

			//printf("\nQtd Habilitados/Candidato = %d/%d \n", qtd_habilitados, tamanho_lc);
			if (qtd_habilitados == 0) {
				copy_int_vetor(rota_original, rota, qtd_n);
				continue;
			}

			for (int h = 0; h < qtd_habilitados; ++h) {

				GENI(vertices, rota, C, vertices_habilitados[h], 5, qtd_n);

				novo_custo_rota = custo(C, rota);

				if (novo_custo_rota < custo_atual_rota) {

					//imprimir_solucao(rota_original, vertices);

					//printf("Vertice v = %d\n", vertice_v);
					//printf("Custo Rota = %.0f\n", custo_atual_rota);
					//printf("Custo Novo= %.0f\n", novo_custo_rota);

					//imprimir_solucao(rota, vertices);
					//printf("Custo Rota: %.0f \n", custo(C, rota));
					//imprimir_solucao_em_arquivo(rota, C, vertices);
					//printf("Parou... Remover Inserir GENIUS");
					//scanf("Parou ! %d", &qtd_n);

					free(rota_original);
					free(rota_pos_remocao);
					free(vertices_habilitados);
					free(copia_rota);
					return 0;

				} else { //Se nao há um vertices melhor, retorna v pra solucao
					//printf("Vertice Habilitado (%d) nao melhorou a rota via GENI!\n", vertice_h);
					copy_int_vetor(rota_pos_remocao, rota, qtd_n);
					//free(vertices_habilitados);
				}
			}

			//Se nao há habilitados, retorna v pra solucao,
			//printf("Sem vertice habilitado!, para substituir o vertice = %d\n", vertice_v);
			copy_int_vetor(rota_original, rota, qtd_n);
			free(vertices_habilitados);
		}

	free(rota_original);
	free(rota_pos_remocao);
	//free(vertices_habilitados);
	free(copia_rota);
	return 0;

}

int movimento_remover_genius_re_insercao_mais_barata(int* vertices_solucao,
		int tamanho_rota, int qtd_n, int* &rota, Matriz* C, Vertice* vertices) {

	/*
	 * Retira e tentar realocar todos os vertices V da solucao
	 *  */
	//printf("Remover GENIUS  e  RE-Insercao V Mais Barata \n");
	int vertice_v = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	//compatibilidade com insercao mais barata;
	int custo_min = -1;
	int custo_max = -1;
	//Matriz com cursos de insercao, Duas dimensoes: 0 - Posicao; 1 - Menor Custo
	Matriz *C_insert = create_matriz(2, qtd_n);
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	copy_int_vetor(rota, rota_original, qtd_n);

	int LC[qtd_n];
	for (int i = 0; i < qtd_n; ++i) {
		LC[i] = 0;
	}

	custo_atual_rota = custo(C, rota);

	for (int i = 0; i < tamanho_rota; ++i) {

		vertice_v = vertices_solucao[i];

		//printf("Vertice V = %d\n", vertice_v);

		//imprimir_solucao(rota, vertices);
		//printf("Custo Rota antes da remocao-genius: %.0f \n", custo(C, rota));

		//printf("Remover Vertice = %d \n", vertice_v);
		//remover v da solucao com um dos tipo de GENIUS;
		int retorno = remover_genius(vertices, rota, vertice_v, C, qtd_n,
				tamanho_rota, vertices_solucao);

		if (retorno == 1) {
			//printf("GENIUS nao conseguiu remover!\n");
			//free(rota_original);
			//free_matriz(C_insert);
			continue;
		}

		//imprimir_solucao(rota, vertices);

		LC[vertice_v] = 1;

		definir_custos_dos_candidatos(vertices, C, rota, LC, C_insert,
				&custo_min, &custo_max, qtd_n);

		inserir_k_na_solucao(rota, vertice_v, get(C_insert, 1, vertice_v));

		novo_custo_rota = custo(C, rota);

		if (novo_custo_rota < custo_atual_rota) {
			//imprimir_solucao(rota_original, vertices);
			//printf("Vertice v = %d\n", vertice_v);
			//printf("Custo Atual Rota = %.0f\n", custo_atual_rota);
			//printf("Novo Custo Rota = %.0f\n", novo_custo_rota);
			//printf("Posicao Novo Custo = %.0f\n", get(C_insert, 1, vertice_v));

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo(C, rota));

			free(rota_original);
			free_matriz(C_insert);
			return 0;

			//scanf("Parou ! %d", &custo_max);

		} else {
			LC[vertice_v] = 0;
			copy_int_vetor(rota_original, rota, qtd_n);
		}
	}

	free(rota_original);
	free_matriz(C_insert);
	return 0;

}

int movimento_remover_genius_inserir_mais_barata(int* vertices_solucao,
		int tamanho_rota, int d, int prize, int qtd_n, int* &rota, Matriz* C,
		Vertice* vertices) {

	/*
	 * Retira um vertice V/T e inserir um novo vertice V/T na solucao
	 *  */
	//printf("Remover GENIUS V e Inserir novo K via MAIS BARATA\n");
	int vertice_v = -1;
	int vertices_vt[1];
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));
	int *copia_rota = (int *) malloc(qtd_n * sizeof(int));
	int *rota_pos_remocao = (int *) malloc(qtd_n * sizeof(int));
	int *vertices_habilitados;

	//Matriz com custos de insercao, Duas dimensoes: 0 - Posicao; 1 - Menor Custo
	Matriz* C_insert = create_matriz(2, qtd_n);

	//compatibilidade com insercao mais barata;
	int k_custo_min = -1; //indice da melhor opcao
	int custo_max = -1;

	copy_int_vetor(rota, rota_original, qtd_n);

	//Criar uma LC com todos os Opcionais que esteja fora da rota
	int tamanho_lc = 0;
	int qtd_habilitados = 0;
	int LC[qtd_n];
	int LC_Habilitados[qtd_n]; //Vertices (V/T) que realmente sao candidatos

	//Criar LC dado um rota
	tamanho_lc = obter_lc_dado_uma_solucao_viavel(vertices, vertices_solucao,
			tamanho_rota, LC, qtd_n);

	custo_atual_rota = custo(C, rota);

	if (tamanho_lc > 0)
		for (int i = 0; i < tamanho_rota; ++i) {

			vertice_v = vertices_solucao[i];
			//printf("Vertice V = %d\n", vertice_v);

			if (strcmp(vertices[vertice_v].tipo, "t") == 0)
				continue;
			/* 1) Retirar vertices V da solucao
			 * 2) Tornar candidato somente os vertice K cujo premio seja suficente para deixar S viável
			 * 3) E, que cubra os vertice que ficaram descobertos com saida de V
			 */

			//imprimir_solucao(rota, vertices);
			//printf("Custo Rota: %.0f \n", custo_atual_rota);
			//printf("\nRemover Vertice  = %d, premio = %d\n", vertice_v, vertices[vertice_v].premio);
			//printf("Premio atual: %d\n", prize_acumulado(vertices, rota, qtd_n));
			//1) Remover V da solucao;
			int retorno = remover_genius(vertices, rota, vertice_v, C, qtd_n,
					tamanho_rota, vertices_solucao);

			if (retorno == 1) {
				//printf("GENIUS nao conseguiu remover o vertice\n");
				continue;
				//free(rota_original);
				//free(C_insert);
				//return 0;
				//scanf("%d", &qtd_n);
			}
			//printf("Vertice Removido!\n");
			//printf("Premio atual: %d\n", prize_acumulado(vertices, rota, qtd_n));
			//printf("Premio Requerido: %d\n", prize-prize_acumulado(vertices, rota, qtd_n));

			copy_int_vetor(rota, rota_pos_remocao, qtd_n);

			vertices_vt[0] = vertice_v;

			//Fazer uma copia para rota para ser usada na obtencao V/T substitutos
			copy_int_vetor(rota_original, copia_rota, qtd_n);

			qtd_habilitados = obter_vertices_vt_habilitados(vertices, C,
					copia_rota, tamanho_rota, qtd_n, vertices_vt,
					vertices_solucao, 1, prize, vertices_habilitados);

			//printf("Qtd Habilitados/Candidato = %d/%d \n", qtd_habilitados, tamanho_lc);

			if (qtd_habilitados == 0) {
				copy_int_vetor(rota_original, rota, qtd_n);
				continue;
			}

			//MAndar os habilitados para pegar o que tiver o melhor custo de insercao (TÁ FEI!!)
			for (int i = 0; i < qtd_n; ++i) {
				LC_Habilitados[i] = 0;
			}
			for (int h = 0; h < qtd_habilitados; ++h) {
				LC_Habilitados[vertices_habilitados[h]] = 1;
			}
			free(vertices_habilitados);

			definir_custos_dos_candidatos(vertices, C, rota_pos_remocao,
					LC_Habilitados, C_insert, &k_custo_min, &custo_max, qtd_n);

			inserir_k_na_solucao(rota_pos_remocao, k_custo_min,
					get(C_insert, 1, k_custo_min));

			novo_custo_rota = custo(C, rota_pos_remocao);

			if (novo_custo_rota < custo_atual_rota) {

				//printf("Vertice v = %d\n", vertice_v);
				//printf("Vertice k = %d, premio = %d\n", k_custo_min, vertices[k_custo_min].premio);
				//printf("Premio atual: %d\n", prize_acumulado(vertices, rota_pos_remocao, qtd_n));

				//printf("Custo Atual Rota = %.0f\n", custo_atual_rota);
				//printf("Novo Custo com k= %.0f\n", novo_custo_rota);
				//printf("Posicao Novo Custo = %.0f\n", get(C_insert, 1, k_custo_min));

				//imprimir_solucao(rota, vertices);
				//printf("Custo Rota: %.0f \n", custo(C, rota));
				//imprimir_solucao_em_arquivo(rota, C, vertices);
				//LC[vertice_v] = 1;
				//LC[k_custo_min] = 0;
				copy_int_vetor(rota_pos_remocao, rota, qtd_n);

				free_matriz(C_insert);
				free(rota_original);
				free(copia_rota);
				free(rota_pos_remocao);
				return 0;
				//scanf("Parou ! %d", &custo_max);

			} else { //Se nao há um vertices melhor, retorna v pra solucao
				//printf("Sem vertice habilitado que apresente melhora!\n");
				copy_int_vetor(rota_original, rota, qtd_n);
			}

		}

	free(rota_original);
	free_matriz(C_insert);
	free(copia_rota);
	free(rota_pos_remocao);
	return 0;

}

void movimento_dupla_remocao_simples_inserir_mais_barata(int* vertices_solucao,
		int tamanho_rota, int d, int prize, int qtd_n, int* &rota, Matriz* C,
		Vertice* vertices) {
	/* DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA: Pegar procurar pares de vertices opcionais (V/T) que possam
	 * substituidos por algum V/T-Candidato
	 */

	int vertices_vt[2];
	int vertice_a = -1;
	int vertice_b = -1;

	int *vertices_habilitados;
	int qtd_habilitados = 0;
	int *copia_rota = (int *) malloc(qtd_n * sizeof(int));
	int *rota_original = (int *) malloc(qtd_n * sizeof(int));

	copy_int_vetor(rota, rota_original, qtd_n);

	//Matriz com custos de insercao, Duas dimensoes: 0 - Posicao; 1 - Menor Custo
	Matriz* C_insert = create_matriz(2, qtd_n);

	//compatibilidade com insercao mais barata;
	int k_custo_min = -1; //indice da melhor opcao
	int custo_max = -1;
	float custo_atual_rota = -1;
	float novo_custo_rota = -1;
	float custo_sem_a_e_b = -1;

	int LC_Habilitados[qtd_n];
	for (int i = 0; i < qtd_n; ++i) {
		LC_Habilitados[i] = 0;
	}

	//imprimir_solucao(rota, vertices);

	for (int i = 0; i < tamanho_rota; ++i) {
		vertice_a = vertices_solucao[i];

		if (strcmp(vertices[vertice_a].tipo, "v/t") != 0)
			continue;

		for (int j = (i + 1); j < tamanho_rota; ++j) {
			vertice_b = vertices_solucao[j];

			if (strcmp(vertices[vertice_b].tipo, "v/t") != 0 or i == j)
				continue;

			vertices_vt[0] = vertice_a;
			vertices_vt[1] = vertice_b;

			//printf("Vertices a remover V/T1 e V/T2: %d ** %d  \n", vertice_a, vertice_b);

			//Fazer uma copia para rota para ser usada na obtencao V/T substitutos
			copy_int_vetor(rota, copia_rota, qtd_n);
			//Pegar os Vertices Opcionais da Rota

			qtd_habilitados = obter_vertices_vt_habilitados(vertices, C,
					copia_rota, tamanho_rota, qtd_n, vertices_vt,
					vertices_solucao, 2, prize, vertices_habilitados);

			if (qtd_habilitados == 0) {
				continue;
			}

			custo_atual_rota = custo(C, rota);
			//printf("Premio atual: %d\n", prize_acumulado(vertices, rota, qtd_n));

			//imprimir_solucao(rota, vertices);
			//printf("Remover Vertice = %d, premio = %d\n", vertice_a, vertices[vertice_a].premio);
			remover_k_da_solucao(rota, vertice_a);
			//imprimir_solucao(rota, vertices);
			//printf("Remover Vertice = %d, premio = %d\n", vertice_b, vertices[vertice_b].premio);
			remover_k_da_solucao(rota, vertice_b);
			//imprimir_solucao(rota, vertices);

			custo_sem_a_e_b = custo(C, rota);
			//printf("Premio atual: %d\n", prize_acumulado(vertices, rota, qtd_n));

			for (int i = 0; i < qtd_n; ++i) {
				LC_Habilitados[i] = 0;
			}
			for (int h = 0; h < qtd_habilitados; ++h) {
				LC_Habilitados[vertices_habilitados[h]] = 1;
				//printf("Vertice Habilitado = %d\n", vertices_habilitados[h]);
			}
			free(vertices_habilitados);

			//INSERIR MAIS BARATA
			definir_custos_dos_candidatos(vertices, C, rota, LC_Habilitados,
					C_insert, &k_custo_min, &custo_max, qtd_n);

			novo_custo_rota = custo_sem_a_e_b + get(C_insert, 0, k_custo_min);

			if (novo_custo_rota < custo_atual_rota) {

				//printf("Vertice k = %d, premio %d\n", k_custo_min, vertices[k_custo_min].premio);
				//printf("Custo v = %.0f\n", custo_atual_vertice_v);
				//printf("Custo k= %.0f\n", novo_custo_vertice_k);
				//printf("Posicao Novo Custo = %.0f\n", get(C_insert, 1, k_custo_min));
				/*printf("Removendo A=%d e B=%d e Inserindo K=%d \n", vertice_a, vertice_b, (int)k_custo_min);
				 imprimir_solucao(rota_original, vertices);
				 imprimir_solucao(rota, vertices);*/
				inserir_k_na_solucao(rota, (int) k_custo_min,
						(int) get(C_insert, 1, k_custo_min));
				//printf("Premio atual: %d\n", prize_acumulado(vertices, rota, qtd_n));
				//imprimir_solucao(rota, vertices);
				//printf("Custo Rota: %.0f \n", custo(C, rota));
				//imprimir_solucao_em_arquivo(rota, C, vertices);
				//LC[vertice_v] = 1;
				//LC[k_custo_min] = 0;
				//printf("Parou ! ok\n");
				//scanf("%d", &custo_max);
				free_matriz(C_insert);
				free(copia_rota);
				free(rota_original);
				return;
			} else { //Se nao há um vertices melhor, retorna v pra solucao
				//printf("Sem vertice habilitado que apresente melhora!\n");
				//imprimir_solucao(rota_original, vertices);
				//imprimir_solucao(rota, vertices);
				//Desmarca como candidato
				//LC_Habilitados[vertices_habilitados[h]] = 0;
				//printf("Parou ! nok\n");
				//scanf("%d", &custo_max);
				copy_int_vetor(rota_original, rota, qtd_n);
			}

		}
	}

	free_matriz(C_insert);
	free(copia_rota);
	free(rota_original);

}

int avaliar_solucao(Vertice *vertices, Matriz *C, int *rota, int prize, float d,
		int qtd_n) {
	/* Avaliar se uma dada solucao é viável, ou sejam
	 * 1) Atinge o PRIZE solicitado
	 * 2) Contem todos os vertices T
	 * 3) Todos os vertices W estejam cobertos
	 */

	//printf("Avaliando a Solução...\n");
	//1 Prize Solicitado
	if (prize_acumulado(vertices, rota, qtd_n) < prize) {
		printf("Restricao de Premio Violada! \n");
		printf("Premio Acumulado/Prize %d / %d \n",
				prize_acumulado(vertices, rota, qtd_n), prize);
		return -1;
	}

	// printf("PREMIO OK!\n");

	if (existe_t_fora_da_solucao(vertices, rota, qtd_n) != 0) {
		printf("Restricao de Vertices Obrigatorios Violada! \n");
		return -1;
	}

	// printf("VERTICES T OK! \n");

	if (existe_w_descoberto(vertices, C, d, rota, qtd_n) != 0) {
		// printf("Restricao de Vertices de Cobertura Violada! \n");
		return -1;
	}

	for (int i = 0; i < qtd_n; ++i) {
		if (strcmp(vertices[i].tipo, "w") == 0 and rota[i] != -1) {
			// printf("Vértice W na rota! \n");
			return -1;
		}

	}

	// printf("VERTICES W OK! \n");

	return (int) custo(C, rota);
}

void obter_vizinho(Vertice *vertices, Matriz* C, int mov, int *&rota,
		int *vertices_solucao, int tamanho_rota, int d, int prize, int qtd_n) {
//	bool debug = false;

//		printf("\n Vizinhanca: %d - %s", mov, nome_movimento(mov).c_str());
	switch (mov) {
	case SHIFT:
		movimento_shift_v2(vertices_solucao, tamanho_rota, qtd_n, rota, C,
				vertices);
		break;
	case SWAP:
		movimento_swap(vertices_solucao, tamanho_rota, qtd_n, rota, C,
				vertices);
		break;
	case OR_OPT:
		movimento_or_opt(vertices_solucao, tamanho_rota, qtd_n, rota, C,
				vertices);
		break;
	case TWO_OPT:
		movimento_two_opt(vertices_solucao, tamanho_rota, qtd_n, rota, C,
				vertices);
		break;
	case THREE_OPT:
		movimento_three_opt(vertices_solucao, tamanho_rota, qtd_n, rota, C,
				vertices);
		break;
	case DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA:
		movimento_dupla_remocao_simples_inserir_mais_barata(vertices_solucao,
				tamanho_rota, d, prize, qtd_n, rota, C, vertices);
		break;
	case REMOVER_RE_INSERCAO_MAIS_BARATA:
		movimento_remover_re_insercao_mais_barata(vertices_solucao,
				tamanho_rota, qtd_n, rota, C, vertices);
		break;
	case REMOVER_INSERCAO_MAIS_BARATA:
		movimento_remover_insercao_mais_barata(vertices_solucao, tamanho_rota,
				d, prize, qtd_n, rota, C, vertices);
		break;
	case REMOVER_RE_INSERIR_GENIUS:
		movimento_remover_re_inserir_genius(vertices_solucao, tamanho_rota,
				qtd_n, rota, C, vertices);
		break;
	case REMOVER_INSERIR_GENIUS:
		movimento_remover_inserir_genius(vertices_solucao, tamanho_rota, d,
				prize, qtd_n, rota, C, vertices);
		break;
	case REMOVER_MAIS_BARATA_RE_INSERIR_GENIUS:
		//Em stand by
		//movimento_remover_mais_barata_re_inserir_genius(vertices_solucao, tamanho_rota, d, prize, qtd_n,	rota, C, vertices);
		break;
	case REMOVER_GENIUS_RE_INSERIR_GENIUS:
		movimento_remover_genius_re_inserir_genius(vertices_solucao,
				tamanho_rota, qtd_n, rota, C, vertices);
		break;
	case REMOVER_GENIUS_INSERIR_GENIUS:
		movimento_remover_genius_inserir_genius(vertices_solucao, tamanho_rota,
				d, prize, qtd_n, rota, C, vertices);
		break;
	case REMOVER_GENIUS_RE_INSERCAO_MAIS_BARATA:
		movimento_remover_genius_re_insercao_mais_barata(vertices_solucao,
				tamanho_rota, qtd_n, rota, C, vertices);
		break;
	case REMOVER_GENIUS_INSERCAO_MAIS_BARATA:
		movimento_remover_genius_inserir_mais_barata(vertices_solucao,
				tamanho_rota, d, prize, qtd_n, rota, C, vertices);
		break;
	}

}

void busca_local(Vertice *vertices, Matriz* C, int *rota, int *s_melhor_vizinho,
		int d, int prize, int qtd_n) {
	//TODO
	/* Alg. de Busca Local: (MRD - Metodo Randomico de Descina) (VND - Variable Neighborhood Descent) (MAtheus e Motta Usaram)
	 * 1) Executar X iteracoes, em cada iteracao gera-se um vizinho (aleatoriamente) de s_inicial.
	 * 2) Ao encotrar vizinho melhor, X volta a 1 (ou 0).
	 * 3) Busca Local se encerra ao executar X iteracoes sem achar um vizinho melhor.
	 */

	/* Alg. de Busca Local: VNRD (Metodo de Descida Randominca em Vizinhanca Variavel) (Matheus Usou)
	 * Usa um conjunto de vizinhanca.
	 * 1) Pega estrura de visinhaca N1;
	 * 2) Gera-se um vizinho usando N1;
	 * 3) Após X iteracoes sem achar vizinho melhor, passa-se pra proxima estrutura de
	 *    vizinhanca N2, N3, Nx;
	 * 4) A qualquer momento que for achado um vizinho melhor  que solucao atual volta-se a X = 0
	 * 5) Ao encerrar X iteracoes, verifica-se o vizinho encontrado e melhor que melhor solucao, se
	 *    positivo volta-se à primeira Estrutura de Vizinanc
	 * 6) A Busca local so termina qdo se percorre X Iteracoes em todas as Estruturas sem Melhora
	 */

	/*
	 * RVND - Random VND
	 * Após X iteracoes sem melhora em uma vizinhança N, seleciona outra randomicamente.
	 */

	for (int i = 0; i < 0; ++i) {
		int LC[qtd_n];
		criar_LC_4_genius(vertices, rota, LC, qtd_n, 0);

		int p = 5;
		//Vertices Candidatos e qtd
		int * vertices_lc = (int *) malloc(sizeof(int));

		int qtd_lc;
		qtd_lc = vertices_de_LC(LC, vertices_lc, qtd_n);
		if (qtd_lc > 0) {
			int vertice_v = criar_LRC_4_genius(vertices, vertices_lc, qtd_lc,
					rota, qtd_n, C, 0);
			GENI(vertices, rota, C, vertice_v, p, qtd_n);
			printf("\n inserindo vértice %d", vertice_v);
		}

	}

	//int *s_inicial_base = (int *) (malloc(qtd_n * sizeof(int)));
	//copy_int_vetor(rota, s_inicial_base, qtd_n);
	copy_int_vetor(rota, s_melhor_vizinho, qtd_n);

	int melhor_custo = custo(C, s_melhor_vizinho);

	//Botar vertices da solucao em um vetor para realizar os sorteios e definir tamanho da solucao
	int* vertices_solucao;

	int tamanho_da_rota = vertices_da_solucao(vertices_solucao,
			s_melhor_vizinho);

	int ordem = 0; //Ordem das vizinhancas, controla o avancar da lista de vizinhancas

	int vizinhanca = 0; //Vizinhanca corrente

	int qtd_vizinhanca = 15;

	int aleatorio = 1; //RVND(1) ou VND(0)

	int* s_vizinho;

	s_vizinho = (int*) (malloc(qtd_n * sizeof(int)));

	copy_int_vetor(s_melhor_vizinho, s_vizinho, qtd_n);

	int vizinhancas[qtd_vizinhanca];

	int contador = 0;

	int random;

	for (int i = 0; i < qtd_vizinhanca; ++i) {
		random = (genrand_int32() % qtd_vizinhanca);
		while (in_int_vetor(random, vizinhancas, contador) == 1)
			random = (genrand_int32() % qtd_vizinhanca);
		vizinhancas[contador++] = random;
	}

	//BUSCA LOCAL: VND
	while (ordem < qtd_vizinhanca) {
		//printf("Vizinhanca: %d \n", vizinhanca);
		//printf("Maria\n");

		//Pegar primeiro melhor vizinho da vizinhanca
		if (aleatorio == 1) {
			vizinhanca = vizinhancas[ordem];
		} else {
			vizinhanca = ordem; //sequencialmente apartir da ZERO
		}

		obter_vizinho(vertices, C, vizinhanca, s_vizinho, vertices_solucao,
				tamanho_da_rota, d, prize, qtd_n);

		/* Se o melhor vizinho da vizinhanca recem explorada for melhor
		 * que melhor vizinho jah encontrado:
		 * 1) Esse vizinho torna-se o s_melhor_vizinho
		 * 2) Volta-se a primeira Estrutura de Vizinhanca
		 */
		/*
		 if (avaliar_solucao(vertices, C, s_vizinho, prize, d, qtd_n) == -1) {
		 printf("Busca Local: Vizinhanca: %d\tNovo custo: %.0f\n", vizinhanca , custo( C, s_vizinho));
		 printf("Solução Inválida\n");
		 scanf("%d", &qtd_n);
		 }*/

		if (custo(C, s_vizinho) < melhor_custo) {
			copy_int_vetor(s_vizinho, s_melhor_vizinho, qtd_n);

			// printf(" --> (novo custo: %.0f)", custo(C, s_vizinho));

			melhor_custo = (int) custo(C, s_melhor_vizinho);

			//Para se achou custo exato
			if (melhor_custo <= valor_otimo)
				return;

			if (vizinhanca == REMOVER_INSERCAO_MAIS_BARATA
					or vizinhanca == REMOVER_INSERIR_GENIUS
					or vizinhanca == REMOVER_GENIUS_INSERIR_GENIUS
					or vizinhanca == REMOVER_GENIUS_INSERCAO_MAIS_BARATA
					or vizinhanca
							== DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA) {
				free(vertices_solucao);
				tamanho_da_rota = vertices_da_solucao(vertices_solucao,
						s_melhor_vizinho);
			}
			ordem = 0;
		} else {
			ordem++;
			//vizinhanca--;
		}
		//printf("Rogerio\n");

	}
	//fclose(file);

	//limpar_solucao(s_melhor_vizinho, vertices, qtd_n, C, d, prize);

	//BUSCA LOCAL: VNRD
	/*
	 //Quantas iteracoes
	 int BL_max_iter = 50;
	 int BL_iter = 0;
	 while (vizinhanca < qtd_vizinhanca){

	 BL_iter = 0;

	 while (BL_iter < BL_max_iter){

	 BL_iter++;

	 // atualiza s_vizinho_encontrado com s_inicial_base
	 copy_int_vetor(s_inicial_base, s_vizinho, qtd_n);

	 //obter_vizinho(vertices, C, vizinhanca, s_vizinho_encontrado, vertices_solucao, tamanho_da_rota, qtd_n);
	 obter_vizinho(vertices, C, vizinhanca, s_vizinho, vertices_solucao, tamanho_da_rota, d, prize, qtd_n);

	 if ( custo(C, s_vizinho) < custo( C, s_inicial_base) ){

	 //printf("Melhora. Vizinhanca: %d, Novo custo: %f\n", vizinhanca, custo(C, s_vizinho));

	 copy_int_vetor(s_vizinho, s_inicial_base, qtd_n);

	 if (vizinhanca == 6){
	 free(vertices_solucao);
	 tamanho_da_rota = vertices_da_solucao(vertices_solucao,	s_inicial_base);
	 }
	 //imprimir_solucao(s_inicial_base);
	 //printf("Vertices: ");
	 for (int i = 0; i < tamanho_da_rota; ++i) {
	 //printf(" %d, ",vertices_solucao[i]);
	 }

	 BL_iter = 0;

	 }

	 }

	 if (custo( C, s_inicial_base) < custo(C, s_melhor_vizinho) ){
	 copy_int_vetor(s_inicial_base, s_melhor_vizinho, qtd_n);

	 if (vizinhanca == 6){
	 free(vertices_solucao);
	 tamanho_da_rota = vertices_da_solucao(vertices_solucao,	s_melhor_vizinho);
	 }

	 vizinhanca = 0;

	 }else {
	 vizinhanca++;
	 }

	 }
	 //FIM RVND
	 */
	//Liberar Memoria
	free(vertices_solucao);
	free(s_vizinho);

}

void busca_local_v2(Vertice *vertices, Matriz* C, int *rota,
		int *s_melhor_vizinho, int d, int prize, int qtd_n) {

	copy_int_vetor(rota, s_melhor_vizinho, qtd_n);

	int melhor_custo = custo(C, s_melhor_vizinho);

	//Botar vertices da solucao em um vetor para realizar os sorteios e definir tamanho da solucao
	int* vertices_solucao;

	int tamanho_da_rota = vertices_da_solucao(vertices_solucao,
			s_melhor_vizinho);

	int ordem_extra = 0; //Ordem das vizinhancas, controla o avancar da lista de vizinhancas

	int ordem_intra = 0;

	int vizinhanca = 0; //Vizinhanca corrente

	int qtd_vizinhanca = 15;
	int qtd_vizinhanca_intra = 10;
	int qtd_vizinhanca_extra = 5;

	int aleatorio = 1; //RVND(1) ou VND(0)

	int* s_vizinho;

	s_vizinho = (int*) (malloc(qtd_n * sizeof(int)));

	copy_int_vetor(s_melhor_vizinho, s_vizinho, qtd_n);

	int vizinhancas[qtd_vizinhanca];
	int intra_rotas[qtd_vizinhanca_intra];
	int extra_rotas[qtd_vizinhanca_extra];

	int contador = 0;
	int contador_intra = 0;
	int contador_extra = 0;

	int random;

	//Aleatoriza todas
	for (int i = 0; i < qtd_vizinhanca; ++i) {
		random = (genrand_int32() % qtd_vizinhanca);
		while (in_int_vetor(random, vizinhancas, contador) == 1)
			random = (genrand_int32() % qtd_vizinhanca);
		vizinhancas[contador++] = random;
	}

	//Agora bota cada uma em seu respectivo vetor
	for (int i = 0; i < qtd_vizinhanca; ++i) {
		if (vizinhancas[i] == 5 || vizinhancas[i] == 7 || vizinhancas[i] == 9
				|| vizinhancas[i] == 12 || vizinhancas[i] == 14) {
			extra_rotas[contador_extra++] = vizinhancas[i];
		} else {
			intra_rotas[contador_intra++] = vizinhancas[i];
		}
	}

	//Explorar IntraRotas inicialmente
	while (ordem_intra < qtd_vizinhanca_intra) {

		vizinhanca = intra_rotas[ordem_intra];

		obter_vizinho(vertices, C, vizinhanca, s_vizinho, vertices_solucao,
				tamanho_da_rota, d, prize, qtd_n);

		if (custo(C, s_vizinho) < melhor_custo) {
			copy_int_vetor(s_vizinho, s_melhor_vizinho, qtd_n);

			printf(" --> (novo custo: %.0f)", custo(C, s_vizinho));

			melhor_custo = (int) custo(C, s_melhor_vizinho);

			//Para se achou custo exato
			if (melhor_custo <= valor_otimo)
				return;

			ordem_intra = 0;
		} else {
			ordem_intra++;
		}

	}
	//tem inicializar as intra-novamente.
	ordem_intra = 0;

	//BUSCA LOCAL: VND (Exploras as Extra e a cada melhora explorar as Intra)
	while (ordem_extra < qtd_vizinhanca_extra) {
		//printf("Vizinhanca: %d \n", vizinhanca);
		//printf("Maria\n");

		//Pegar primeiro melhor vizinho da vizinhanca
		if (aleatorio == 1) {
			vizinhanca = extra_rotas[ordem_extra];
		} else {
			vizinhanca = ordem_extra; //sequencialmente apartir da ZERO
		}

		obter_vizinho(vertices, C, vizinhanca, s_vizinho, vertices_solucao,
				tamanho_da_rota, d, prize, qtd_n);

		//Se vizinhanca Extra encontrou vizinho melhor
		if (custo(C, s_vizinho) < melhor_custo) {
			copy_int_vetor(s_vizinho, s_melhor_vizinho, qtd_n);

			printf(" --> (novo custo: %.0f)", custo(C, s_vizinho));

			melhor_custo = (int) custo(C, s_melhor_vizinho);

			//Para se achou custo exato
			if (melhor_custo <= valor_otimo)
				return;

			if (vizinhanca == REMOVER_INSERCAO_MAIS_BARATA
					or vizinhanca == REMOVER_INSERIR_GENIUS
					or vizinhanca == REMOVER_GENIUS_INSERIR_GENIUS
					or vizinhanca == REMOVER_GENIUS_INSERCAO_MAIS_BARATA
					or vizinhanca
							== DUPLA_REMOCAO_SIMPLES_INSERIR_MAIS_BARATA) {
				free(vertices_solucao);
				tamanho_da_rota = vertices_da_solucao(vertices_solucao,
						s_melhor_vizinho);
			}
			ordem_extra = 0;

			//Explorar IntraRotas
			while (ordem_intra < qtd_vizinhanca_intra) {

				vizinhanca = intra_rotas[ordem_intra];

				obter_vizinho(vertices, C, vizinhanca, s_vizinho,
						vertices_solucao, tamanho_da_rota, d, prize, qtd_n);

				if (custo(C, s_vizinho) < melhor_custo) {
					copy_int_vetor(s_vizinho, s_melhor_vizinho, qtd_n);

					printf(" --> (novo custo: %.0f)", custo(C, s_vizinho));

					melhor_custo = (int) custo(C, s_melhor_vizinho);

					//Para se achou custo exato
					if (melhor_custo <= valor_otimo)
						return;

					ordem_intra = 0;
				} else {
					ordem_intra++;
				}

			}
			//tem inicializar as intra-novamente.
			ordem_intra = 0;

			//Fim-Explorar IntraRotas

		} else {
			ordem_extra++;
			//vizinhanca--;
		}
		//printf("Rogerio\n");

	}

	//Liberar Memoria
	free(vertices_solucao);
	free(s_vizinho);

	limpar_solucao(s_melhor_vizinho, vertices, qtd_n, C, d, prize);

}

void gravar_dados_solucao(int *solucao, Vertice *vertices, Matriz *C,
		int qtd_n) {

}

int add_to_elite(Matriz *conjunto_elite, int *solucao, int qtd_n,
		int elite_size, int elite_element, float custo_solucao,
		float *custo_pior_elite) {

	int qtd_vertices = 0;

	//Custo da solucao Elite conjunto_elite[qtd_n]
	//Entra no conjunto elite se for melhor que a pior e for diferente das demais
	printf("Custo solucao candidata: %.0f \n", custo_solucao);
	printf("Pior Custo: %.0f \n", *custo_pior_elite);

	if (elite_element == 0)
		*custo_pior_elite = custo_solucao;

	if (elite_element < elite_size) {

		//Verifica se o custo da solução candidato é igual a de alguma elite
		for (int e = 0; e < elite_element; ++e) {
			if (get(conjunto_elite, e, qtd_n) != *custo_pior_elite
					&& get(conjunto_elite, e, qtd_n) == custo_solucao) {

				//Nao entrar no conjunto elite
				return elite_element;
			}
		}

		printf("Criando conjunto ELITE \n");
		printf("Elite Elements Inserted: %d \n", elite_element);

		for (int i = 0; i < qtd_n; ++i) {
			//Matriz elite, cada dimensao é uma solucao
			setMatriz(conjunto_elite, elite_element, i, solucao[i]);
			if (solucao[i] != -1)
				qtd_vertices++;
		}
		setMatriz(conjunto_elite, elite_element, qtd_n, custo_solucao);
		setMatriz(conjunto_elite, elite_element, qtd_n + 1, qtd_vertices);

		//Pegar pior custo
		for (int e = 0; e < elite_element; ++e) {
			if (get(conjunto_elite, e, qtd_n) > *custo_pior_elite)
				*custo_pior_elite = get(conjunto_elite, e, qtd_n);
		}
		return ++elite_element;
	} else if (custo_solucao < *custo_pior_elite) { //Melhor que pior e diferente das demais

		//Verifica se o custo da solução candidato é igual a de alguma elite
		for (int e = 0; e < elite_size; ++e) {
			if (get(conjunto_elite, e, qtd_n) != *custo_pior_elite
					&& get(conjunto_elite, e, qtd_n) == custo_solucao) {
				//Nao entrar no conjunto elite
				return elite_element;
			}
		}

		//Substituir a de pior custo
		for (int e = 0; e < elite_size; ++e) {
			if (get(conjunto_elite, e, qtd_n) == *custo_pior_elite) {

				printf("Atualizando conjunto ELITE \n");

				//Substituir
				for (int i = 0; i < qtd_n; ++i) {
					//Matriz elite, cada dimensao é uma solucao
					setMatriz(conjunto_elite, e, i, solucao[i]);
					if (solucao[i] != -1)
						qtd_vertices++;
				}
				setMatriz(conjunto_elite, e, qtd_n, custo_solucao);
				setMatriz(conjunto_elite, e, qtd_n + 1, qtd_vertices);

				//parar
				break;

			}
		} //fim substituir

		//Atualizar pior custo
		//pegar pior custo
		*custo_pior_elite = get(conjunto_elite, 0, qtd_n);
		for (int e = 0; e < elite_size; ++e) {
			printf("Custo Elite %d --> %.0f tamanho = %d \n", e,
					get(conjunto_elite, e, qtd_n),
					(int) get(conjunto_elite, e, qtd_n + 1));
			if (get(conjunto_elite, e, qtd_n) > *custo_pior_elite)
				*custo_pior_elite = get(conjunto_elite, e, qtd_n);
		}
		printf("Pior Custo %.0f \n", *custo_pior_elite);

	}

	return elite_element;

}

//int matriz_BL[MAX][MAX];
int maior = 0;
int it_maior = 0;
int antigo_cont;

float grasp_basico(vector<Vertex> vsol, int &cont, vector<edgeSol>& arestasFreq,
		int &itExecutadas, double gamma, double &pE, int **& matriz_BL,
		Vertice *vertices, Matriz *C, float alfa, int grasp_max,
		float tempo_max, int qtd_t, int prize, int seed, float d, int qtd_n,
		char *nome_instancia, int premio_solicitado, int *&s_melhor_rota,
		int numIt, int*& subProbSize) {
	antigo_cont = qtd_n * (qtd_n - 1) / 2;

	int rota[qtd_n];
	int *s_melhor_vizinho_rota = (int *) malloc(qtd_n * sizeof(int));
	float melhor_custo = 1.0E10;
	float custo_acumulado_bl = 0;
	float custo_acumulado_si = 0;
	float custo_filtro = 0;
	int contador_iteracoes = 0;
	float custo_solucao_inicial = 0;
	init_genrand(seed);

	//Mineração de Dados
	int minerou = 0;
	DMPadrao *padroes;
	int elite_size = 10;
	int elite_element = 0; //Qual elemento corrente [0..elite_size]
	float custo_pior_elite = 9999999999; //custo da pior solucao elite
	Matriz* conjunto_elite = create_matriz(elite_size, qtd_n + 2); //qtd_n(custo), nqtd_n+1(qtd_vertices)
	double cont_zero = 0.0;
	double arestas = 0.0;
	arestas = qtd_n * (qtd_n - 1) / 2;

	//Grasp_Filtro
	float filtra_em = 0; // valor percentual (exem 0.2) apartir do qual aplicar o filtro, valor 0 inativa

	//Grasp_diversifica
	//Diversificação das Soluções
	//A t iteracoes sem melhora aumentar valor do alfa em 0.1
	int diversificar = 0; //valor 1 ativa
	float alfa_inicial = alfa;
	int grasp_max_10porcento = grasp_max / 20;
	int qtd_iter_sem_melhora = 0;

	// printf("arestas = %f\n", arestas);
	// printf("qtd_n = %d\n", qtd_n);

//	if (debug)
//		printf("Valor Otimo: %d \n", valor_otimo);

	//Inicar contagem de tempo
	clock_t start = 0;
	clock_t end = 0;
	double tempo_total = 0;

	if (tempo_max > 0) {
		start = clock();
		grasp_max = 1000000000;
//		if (debug)
//			printf("Critério tempo ativado!!");
	}

	int qnd_vt = 0;
	double limPrimal = 10000.0;

	//Iteracoes GRASP
	int iteracao_final = 0;
	int parada = gamma * numIt;

	PR pr(vertices, C, 5, qtd_n); // instancia do PR para fazer uso apenas do pool

//	printf("it. \t w \t t \t prize \n");

	while (itExecutadas < numIt) {
		++itExecutadas;

#if debug
		cout << "\n\nITERAÇÃO GRASP: " << itExecutadas << "\n\n";
#endif
//		printf("%d \t", itExecutadas);

		//fase de construção
		construir_solucao_inicial(vsol, vertices, C, alfa, rota, qtd_t, prize,
				d, qtd_n, GENIUS, minerou, padroes, subProbSize);
		custo_solucao_inicial = custo(C, rota);

//		printf("%0.f\t", custo_solucao_inicial);

		//imprimir_solucao(rota, vertices);

		//atualiza melhor solucao
		if (custo_solucao_inicial < melhor_custo) {
			melhor_custo = custo_solucao_inicial;
			copy_int_vetor(rota, s_melhor_rota, qtd_n);
		}



		//fase de busca local
		busca_local(vertices, C, rota, s_melhor_vizinho_rota, d, prize, qtd_n);

		float custo_vizinho = custo(C, s_melhor_vizinho_rota);


//		printf("%0.f\t", custo_vizinho);

		//atualiza melhor solucao
		if (custo_vizinho < melhor_custo) {
			copy_int_vetor(s_melhor_vizinho_rota, s_melhor_rota, qtd_n);
			melhor_custo = custo_vizinho;

		}

		/*
		 * Atualização do E0
		 */



		int p = 0;
		while (rota[p] == -1) {
			p++;
		}

		bool E0Mudou = false;
		int antecessor = p;
		int sucessor = rota[antecessor];

		//checagem de estabilidade de E0
		if (matriz_BL[antecessor][sucessor] == 0
				or matriz_BL[sucessor][antecessor] == 0) {
			E0Mudou = true;
			matriz_BL[antecessor][sucessor] = 1;
			matriz_BL[sucessor][antecessor] = 1;
//			printf("E0 mudou \n");
			cont_zero++;
		}

		while (sucessor != p) {
			antecessor = sucessor;
			sucessor = rota[sucessor];

			//checagem de estabilidade de E0
			if (matriz_BL[antecessor][sucessor] == 0
					or matriz_BL[sucessor][antecessor] == 0) {
				E0Mudou = true;
				matriz_BL[antecessor][sucessor] = 1;
				matriz_BL[sucessor][antecessor] = 1;
				cont_zero++;
//				printf("E0 mudou \n");
			}

		}

		//contagem de iterações para estabilidade
		if (E0Mudou) {
			cont = 0; // se houve mudança em E0, o contador é resetado
		} else {
			cont++;
		}

		/*
		 * Fim da atualização
		 */

		antigo_cont = cont_zero;

//		printf("cont: %d\n", cont);
//		printf("parada em: %d\n", parada);

		if (cont == parada) { //E0 estabilizou, o grasp será encerrado

//			pr.extrairFreq(arestasFreq); //extrair as arestas mais frequentes do pool e ordenadas de acordo com a frequencia
			cont = 0;
			break;
		}

		/****************************************************
		 *  				Path-Relinking
		 ****************************************************/

		//try to update elite set with 's_melhor_vizinho_rota'
//		pr.updatePool(s_melhor_vizinho_rota, custo_vizinho);

#if debug
		printf("elite set was updated. Size: %d\n", pr.getPoolSize());
#endif

	}

////preenchendo matriz solucao
	int t = 0;
	while (s_melhor_rota[t] == -1) {
		t++;
	}

	int ant = t;
	int suc = s_melhor_rota[ant];

	edgeSol e;
	e.o = ant;
	e.d = suc;
	e.custo = get(C, ant, suc);
//	e.frequencia = pr.frequencias[ant][suc];

	arestasFreq.push_back(e);

//	sol[ant][suc] = 1;
//	sol[suc][ant] = 1;
	while (suc != t) {
		ant = suc;
		suc = s_melhor_rota[suc];
//		sol[ant][suc] = 1;
//		sol[suc][ant] = 1;
		e.o = ant;
		e.d = suc;
		e.custo = get(C, ant, suc);
//		e.frequencia = pr.frequencias[ant][suc];
		arestasFreq.push_back(e);

	}

//	printf("tam da sol do grasp: %d \n", arestasFreq.size());

//	printf("iterações: %d \n", iteracao);

//	printf("%.0f \n", cont_zero);

//	if (debug)
//		printf("*****************************\n");
//	if (debug)
//		printf("***** GRASP FINALIZADO *****\n");
//	if (debug)
//		printf("-----------------------------------------\n\n");

//	imprimir_solucao(s_melhor_rota, vertices);
//
//	if (debug)
//		printf("* Melhor Solucao * \n\n");
//	if (debug)
//		imprimir_solucao(s_melhor_rota, vertices);
//	printf("\nPRIZE Coletado / Solicitado: %d / %d\n",
//			prize_acumulado(vertices, s_melhor_rota, qtd_n), prize);

//	if (avaliar_solucao(vertices, C, s_melhor_rota, prize, d, qtd_n) != -1) {
//		if (debug)
//			printf("Solução Válida !\n");
//	}
//

//	free(s_melhor_vizinho_rota);
//	free_matriz(conjunto_elite);

//	if (minerou == 1)
//		free(padroes);

	return melhor_custo;
}

void print_matriz(int m[MAX][MAX], int n) {
	float cont_zero;
	float arestas;
	arestas = n * (n - 1) / 2;

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			printf("[%d][%d] = %d \t", i, j, m[i][j]);
			if (m[i][j] == 0)
				cont_zero++;
		}
		printf("\n");
	}

	printf("zeros = %.2f \n", cont_zero / arestas);
}

void preenche_matriz(int m[MAX][MAX], int* rota) {
	int p = 0;
	while (rota[p] == -1) {
		p++;
	}

//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = rota[antecessor];
	m[antecessor][sucessor] = 1;
	m[sucessor][antecessor] = 1;
	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = rota[sucessor];
		m[antecessor][sucessor] = 1;
		m[sucessor][antecessor] = 1;
	}
}

//string solutionToString(int* rota) {
//	string str;
//	int p = 0;
//	while (rota[p] == -1) {
//		p++;
//	}

////p e primeiro vertice pertencente a solucao
//	int antecessor = p;
//	int sucessor = rota[antecessor];
//	str += intToString(antecessor) + " ";
//	str += intToString(sucessor) + " ";
//	while (sucessor != p) {
//		antecessor = sucessor;
////		str += intToString(antecessor) + ",";
//		sucessor = rota[sucessor];
//		if (sucessor == p)
//			str += intToString(sucessor);
//		else
//			str += intToString(sucessor) + " ";
//	}
//	return str;
//}

//string intToString(int x) {
//	stringstream ss;
//	ss << x;
////	string s;
////	s.compare()
//	return ss.str();
//}
