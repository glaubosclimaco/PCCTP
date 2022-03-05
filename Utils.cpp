/*
 * Utils.cpp
 *
 *	Neste arquivo devera ter rotinas auxiliares
 *
 *	1) Carregar instancias
 *
 *  Created on: 20/05/2012
 *      Author: Rogerio
 */

#include "Utils.h"

FILE* abrir_arquivo(char *file_name) {

	FILE *fp;

	fp = fopen(file_name, "rt");

	if (fp == NULL) {
		printf("Erro na abertura do arquivo!\n");
	} else {
		printf("Arquivo aberto com sucesso!\n");
	}

	return (fp);

}

//Funcao para comparar dois numero reais
int comp_int(const void* p1, const void* p2) {

	//convert ponteiros genericos em pontos de int
	int *f1 = (int*) p1;
	int *f2 = (int*) p2;

	//faz a comparados dos valores que os ponteiros apontam
	if (*f1 < *f2)
		return -1;
	else if (*f1 > *f2) {
		return 1;
	} else {
		return 0;
	}

}

int in_int_vetor(int a, int *v, int tam) {
	int achou = 0;

	if (!tam) {
		return achou;
	}

	for (int i = 0; i < tam; i++) {
		if (a == v[i]) {
			//achou = 1;
			return 1;
		}
	}
	return achou;
}

int copy_int_vetor(int *vetor_from, int *vetor_to, int size) {

	for (int i = 0; i < size; ++i) {
		vetor_to[i] = vetor_from[i];
	}

	return 0;

}

int get_indice_do_menor_valor(int *vetor, int size, int *indices_excluidos,
		int size_excluidos) {

	int menor = -1;
	int indice = -1;
	int cont = 0;
	while (menor == -1) {
		if ((in_int_vetor(cont, indices_excluidos, size_excluidos) == 0)
				&& (vetor[cont] > 0)) {
			menor = vetor[cont];
			indice = cont;
		}
		cont++;
	}

	for (int i = 0; i < size; i++) {
		if (in_int_vetor(i, indices_excluidos, size_excluidos) == 0) {
			if (vetor[i] < menor && vetor[i] > 0) {
				menor = vetor[i];
				indice = i;
			}
		}
	}

	return indice;
}

int get_indice_do_maior_valor(int *vetor, int size, int *indices_excluidos,
		int size_excluidos) {

	int maior = -1;
	int indice = -1;
	int cont = 0;
	while (maior == -1) {
		if (in_int_vetor(cont, indices_excluidos, size_excluidos) == 0) {
			maior = vetor[cont];
			indice = cont;
		}
		cont++;
	}

	for (int i = 0; i < size; i++) {
		if (in_int_vetor(i, indices_excluidos, size_excluidos) == 0) {
			if ((vetor[i] > maior)) {
				maior = vetor[i];
				indice = i;
			}
		}
	}

	return indice;
}

Matriz* create_matriz(int linhas, int colunas) {
	Matriz* mat = (Matriz*) malloc(sizeof(Matriz));
	mat->linha = linhas;
	mat->coluna = colunas;
	mat->v = (float*) malloc(linhas * colunas * sizeof(float));
	/*for (int i = 0; i < linhas; i++) {
	 mat->v[i] = (float*) malloc(colunas*sizeof(float));
	 }*/
	return mat;
}

float get(Matriz* mat, int i, int j) {
	/*
	 if (i < 0 || i >= mat->linha || j < 0 || j >= mat->coluna ){
	 printf("Linhas: %d e Colunas: %d\n", mat->linha, mat->coluna);
	 printf("i: %d e j: %d \n",i,j);
	 printf("Acesso inválido em Matriz! ");
	 exit(1);
	 }
	 */
	//return mat->v[i][j];
	int k;
	k = i * mat->coluna + j;
	return mat->v[k];
}

void setMatriz(Matriz* mat, int i, int j, float valor) {

	if (i < 0 || i >= mat->linha || j < 0 || j >= mat->coluna) {
		printf("Acesso inválido em Matriz!");
	}

	//mat->v[i][j] = valor;
	int k;
	k = i * mat->coluna + j;
	mat->v[k] = valor;

}

void free_matriz(Matriz* mat) {
	/*
	 printf("a\n");
	 for (int i = 0; i < mat->linha; ++i){
	 printf("b[%d]\n",i);
	 free(mat->v[i]);
	 }
	 printf("c\n");
	 free(mat->v);
	 printf("d\n");
	 //free(mat);
	 printf("e\n");
	 */
	free(mat->v);
	free(mat);
}

int qtd_vertices_em_s(int *s, int qtd_n) {
	int qtd = 0;

	for (int i = 0; i < qtd_n; ++i) {
		if (s[i] != -1)
			qtd++;
	}

	return qtd;
}

int qtd_vertices_em_lc(int *LC, int qtd_lc) {
	int qtd = 0;

	for (int i = 0; i < qtd_lc; ++i) {
		if (LC[i] == 1)
			qtd++;
	}

	return qtd;
}

int vertices_de_LC(int *LC, int * &vertices_lc, int qtd_n) {
	//vertices em LC
	int vertices[qtd_n];

	int qtd_lc = 0;

	for (int i = 0; i < qtd_n; i++) {

		if (LC[i] == 1) {
			vertices[qtd_lc] = i;
			qtd_lc++;
		}

	}

//	printf("Tamanho de LC: %d, V em LC: ", qtd_lc);

	//Criar o vetor no tamanho certo;
	free(vertices_lc);

	vertices_lc = (int *) malloc(qtd_lc * sizeof(int));

	//botar os vertices em *vertices_lc
	for (int i = 0; i < qtd_lc; ++i) {
		vertices_lc[i] = vertices[i];
		//	printf("%d, ", vertices_lc[i]);
	}

//	printf("\n");

	return qtd_lc;
}

int get_premio_vertice(int vertice, Vertice *vertices_tsp) {

	int parou = 0;
	int i = 0;
	int premio = 0;

	while (parou == 0) {
		if (vertices_tsp[i].id == vertice) {
			premio = vertices_tsp[i].premio;
			parou = 1;
		} else {
			i++;
		}
	}

	return premio;

}

/*
 void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
 {
 long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
 struct rusage ptempo;
 //TODO

 getrusage(0,&ptempo);

 seg_CPU = ptempo.ru_utime.tv_sec;
 mseg_CPU = ptempo.ru_utime.tv_usec;
 seg_sistema = ptempo.ru_stime.tv_sec;
 mseg_sistema = ptempo.ru_stime.tv_usec;

 *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
 *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
 }
 */

int antecessor_de(int *s, int vertice) {
	int i = 0;
	int achou = 0;

	while (achou == 0) {
		if (s[i] == vertice) {
			achou = 1;
		} else {
			i++;
		}
	}

	return i;
}

int comp_vizinhos(const void* p1, const void* p2) {

	//convert ponteiros generics em pontos de Vizinho
	Vizinho *f1 = (Vizinho*) p1;
	Vizinho *f2 = (Vizinho*) p2;

	float d1 = f1->distancia;
	float d2 = f2->distancia;

	//faz a comparados dos valores que os ponteiros apontam
	if (d1 < d2)
		return -1;
	else if (d1 > d2) {
		return 1;
	} else {
		return 0;
	}

}

float get_custo_incremental_genius_v2(int *rota, int qtd_n, int vertice,
		Matriz *C) {

	//Menos somatório de custos entre os p mais próximos
	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);
	float menor_custo = 100000000000000;
	float custo = 0;
	int vertice_k = 0;

	int p = 5;
	if (p > tamanho_da_rota) {
		p = tamanho_da_rota;
	}

	NP* NP_S = (NP *) (malloc(qtd_n * sizeof(NP)));
	NP_S[vertice].vizinhos = (int *) malloc(p * sizeof(int));
	int * vertices_k;
	int qtd_k = 0;

	pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C, NP_S,
			vertice, p);

	//Vértice i e j, pertencem a NP deo vertice
	//Vértice k, pertece ao caminho j até i e pertence a NP de i+1
	//Para cada combinação de i e j
	for (int i = 0; i < p; i++) {
		int vertice_i = NP_S[vertice].vizinhos[i];
		//Gerar lista NP para sucessor de NP
		NP_S[rota[vertice_i]].vizinhos = (int *) malloc(p * sizeof(int));
		pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C, NP_S,
				rota[vertice_i], p);

		for (int j = i + 1; j < p; ++j) {
			int vertice_j = NP_S[vertice].vizinhos[j];

//			printf("\nVertice i[%d] = %d; vertice j[%d] = %d", i,vertice_i, j,vertice_j);

			//Se adjacentes (custo: Civ + Cvj - Cij
			if (rota[vertice_i] == vertice_j or rota[vertice_j] == vertice_i) {
				//printf("\nAdjcentes");
				custo = get(C, NP_S[vertice].vizinhos[i], vertice)
						+ get(C, vertice, NP_S[vertice].vizinhos[j])
						- get(C, NP_S[vertice].vizinhos[i],
								NP_S[vertice].vizinhos[j]);

				//printf("--> Custo: %.0f <--", custo);

				if (custo < menor_custo)
					menor_custo = custo;

			} else { //se não for adjacentes ver lógica do GENIUS

				//pegar vértice K
				qtd_k = obter_vertices_k(rota, vertice_i, vertice_j, NP_S,
						vertices_k, p);

				if (qtd_k == 0) { //retornar vertices que tenhoa os vizinhos mais proximos
					//printf("\nK - NOT OK");

					continue;

				} else { //Ver logica do GENIUS
					//printf("\nK - OK (%d)", qtd_k);

					for (int k = 0; k < qtd_k; ++k) {

						vertice_k = vertices_k[k];
						if (vertice_i == vertice_k or vertice_j == vertice_k
								or rota[vertice_k] == rota[vertice_j])
							continue;

						float custo_in = get(C, vertice_i, vertice)
								+ get(C, vertice, vertice_j)
								+ get(C, rota[vertice_i], vertice_k)
								+ get(C, rota[vertice_k], rota[vertice_j]);
						float custo_out = get(C, vertice_i, rota[vertice_i])
								+ get(C, vertice_j, rota[vertice_j])
								+ get(C, vertice_k, rota[vertice_k]);

						custo = custo_in - custo_out;

						if (custo < menor_custo)
							menor_custo = custo;

						// printf("--> Custo: %.0f <--", custo);
					}

					free(vertices_k);

				}

			}

		}
		//A cada iteracao liberar os vizinhos do i+1
		free(NP_S[rota[vertice_i]].vizinhos);
	}

	free(NP_S[vertice].vizinhos);
	free(NP_S);

	free(vertices_solucao);

	return menor_custo;

}

float get_custo_incremental_genius(int *rota, int qtd_n, int vertice,
		Matriz *C) {

	//Menos somatório de custos entre os p mais próximos
	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

	int p = 5;
	if (p > tamanho_da_rota) {
		p = tamanho_da_rota;
	}

	NP* NP_S = (NP *) (malloc(qtd_n * sizeof(NP)));
	NP_S[vertice].vizinhos = (int *) malloc(p * sizeof(int));

	pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C, NP_S,
			vertice, p);

	float distancia_2_mais_proximos = (get(C, vertice,
			NP_S[vertice].vizinhos[0])
			+ get(C, vertice, NP_S[vertice].vizinhos[1]));

	free(NP_S[vertice].vizinhos);
	free(NP_S);

	free(vertices_solucao);

	return distancia_2_mais_proximos;
}

int criar_LRC_4_genius(Vertice *vertices, int *vertices_lc, int qtd_lc,
		int *rota, int qtd_n, Matriz *C, float alfa) {

	//aproveita a estrutura vizinho para ordenar os curso incrementais e gerar LRC
	Vizinho candidatos[qtd_lc];
	int tam_lrc = 0;

	//Primeiro passo, benefício de cada vértice candidato

	//printf("Criar LRC!!\n");
	for (int i = 0; i < qtd_lc; ++i) {
		candidatos[i].vertice = vertices_lc[i];
		candidatos[i].distancia = get_custo_incremental_genius_v2(rota, qtd_n,
				vertices_lc[i], C);
		//*vertices[vertices_lc[i]].premio;
	}

	//Ordenar as distancias (digo custos incrementais indiretos)
	qsort(candidatos, qtd_lc, sizeof(Vizinho), comp_vizinhos);

	float min = candidatos[0].distancia;
	float max = candidatos[qtd_lc - 1].distancia;

	float custo_de_corte = (min + (alfa * (max - min)));

//	printf("qtd_lc: %d", qtd_lc);
//	printf("Min: %.0f \n", min);
//	printf("Max: %.0f \n", max);
//	printf("Alfa: %.1f \n", alfa);
//	printf("LRC (Corte: %.0f): \n", custo_de_corte);

	for (int i = 0; i < qtd_lc; i++) {

		if (candidatos[i].distancia <= custo_de_corte) {
			//printf("%d, ", candidatos[i].vertice);
			tam_lrc++;
		}

	}
	//printf("\n");

	//Sortear

//	printf("tam_lrc: %d \n", tam_lrc);

	int pos_lcr = (genrand_int32() % tam_lrc);

	//printf(" Sorteado: %d (%d / %d)\n", candidatos[pos_lcr].vertice, tam_lrc, qtd_lc);

	return candidatos[pos_lcr].vertice;
}

void criar_LC_4_genius(Vertice *vertices, int *rota, int *LC, int qtd_n,
		int somente_t) {
	/* LC:
	 * Vetor Binário, onde o indice identifica o vertice e valor é binário,
	 * onde 1 indica vertice candidato
	 * a,  b e c sao os tres vertices iniciais de s_inicial, escolhidos aleatoriamente entre os T
	 */

	int *vertices_w;

	int qtd_wd = 0;

	if (somente_t == 0) {
		qtd_wd = vertices_w_descobertos(vertices, rota, qtd_n, vertices_w);
	}

	for (int i = 0; i < qtd_n; ++i) {

		//Avanca se i ja esta na rota
		//if (rota[i] != -1) continue;

		LC[i] = 0;

		//Qq T sempre é candidato
		if ((strcmp(vertices[i].tipo, "t") == 0) and rota[i] == -1) {
			LC[i] = 1;
		}

		if (somente_t == 0) {

			if ((strcmp(vertices[i].tipo, "v/t") == 0) and (rota[i] == -1)) {
				//Se houver w descobertos, inserir como candidatos somente vertices V/T que os cubra
				if (qtd_wd > 0) {
					LC[i] = 0;
					for (int w = 0; w < qtd_wd; ++w) {
						if (in_int_vetor(vertices_w[w], vertices[i].cobre,
								vertices[i].qtd_cobre) == 1) {
							LC[i] = 1;
							break;
						}
					}
				} else { // Se nao houver w descobertos, todos os V/T serão candidatos
					LC[i] = 1;
				}

			}
		}
	}
	/*
	 LC[a] = 0;
	 LC[b] = 0;
	 LC[c] = 0;
	 */

	if (qtd_wd > 0)
		free(vertices_w);

}

int vertices_da_solucao(int*& vertices_solucao, int* &solucao) {
	vertices_solucao = (int*) (malloc(2 * sizeof(int)));
	int p = 0;
	while (solucao[p] == -1)
		p++;
	//p eh o primeiro vertice pertencente a solucao
	int ant = p;
	int pos = solucao[ant];
	vertices_solucao[0] = ant;
	vertices_solucao[1] = pos;
	int tamanho_da_rota = 2;
	while (pos != p) {
		ant = pos;
		pos = solucao[pos];
		if (pos != p) {
			tamanho_da_rota++;
			vertices_solucao = (int*) (realloc(vertices_solucao,
					(tamanho_da_rota) * sizeof(int)));
			vertices_solucao[tamanho_da_rota - 1] = pos;
		}
	}
	return tamanho_da_rota;
}

//Funcao para comprar distancia entre vertices
int comp_inteiro(const void* p1, const void* p2) {

	//convert ponteiros generics em pontos de int
	int *f1 = (int *) p1;
	int *f2 = (int *) p2;

	//faz a comparacao dos valores que os ponteiros apontam
	if (*f1 < *f2)
		return -1;
	else if (*f1 > *f2) {
		return 1;
	} else {
		return 0;
	}

}

//Funcao para comprar dois numero reais
int comp_reais(const void* p1, const void* p2) {

	//convert ponteiros generics em pontos de float
	float *f1 = (float*) p1;
	float *f2 = (float*) p2;

	//faz a comparados dos valores que os ponteiros apontam
	if (*f1 < *f2)
		return -1;
	else if (*f1 > *f2) {
		return 1;
	} else {
		return 0;
	}

}

//Funcao para comprar dois vizinhos
float custo(Matriz* C, int *solucao) {

	int primeiro_vertice = 0;
	while (solucao[primeiro_vertice] == -1)
		primeiro_vertice++;

	int vertice_ant = primeiro_vertice;
	int vertice_pos = solucao[vertice_ant];

	float sum_distancia = get(C, vertice_ant, vertice_pos);

	while (vertice_pos != primeiro_vertice) {

		vertice_ant = vertice_pos;
		vertice_pos = solucao[vertice_pos];
		sum_distancia += get(C, vertice_ant, vertice_pos);

	}

	return sum_distancia;

}

float pegar_tempo_matheus(char nome_instancia[30], int seed) {

	char instancia[30];
	char ctoken[30];
	int semente = -1;
	float tempo;

	std::ifstream get("tempos_matheus.txt");

	for (int i = 0; i < 1440; ++i) {

		get >> ctoken;

		strcpy(instancia, ctoken);

		get >> ctoken;

		semente = atoi(ctoken);

		get >> ctoken;

		tempo = atof(ctoken);

		/*
		 printf("\nInstancia: %s", instancia);
		 printf("\tSemente: %d", semente);
		 printf("\ttempo: %f", tempo);
		 */
		if ((strcmp(nome_instancia, instancia) == 0) and (seed == semente))
			return tempo;
	}

	return 0;
}

int prize_acumulado(Vertice *vertices, int *rota, int qtd_n) {
	int prize = 0;
	/*
	 for (int i = 0; i < qtd_n; i++) {
	 if ((strcmp(vertices[i].tipo,"w") != 0) && (in_int_vetor(i, s_inicial, qtd_n) == 1)){
	 prize += vertices[i].premio;
	 }
	 }*/

	for (int i = 0; i < qtd_n; ++i) {
		if (rota[i] != -1)
			prize += vertices[i].premio;
	}

	return prize;
}

int existe_t_fora_da_solucao(Vertice *vertices, int *rota, int qtd_n) {

	for (int i = 0; i < qtd_n; i++) {
		//if ((strcmp(vertices[i].tipo,"t") == 0) && (in_int_vetor(i, s_inicial, qtd_n) == 0) ){ // Se Obrigatorio T.
		if ((strcmp(vertices[i].tipo, "t") == 0) && (rota[i] == -1)) { // Se Obrigatorio T.
			//printf("T fora da solucao: %d\n", i);
			return 1;
		}
	}

	return 0;
}

int existe_w_descoberto(Vertice *vertices, Matriz* C, float d, int *rota,
		int qtd_n) {
	/*
	 * Para qualquer i pertecente W, deve haver pelo menos um j Pertecente a V,
	 * cujo Cij <= D.
	 */

	for (int i = 0; i < qtd_n; ++i) {
		int descoberto = 1;
		if (strcmp(vertices[i].tipo, "w") == 0 and vertices[i].reduzido == 0) { //se for W e nao tiver sido reduzido
			for (int j = 0; j < qtd_n; j++) {
				if ((strcmp(vertices[j].tipo, "t") == 0)
						or strcmp(vertices[j].tipo, "v/t") == 0) { //se for V
						//if ((get(C, i, j) <= d) && (in_int_vetor(j, rota, qtd_n)==1)){ // j E V cobre i E j está na solução
					if ((get(C, i, j) <= d) && (rota[j] != -1)) { // j E V cobre i, j está na solução
						descoberto = 0;
					}
				}
			}
			if (descoberto == 1) {
				////printf("Vertice W=%d esta descoberto\n", i);
				return 1;
			}
		}
	}

	////printf("Todos os vertices W ja cobertos\n");
	return 0;
}

int checarSolucao(Vertice *vertices, Matriz *C, int *rota, int prize, float d,
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
		return 0;
	}

	// printf("PREMIO OK!\n");

	if (existe_t_fora_da_solucao(vertices, rota, qtd_n) != 0) {
		printf("Restricao de Vertices Obrigatorios Violada! \n");
		return 0;
	}

	// printf("VERTICES T OK! \n");

	if (existe_w_descoberto(vertices, C, d, rota, qtd_n) != 0) {
		printf("Restricao de Vertices de Cobertura Violada! \n");
		return -1;
	}

	for (int i = 0; i < qtd_n; ++i) {
		if (strcmp(vertices[i].tipo, "w") == 0 and rota[i] != -1) {
			printf("Vértice W na rota! \n");
			return 0;
		}

	}

	// printf("VERTICES W OK! \n");

	printf("SOLUCAO VIAVEL! \n");
	return 1;
}
