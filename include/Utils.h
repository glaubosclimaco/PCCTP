/*
 * Ultis.h
 *
 *  Created on: 31/08/2017
 *      Author: user
 */

#ifndef ULTIS_H_
#define ULTIS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
//#include <sys/resource.h> //Nao existe no windows :(
#include "random_provider.h"
#include <fstream>
#include <limits.h>


typedef struct np {
	int vertice; //identificador do vertices
	int *vizinhos; //Vizinhos mais proximos
} NP;

typedef struct matriz {
	int linha;
	int coluna;
	float* v;
} Matriz;

//para analise
//int contw = 0, contt = 0, contp = 0;
//

FILE* abrir_arquivo(char* file_name);

int comp_int(const void* p1, const void* p2);

int in_int_vetor(int a, int *v, int tam);

int vertices_da_solucao(int*& vertices_solucao, int*& solucao);

int obter_vertices_k(int *rota, int vertice_i, int vertice_j, NP *NP_lista,
		int *&vertices_k, int p);

//Tipo Estruturado para vertice
typedef struct vertice {
	int id; //precisei..entao botei., guarda o numero do vertices (0, 1, 2, ... N)
	int x;
	int y;
	const char tipo[4];
	int premio;
	int *cobre; //conjunto de vertices W, que este vertice V cobre
	int *coberto_por; //conjunto de vertices V, que cobre o este W
	int qtd_cobre;
	int qtd_coberto_por;
	int reduzido;
} Vertice;

typedef struct edgeSol {
	int o, d;
	double custo = 0.0;
	double premio = 0.0; //premio acumulado dos vértices da aresta
	unsigned frequencia = 0;
} edgeSol;

//Tipo usado para gerar NP
typedef struct vizinho {
	int vertice;
	float distancia;
} Vizinho;

void pegar_vertices_proximos_v2(int *vertices_solucao, int qtd_vertices_solucao,
		Matriz *C, NP *NP_, int vertice, int p);
int vertices_w_descobertos(Vertice *vertices, int *rota, int qtd_n,
		int *&vertices_w);
int copy_int_vetor(int *vetor_from, int *vetor_to, int size);
void setMatriz(Matriz* mat, int i, int j, float valor);
float custo(Matriz* C, int *solucao);
int checarSolucao(Vertice *vertices, Matriz *C, int *rota, int prize, float d,
		int qtd_n);
int prize_acumulado(Vertice *vertices, int *rota, int qtd_n);
int existe_t_fora_da_solucao(Vertice *vertices, int *rota, int qtd_n);
int existe_w_descoberto(Vertice *vertices, Matriz* C, float d, int *rota,
		int qtd_n);


#endif /* ULTIS_H_ */
