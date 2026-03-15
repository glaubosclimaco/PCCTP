/*
 * mineracao.h
 *
 *  Created on: 02/04/2013
 *      Author: rogerio
 */

#ifndef MINERACAO_H_
#define MINERACAO_H_

#include "Utils.h"

//Componente Conexa CC
typedef struct componente_conexa{
	int *aresta_id; //ids das arestas da cc no padrão
	int tamanho;  //qtd de arestas
}ComponenteConexa;

//Aresta: representa uma aresta em padrão minerado
typedef struct Aresta{
	int id; //identifica unicamente uma aresta
	int vertice_origem;
	int vertice_destino;
	int analisada; // usado no processo de identificacao de CC [ 0 - nao; 1 - sim ]
}Aresta;

//Padrao
typedef struct padrao{
	int id;
	int tamanho; //qtd de areastas
	int suporte;
	Aresta *arestas;
	ComponenteConexa *ccs; //conjunto de componentes conexas (CC)
	int qtd_cc;
	int cc_atual; //qual CC usar no processo de solucao inicial
}DMPadrao;

//Tipo Estruturado para Matriz e funcoes relacionadas
//typedef struct matriz{
//	int linha;
//	int coluna;
//	float* v;
//}Matriz;

Matriz* create_matriz(int linhas, int colunas);
float get(Matriz* mat, int i, int j);
void set(Matriz* mat, int i, int j, float valor);
void free_matriz(Matriz* mat);

void encode(Matriz* conjunto_elite, int elite_size,int qtd_n);
int decode(int qtd_n, DMPadrao *&padroes);
int mining(Matriz* conjunto_elite, int elite_size,int qtd_n, DMPadrao *&padroes);

//Verificar se vertice_v é vertice de origem no padrão passado
Aresta in_padrao(int vertice_v, DMPadrao padrao);
Aresta in_padrao2(int vertice_v, DMPadrao padrao);

int in_componente_conexa(int vertice_v, DMPadrao padrao, ComponenteConexa &cc);

//Se esta em outra cc diferente de cc
int in_componente_conexa2(int vertice_v, DMPadrao padrao, int cc_index);


#endif /* MINERACAO_H_ */
