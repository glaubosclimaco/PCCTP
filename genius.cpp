/*
 * genius.cpp
 *
 * Implementacao da Heuristica GENIUS
 *
 *  Created on: 26/06/2012
 *      Author: rogerio
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Utils.h"
#include <vector>
#include <algorithm>
#include "./model/Graph.h"

//Numero Aleatorios
#include "random_provider.h"

#include "mineracao.h"

//Tipo Estruturado para vértice
//typedef struct vertice{
//	int id; //precisei..entao botei., guarda o numero do vertices (0, 1, 2, ... N)
//	int x;
//	int y;
//	char tipo [4];
//	int premio;
//	int *cobre; //conjunto de vertices W, que este vertice V cobre
//	int *coberto_por; //conjunto de vertices V, que cobre o este W
//	int qtd_cobre;
//	int qtd_coberto_por;
//	int reduzido;
//}Vertice;

////Tipo estruturado para os Conjuntos NP de genius
//typedef struct np{
//	int vertice; //identificador do vertices
//	int *vizinhos; //Vizinhos mais proximos
//}NP;
//
////Tipo usado para gerar NP
//typedef struct vizinho{
//	int vertice;
//	float distancia;
//}Vizinho;

//from grasp
int existe_w_descoberto(Vertice *vertices, Matriz* C, float d, int *s_inicial,
		int qtd_n);
int existe_t_fora_da_solucao(Vertice *vertices, int *s_inicial, int qtd_n);
int prize_acumulado(Vertice *vertices, int *s_inicial, int qtd_n);
void imprimir_solucao(int *solucao, Vertice* vertices);
void inverter_sentido_rota(int*& s_vizinho_encontrado, int vertice_b,
		int vertice_c);
void inserir_k_na_solucao(int *&solucao, int k, int pos);

//from utils
int vertices_de_LC(int *LC, int *&vertices_lc, int qtd_n);
void criar_LC_4_genius(Vertice *vertices, int *rota, int *LC, int qtd_n,
		int somente_t);
int criar_LRC_4_genius(Vertice *vertices, int *vertices_lc, int qtd_lc,
		int *rota, int qtd_n, Matriz *C, float alfa);
int vertices_da_solucao(int*& vertices_solucao, int*& solucao);
int comp_reais(const void* p1, const void* p2);
int in_int_vetor(int a, int *v, int tam);
float custo(Matriz* C, int *solucao);
int copy_int_vetor(int *vetor_from, int *vetor_to, int size);
int antecessor_de(int *s, int vertice);
void remover_k_da_solucao(int *&rota, int k);
void inserir_k_na_rota_genius_like(Vertice* vertices, int *&rota, Matriz *C,
		int vertice_v, int vertice_i, int p, int qtd_n);
int comp_vizinhos(const void* p1, const void* p2);

//Sentidos de Insercao em GENIUS
enum SENTIDO_INSERCAO {
	HORARIO = 0, ANTI_HORARIO = 1,
};

void pegar_vertices_proximos_v2(int *vertices_solucao, int qtd_vertices_solucao,
		Matriz *C, NP *NP_, int vertice, int p) {
	//Para cada vértice ordenar a distância dele para todos os demais em seguida pegar o p primeiros que estejam na solucao

	Vizinho vizinhos[qtd_vertices_solucao];

	for (int i = 0; i < qtd_vertices_solucao; i++) {
		vizinhos[i].vertice = vertices_solucao[i];
		vizinhos[i].distancia = get(C, vertice, vertices_solucao[i]);
	}

	//Ordenar as distancias
	qsort(vizinhos, qtd_vertices_solucao, sizeof(Vizinho), comp_vizinhos);

	int inicio = 0;
	if (vizinhos[0].vertice == vertice)
		inicio = 1;

	for (int i = 0; i < p; ++i) {
		NP_[vertice].vizinhos[i] = vizinhos[i + inicio].vertice;
	}

}

void pegar_vertices_proximos(int *vertices_solucao, int qtd_vertices_solucao,
		Matriz *C, NP *NP_, int vertice, int p) {
	//Pegar o P vertices em S, mais proximos de V

	/* Estrategia inicial: Ordenar as distancias dos vertices para v e depois pegar a P distancia e adicionar todos que estiver a P ou menos de distancia
	 *
	 * Sugestao: Logo apos carregar a instancia, para vertice V, ordenar todos demais vertices e botar em vetor para facilitar essa tarefa
	 *
	 * */

	////printf("Pegando vizinhos de (%d)..\n", vertice);
	/*
	 //printf("tamanho da rota: %d\n", qtd_vertices_solucao);
	 for (int i = 0; i < qtd_vertices_solucao; ++i) {
	 //printf(" --> %d", vertices_solucao[i]);
	 }
	 */

	int v = vertice;
	int real_p = p; //qtd real de vizinho a ser pego
	int posicao_x = -1; //Posicao da distancia limite, ou seja, os p_vizinhos serao que tem distancia menor ou igual a da posicao_x em distancias ordenadas
	////printf("\n");
	if (in_int_vetor(v, vertices_solucao, qtd_vertices_solucao) == 0) {
		////printf("Vertice Candidato!\n");
		if (qtd_vertices_solucao < p) {
			////printf("c\n");
			real_p = qtd_vertices_solucao;
		}
	} else {
		////printf("Vertice (%d) ta na rota \n", v);
		if (qtd_vertices_solucao <= p) {
			////printf("P == a qtd de vertices da solucao \n");
			real_p = qtd_vertices_solucao - 1;
		}
	}
	////printf("Real P: %d\n", real_p);
	posicao_x = real_p - 1;
	int qtd_distancias = (
			in_int_vetor(v, vertices_solucao, qtd_vertices_solucao) == 1 ?
					qtd_vertices_solucao - 1 : qtd_vertices_solucao);

	float distancias[qtd_distancias];

	//Pegar as distancias de todos os vertices da solucao para V
	////printf("Distancias (qtd=%d): ", qtd_distancias);
	int j = 0;
	for (int i = 0; i < qtd_vertices_solucao; ++i) {
		if (v != vertices_solucao[i]) {
			distancias[j] = get(C, v, vertices_solucao[i]);
			////printf("V(%d) %.1f",vertices_solucao[i], distancias[j]);
			j++;
		}
	}

	//Ordenar as distancias
	qsort(distancias, qtd_distancias, sizeof(float), comp_reais);

	////printf("\nDistancias ordenadas(qtd=%d): ", qtd_distancias);
	/*
	 for (int i = 0; i < qtd_distancias; ++i) {
	 //printf(" --> %.1f", distancias[i]);
	 }
	 //printf("\n");
	 */
	//Pegar os P vertices que esteja ate a P distancia
	int qtd_inseridos = 0;
	/*
	 //printf("Qtd de vizinhos: %d\n", real_p);
	 //printf("Distancia limite: %.1f\n", distancias[posicao_x]);
	 */
	//Botar em NPS somente os vertices que esteja com ate P distancia
	for (int i = 0; i < qtd_vertices_solucao; ++i) {
		if (v != vertices_solucao[i]) {
			if (get(C, v, vertices_solucao[i]) <= distancias[posicao_x]) {
				NP_[v].vizinhos[qtd_inseridos] = vertices_solucao[i];
				////printf("add vizinho: %d ", vertices_solucao[i]);
				qtd_inseridos++;
				if (qtd_inseridos == real_p) {
					break;
				}
			}
		}
	}

}

void gerar_NP_lista(int * s_inicial, Matriz *C, NP* &NP_S, int p, int qtd_n,
		int *vertices_solucao, int tamanho_da_rota) {

	////printf("Gerando NP (%d): \n", qtd_n);

	for (int i = 0; i < qtd_n; i++) {
		////printf("\n NP para vertice: %d --> ", i);
		NP_S[i].vertice = i;
		NP_S[i].vizinhos = (int *) malloc(p * sizeof(int));
		/* Pegar os P vertices em S mais proximo de NP_S[i].vertice
		 * e armazenar em NP_S[i].vizinhos
		 */

		pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C, NP_S,
				i, p);

	}
	////printf("NP gerada!\n");

}

void imprimir_np_lista(int qtd_n, int p, int* vertices_solucao,
		int tamanho_da_rota, NP* NP_lista) {

	for (int i = 0; i < qtd_n; ++i) {
		int real_p = p;
		if (in_int_vetor(i, vertices_solucao, tamanho_da_rota) == 0) {
			//printf("Vertice Candidato!\n");
			if (tamanho_da_rota < p) {
				//printf("c\n");
				real_p = tamanho_da_rota;
			}
		} else {
			//printf("Vertice (%d) ta na rota \n", v);
			if (tamanho_da_rota <= p) {
				//printf("P == a qtd de vertices da solucao \n");
				real_p = tamanho_da_rota - 1;
			}
		}

		printf("Vizinhos de %d: ", i);
		for (int j = 0; j < real_p; ++j) {
			printf("%d, ", NP_lista[i].vizinhos[j]);
		}
		printf("\n");
	}

}

int vertices_no_caminho(int *rota, int vertice_i, int vertice_j,
		int *&vertices_no_caminho) {

	int vertice_atual = rota[vertice_i];
	int qtd = 0;

	if (vertice_i == vertice_j)
		return 0;

	while (vertice_atual != vertice_j) {
		//Entao pertece ao conjunto_k

		if (qtd == 0) {
			vertices_no_caminho = (int *) malloc(sizeof(int));
		} else {
			vertices_no_caminho = (int *) realloc(vertices_no_caminho,
					(qtd + 1) * sizeof(int));
		}
		vertices_no_caminho[qtd] = vertice_atual;
		qtd++;

		vertice_atual = rota[vertice_atual];
	}

	return qtd;

}

int obter_vertices_k(int *rota, int vertice_i, int vertice_j, NP *NP_lista,
		int *&vertices_k, int p) {

	int sucessor_vertice_i = rota[vertice_i];
	int qtd_k = 0;
	int *vertices_caminho_ji;

	//printf("Obtendo vertices K! \n");
	//printf("valor de p = %d \n", p);

	/*printf("Vertice NP( i+1=%d)-->", sucessor_vertice_i);
	 for (int i = 0; i < p; ++i) {
	 printf("%d, ", NP_lista[sucessor_vertice_i].vizinhos[i]);
	 }
	 printf("\n");
	 */
	int qtd_vertices_no_caminho_ji = vertices_no_caminho(rota, vertice_j,
			vertice_i, vertices_caminho_ji);

	/*
	 printf("\nVertices entre j=%d e i = %d --> ", vertice_j, vertice_i);
	 for (int i = 0; i < qtd_vertices_no_caminho_ji; ++i) {
	 printf("%d, ", vertices_caminho_ji[i]);
	 }
	 printf("\n");
	 */
	if (qtd_vertices_no_caminho_ji == 0) {
		return 0;
	}

	for (int i = 0; i < p; ++i) {
		if (in_int_vetor(NP_lista[sucessor_vertice_i].vizinhos[i],
				vertices_caminho_ji, qtd_vertices_no_caminho_ji) == 1) {
			qtd_k++;
		}
	}

	if (qtd_k > 0) {
		vertices_k = (int *) malloc(qtd_k * sizeof(int));
		int cont = 0;
		for (int i = 0; i < p; ++i) {
			if (in_int_vetor(NP_lista[sucessor_vertice_i].vizinhos[i],
					vertices_caminho_ji, qtd_vertices_no_caminho_ji) == 1) {
				/*
				 if (qtd_k == 0){
				 vertices_k = (int *) malloc(sizeof(int));
				 }else {
				 vertices_k = (int *) realloc(vertices_k, (qtd_k+1)*sizeof(int));
				 }
				 vertices_k[qtd_k] = NP_lista[sucessor_vertice_i].vizinhos[i];
				 */

				vertices_k[cont] = NP_lista[sucessor_vertice_i].vizinhos[i];
				cont++;

			}
		}

	}

	free(vertices_caminho_ji);

	return qtd_k;
}

int obter_vertices_l(int *rota, int vertice_i, int vertice_j, NP *NP_lista,
		int *&vertices_l, int p) {

	int sucessor_vertice_j = rota[vertice_j];
	//int sucessor_vertice_i = rota[vertice_i];
	//int vertice_atual = rota[vertice_i];
	int qtd_l = 0;
	//printf("Obtendo vertices L! \n");
	////imprimir_solucao(rota, vertices);
	//printf("valor de p = %d \n", p);
	/*printf("Vertice NP( j+1=%d)-->", sucessor_vertice_j);
	 for (int i = 0; i < p; ++i) {
	 printf("%d, ", NP_lista[sucessor_vertice_j].vizinhos[i]);
	 }
	 printf("\n");
	 */

	//imprimir_solucao(rota, vertices);
	//printf("\nVertices entre i=%d e j = %d --> ", vertice_i,  vertice_j);
	int *vertices_caminho_ij;
	int qtd_vertices_no_caminho_ij = vertices_no_caminho(rota, vertice_i,
			vertice_j, vertices_caminho_ij);

	if (qtd_vertices_no_caminho_ij == 0) {
		return 0;
	}

	for (int i = 0; i < p; ++i) {
		if (in_int_vetor(NP_lista[sucessor_vertice_j].vizinhos[i],
				vertices_caminho_ij, qtd_vertices_no_caminho_ij) == 1) {

			if (qtd_l == 0) {
				vertices_l = (int *) malloc(sizeof(int));
			} else {
				vertices_l = (int *) realloc(vertices_l,
						(qtd_l + 1) * sizeof(int));
			}

			vertices_l[qtd_l] = NP_lista[sucessor_vertice_j].vizinhos[i];

			qtd_l++;

		}
	}

	/* Old Version
	 while (vertice_atual != vertice_j){

	 if ((in_int_vetor(vertice_atual, NP_lista[sucessor_vertice_j].vizinhos, p) == 1) and (vertice_atual != sucessor_vertice_i)){
	 //Entao pertece ao conjunto_l
	 //printf("%d (L=%d), ",vertice_atual, qtd_l);
	 vertices_l[qtd_l] = vertice_atual;
	 qtd_l++;
	 vertices_l = (int *) realloc(vertices_l, (qtd_l+1)*sizeof(int));
	 }else {
	 //printf("%d, ",vertice_atual);
	 }
	 vertice_atual = rota[vertice_atual];

	 }*/
	free(vertices_caminho_ij);
	return qtd_l;
}

int obter_vertices_k_us1(int *rota, int vertice_i, int vertice_j, NP *NP_lista,
		int *&vertices_k, int p) {

	//printf("Obtendo Vertices K! \n");

	int sucessor_vertice_i = rota[vertice_i];
	int antecessor_vertice_i = antecessor_de(rota, vertice_i);
	int antecessor_vertice_j = antecessor_de(rota, vertice_j);

	int qtd_k = 0;

	if (sucessor_vertice_i == antecessor_vertice_j) {
		return 0;
	}

	////imprimir_solucao(rota, vertices);
	//printf("Vertice NP( i-1=%d)-->", antecessor_vertice_i);
	/*for (int i = 0; i < p; ++i) {
	 printf("%d, ", NP_lista[antecessor_vertice_i].vizinhos[i]);
	 }
	 */
	//printf("Vertices i (%d) e vertices j (%d)\n", vertice_i, vertice_j);
	//printf("\nVertices entre i+1(%d) e j-1(%d): ", sucessor_vertice_i, antecessor_vertice_j);
	//New Version
	int *vertices_caminho;
	int qtd_vertices_no_caminho = vertices_no_caminho(rota, sucessor_vertice_i,
			antecessor_vertice_j, vertices_caminho);

	if (qtd_vertices_no_caminho == 0) {
		return 0;
	}

	for (int i = 0; i < p; ++i) {
		if (in_int_vetor(NP_lista[antecessor_vertice_i].vizinhos[i],
				vertices_caminho, qtd_vertices_no_caminho) == 1) {

			if (qtd_k == 0) {
				vertices_k = (int *) malloc(sizeof(int));
			} else {
				vertices_k = (int *) realloc(vertices_k,
						(qtd_k + 1) * sizeof(int));
			}
			vertices_k[qtd_k] = NP_lista[antecessor_vertice_i].vizinhos[i];

			qtd_k++;

		}
	}

	free(vertices_caminho);

	/* Old Version

	 while (vertice_atual != antecessor_vertice_j){
	 ////printf("%d, ",vertice_atual);
	 if (in_int_vetor(vertice_atual, NP_lista[antecessor_vertice_i].vizinhos, p) == 1){
	 //Entao pertece ao conjunto_k
	 vertices_k[qtd_k] = vertice_atual;
	 qtd_k++;
	 vertices_k = (int *) realloc(vertices_k, (qtd_k+1)*sizeof(int));
	 }
	 vertice_atual = rota[vertice_atual];
	 }
	 */

	return qtd_k;
}

int obter_vertices_k_us2(int *rota, int vertice_i, int vertice_j, NP *NP_lista,
		int *&vertices_k, int p) {

	//printf("Obtendo Vertices K - US2! \n");

	int sucessor_vertice_i = rota[vertice_i];
	int sucessor_vertice_j = rota[vertice_j];
	int antecessor_vertice_i = antecessor_de(rota, vertice_i);
	int antecessor_vertice_j = antecessor_de(rota, vertice_j);
	int *vertices_caminho_ji;
	int qtd_k = 0;

	if (sucessor_vertice_i == antecessor_vertice_j) {
		return 0;
	}

	////imprimir_solucao(rota, vertices);
	//printf("Vertice NP( i-1=%d)-->", antecessor_vertice_i);
	/*for (int i = 0; i < p; ++i) {
	 printf("%d, ", NP_lista[antecessor_vertice_i].vizinhos[i]);
	 }
	 */
	//printf("Vertices i (%d) e vertices j (%d)\n", vertice_i, vertice_j);
	//printf("\nVertices entre j+1(%d) e i-2(%d): ", sucessor_vertice_j, antecessor_de(rota, antecessor_vertice_i));
	int qtd_no_caminho = vertices_no_caminho(rota, sucessor_vertice_j,
			antecessor_de(rota, antecessor_vertice_i), vertices_caminho_ji);

	if (qtd_no_caminho == 0) {
		return 0;
	}

	/*
	 for (int i = 0; i < qtd_no_caminho; ++i) {
	 printf("%d, ", vertices_caminho_ji[i]);
	 }
	 printf("\n");
	 */

	for (int i = 0; i < p; ++i) {
		//CORRIGIR: e o antecessor de i !
		//if (in_int_vetor(NP_lista[sucessor_vertice_i].vizinhos[i],vertices_caminho_ji , qtd_no_caminho) == 1){
		if (in_int_vetor(NP_lista[antecessor_vertice_i].vizinhos[i],
				vertices_caminho_ji, qtd_no_caminho) == 1) {

			if (qtd_k == 0) {
				vertices_k = (int *) malloc(sizeof(int));
			} else {
				vertices_k = (int *) realloc(vertices_k,
						(qtd_k + 1) * sizeof(int));
			}
			//CORRIGIR: e o antecessor de i !
			//vertices_k[qtd_k] = NP_lista[sucessor_vertice_i].vizinhos[i];
			vertices_k[qtd_k] = NP_lista[antecessor_vertice_i].vizinhos[i];

			qtd_k++;

		}
	}

	free(vertices_caminho_ji);

	return qtd_k;
}

int obter_vertices_l_us2(int *rota, int vertice_k, int vertice_j, NP *NP_lista,
		int *&vertices_l, int p) {

	//printf("Obtendo Vertices L - US2! \n");

	int sucessor_vertice_k = rota[vertice_k];
	//int antecessor_vertice_k = antecessor_de(rota, vertice_k);
	int *vertices_caminho_jk;
	int qtd_l = 0;

	////imprimir_solucao(rota, vertices);
	//printf("Vertice NP( k+1=%d)-->", sucessor_vertice_k);
	/*for (int i = 0; i < p; ++i) {
	 printf("%d, ", NP_lista[sucessor_vertice_k].vizinhos[i]);
	 }
	 */
	//printf("Vertices k (%d) e vertices j (%d)\n", vertice_k, vertice_j);
	//printf("\nVertices entre j(%d) e k-1(%d): ", vertice_j, antecessor_vertice_k);
	int qtd_no_caminho = vertices_no_caminho(rota, vertice_j,
			antecessor_de(rota, vertice_k), vertices_caminho_jk);

	if (qtd_no_caminho == 0) {
		return 0;
	}
	/*
	 for (int i = 0; i < qtd_no_caminho; ++i) {
	 printf("%d, ", vertices_caminho_jk[i]);
	 }
	 printf("\n");
	 */

	for (int i = 0; i < p; ++i) {
		if (in_int_vetor(NP_lista[sucessor_vertice_k].vizinhos[i],
				vertices_caminho_jk, qtd_no_caminho) == 1) {

			if (qtd_l == 0) {
				vertices_l = (int *) malloc(sizeof(int));
			} else {
				vertices_l = (int *) realloc(vertices_l,
						(qtd_l + 1) * sizeof(int));
			}
			vertices_l[qtd_l] = NP_lista[sucessor_vertice_k].vizinhos[i];

			qtd_l++;

		}
	}

	free(vertices_caminho_jk);

	return qtd_l;
}

int inserir_tipo_1(Vertice *vertices, int *rota, int vertice_v, int vertice_i,
		int vertice_j, NP *NP_lista, int qtd_n, int p, Matriz *C,
		int *&rota_obtida) {
	//printf("\nInserir TIPO 1 (p = %d): \n", p);
	int *nova_rota = (int *) malloc(qtd_n * sizeof(int));
	int *melhor_rota = (int *) malloc(qtd_n * sizeof(int));
	int melhor_custo = 1000000000;

	//Fazer um copia da solucao_inicial
	copy_int_vetor(rota, nova_rota, qtd_n);

	//Se vertice_i e vertice_j sao adjacentes, inserir v entre eles
	if (rota[vertice_i] == vertice_j) {
		nova_rota[vertice_i] = vertice_v;
		nova_rota[vertice_v] = vertice_j;
		copy_int_vetor(nova_rota, rota_obtida, qtd_n);
		//printf("Vertices adjacentes! a \n");
		free(nova_rota);
		free(melhor_rota);
		return 0;

	} else if (rota[vertice_j] == vertice_i) {
		nova_rota[vertice_j] = vertice_v;
		nova_rota[vertice_v] = vertice_i;
		copy_int_vetor(nova_rota, rota_obtida, qtd_n);
		//printf("Vertices adjacentes! b \n");
		free(nova_rota);
		free(melhor_rota);
		return 0;
	}

	//Vertices K, pertence a NP(sucessor_i) e a rota_parcial entre vertice_j e vertice_i
	int * vertices_k;	// = (int *) malloc(sizeof(int));

	//Obter o vertices k
	////printf("Obter vertices K!\n");
	int qtd_k = obter_vertices_k(rota, vertice_i, vertice_j, NP_lista,
			vertices_k, p);
	//printf("Vertices K Obtidos (qtd=%d)!\n", qtd_k);
	/*for (int i = 0; i < qtd_k; i++) {
	 printf("%d ", vertices_k[i]);
	 }*/

	if (qtd_k == 0) { //ou seja, nao ha vertices entre j e i que tb pertecam a NP(i+1)
		//printf("Qtd_k = 0 !\n;");
		free(nova_rota);
		free(melhor_rota);
		//free(vertices_k);
		//rota_obtida = NULL;
		return 1;
	}

	int sucessor_vertice_k = -1;
	int vertice_k = -1;
	int sucessor_vertice_i = rota[vertice_i];
	int sucessor_vertice_j = rota[vertice_j];

	//Para cada vertice_k
	for (int k = 0; k < qtd_k; ++k) {

		vertice_k = vertices_k[k];
		if (vertice_i == vertice_k or vertice_j == vertice_k)
			continue;

		//printf("\nUsando vertice_k = %d \n", vertice_k);
		sucessor_vertice_k = rota[vertice_k];

		//printf("V = %d, I = %d, J = %d, K = %d, J+1 = %d, K+1 = %d, I+1 = %d\n", vertice_v, vertice_i, vertice_j, vertice_k, sucessor_vertice_j, sucessor_vertice_k, sucessor_vertice_i);

		//imprimir_solucao(nova_rota, vertices);

		inverter_sentido_rota(nova_rota, sucessor_vertice_i, vertice_j);
		inverter_sentido_rota(nova_rota, sucessor_vertice_j, vertice_k);

		nova_rota[vertice_i] = vertice_v;
		nova_rota[vertice_v] = vertice_j;

		nova_rota[sucessor_vertice_i] = vertice_k;
		nova_rota[sucessor_vertice_j] = sucessor_vertice_k;
		//imprimir_solucao(nova_rota, vertices);
		//scanf("%d", &i);

		if (custo(C, nova_rota) < melhor_custo) {
			melhor_custo = custo(C, nova_rota);
			copy_int_vetor(nova_rota, melhor_rota, qtd_n);
		}

		copy_int_vetor(rota, nova_rota, qtd_n);
	}

	copy_int_vetor(melhor_rota, rota_obtida, qtd_n);

	free(nova_rota);
	free(melhor_rota);
	free(vertices_k);

	if (melhor_custo == 1000000000) {
		return 1;
	} else {
		return 0;
	}

}

int inserir_tipo_2(int *rota, int vertice_v, int vertice_i, int vertice_j,
		NP *NP_lista, int qtd_n, int p, Matriz *C, int *&rota_obtida) {

	//printf("\nInserir TIPO 2: \n");
	//printf("Vertice V, I e J = %d, %d e %d\n", vertice_v, vertice_i, vertice_j);

	int *nova_rota = (int *) malloc(qtd_n * sizeof(int));
	int *melhor_rota = (int *) malloc(qtd_n * sizeof(int));
	int melhor_custo = 1000000000;

	//Vertices K, pertence a NP(sucessor_i) e a rota_parcial entre vertice_j e vertice_i

	//Fazer um copia da solucao_inicial
	copy_int_vetor(rota, nova_rota, qtd_n);

	//Se vertice_i e vertice_j sao adjacentes, inserir v entre eles
	if (rota[vertice_i] == vertice_j) {
		nova_rota[vertice_i] = vertice_v;
		nova_rota[vertice_v] = vertice_j;
		copy_int_vetor(nova_rota, rota_obtida, qtd_n);
		free(nova_rota);
		free(melhor_rota);
		return 0;

	} else if (rota[vertice_j] == vertice_i) {
		nova_rota[vertice_j] = vertice_v;
		nova_rota[vertice_v] = vertice_i;
		copy_int_vetor(nova_rota, rota_obtida, qtd_n);
		free(nova_rota);
		free(melhor_rota);
		return 0;
	}

	int * vertices_k;	// = (int *) malloc(sizeof(int));

	//Obter vertices k: Vertices entre j e i que tb pertecam a NP(i+1)
	////printf("Obter vertices K!\n");
	int qtd_k = obter_vertices_k(rota, vertice_i, vertice_j, NP_lista,
			vertices_k, p);
	////printf("Vertices K Obtidos (qtd=%d)!\n", qtd_k);

	if (qtd_k == 0) {
		//rota_obtida = NULL;
		free(nova_rota);
		free(melhor_rota);
		//free(vertices_k);
		return 1;
	}

	int antecessor_vertice_k = -1;
	int vertice_k = -1;
	int antecessor_vertice_l = -1;
	int vertice_l = -1;
	int sucessor_vertice_i = rota[vertice_i];
	int sucessor_vertice_j = rota[vertice_j];

	int * vertices_l;	 // = (int *) malloc(sizeof(int));
	//Obter vertices l: Vertices entre i e j
	int qtd_l = obter_vertices_l(rota, vertice_i, vertice_j, NP_lista,
			vertices_l, p);
	/*printf("Vertices L Obtidos (qtd=%d) --> \n", qtd_l);
	 for (int i = 0; i < qtd_l; ++i) {
	 printf("%d ", vertices_l[i]);
	 }
	 printf("\n");
	 */
	if (qtd_l == 0) {
		//rota_obtida = NULL;
		free(nova_rota);
		free(melhor_rota);
		free(vertices_k);
		//free(vertices_l);
		return 1;
	}

	//Para cada vertice_k
	for (int k = 0; k < qtd_k; ++k) {
		vertice_k = vertices_k[k];
		antecessor_vertice_k = antecessor_de(rota, vertice_k);
		//printf("Usando vertice_k = %d (antecessor=%d)\n", vertice_k, antecessor_vertice_k);

		if (vertice_k == vertice_j or vertice_k == sucessor_vertice_j) {

			//Remover ??
			/*
			 if (k == qtd_k-1){
			 rota_obtida = NULL;
			 free(vertices_l);
			 free(vertices_k);
			 free(nova_rota);
			 free(melhor_rota);
			 return;
			 }*/
			// printf("K == J ou ao Sucesso de J !!\n");
			continue;
		}

		for (int l = 0; l < qtd_l; ++l) {

			vertice_l = vertices_l[l];
			antecessor_vertice_l = antecessor_de(rota, vertice_l);
			//printf("Usando vertice_l = %d (antecessor=%d)\n", vertice_l, antecessor_vertice_l);

			if (vertice_l == vertice_i or vertice_l == sucessor_vertice_i
					or vertice_l == sucessor_vertice_j) {
				continue;
			}

			////imprimir_solucao(nova_rota);
			//printf("vertice_v: %d\n", vertice_v);
			//printf("vertice_i: %d\n", vertice_i);
			//printf("vertice_j: %d\n", vertice_j);
			//printf("vertice_l: %d\n", vertice_l);
			//printf("vertice_k: %d\n", vertice_k);

			inverter_sentido_rota(nova_rota, vertice_l, vertice_j);
			inverter_sentido_rota(nova_rota, sucessor_vertice_i,
					antecessor_vertice_l);
			//inserir v entre i e j
			nova_rota[vertice_i] = vertice_v;
			nova_rota[vertice_v] = vertice_j;
			//rearranjar solucao
			nova_rota[vertice_l] = sucessor_vertice_j;
			nova_rota[antecessor_vertice_k] = antecessor_vertice_l;
			nova_rota[sucessor_vertice_i] = vertice_k;

			//printf("a\n");
			////imprimir_solucao(nova_rota);
			//printf("b\n");
			if (custo(C, nova_rota) < melhor_custo) {
				melhor_custo = custo(C, nova_rota);
				copy_int_vetor(nova_rota, melhor_rota, qtd_n);
			}

			copy_int_vetor(rota, nova_rota, qtd_n);

		}
	}

	////imprimir_solucao(rota_obtida);

	//printf("Fim inserir Tipo 2");

	if (melhor_custo == 1000000000) {
		//rota_obtida = NULL;
	} else {
		copy_int_vetor(melhor_rota, rota_obtida, qtd_n);
	}

	free(vertices_l);
	free(vertices_k);
	free(nova_rota);
	free(melhor_rota);

	if (melhor_custo == 1000000000) {
		return 1;
	} else {
		return 0;
	}

}

//Fase de Insercao (GENI)
void GENI(Vertice* vertices, int *&rota, Matriz *C, int vertice_v, int p,
		int qtd_n) {

	//printf("\n<--------------------- Inicio GENI(v = %d) -------------------->\n", vertice_v);
	////imprimir_solucao(s_inicial);

	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);
	/*
	 //printf("Tamanho da rota: %d Vertices da Rota: ", tamanho_da_rota);
	 for (int i = 0; i < tamanho_da_rota; ++i) {
	 //printf(" --> %d", vertices_solucao[i]);
	 }
	 //printf("\n");
	 */
	//int tamanho_lc = qtd_lc;
	// Criar lista NP, que pegar os P vertices em S mais proximo de cada v, tanto de S quanto de LC
	NP* NP_lista = (NP *) (malloc(qtd_n * sizeof(NP)));
	//printf("Criar NP - Lista\n");
	//gerar_NP_lista(rota, C, NP_lista, p, qtd_n, vertices_solucao, tamanho_da_rota);
	//printf("NP - Lista Criada\n");
	//imprimir_np_lista(qtd_n, p, vertices_solucao, tamanho_da_rota, NP_lista);
	int *rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
	int *melhor_rota = (int *) (malloc(qtd_n * sizeof(int)));
	int melhor_custo = 1000000000;
	int retorno_insercao;
	int vertice_i = -1, vertice_j = -1;

	if (p > tamanho_da_rota) {
		p = tamanho_da_rota;
	}

	//Vizinhos de V
	NP_lista[vertice_v].vertice = vertice_v;
	NP_lista[vertice_v].vizinhos = (int *) malloc(p * sizeof(int));
	pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C, NP_lista,
			vertice_v, p);

	// Inserir v entre i e j pertencente a NP(v) (adjacentes ou não)
	for (int i = 0; i < p; ++i) {
		/* Pegar os P vertices em S mais proximo de NP_S[i].vertice
		 * e armazenar em NP_S[i].vizinhos
		 */

		for (int j = 0; j < p; ++j) {
			if (i == j)
				continue;

			//Inverter sentido da rota
			for (int sentido = 0; sentido < 2; ++sentido) {
				//printf("I e J = %d e %d \n", i, j);
				vertice_i = NP_lista[vertice_v].vizinhos[i];
				vertice_j = NP_lista[vertice_v].vizinhos[j];
				//printf("Tentar insercao de V (%d), entre os vertices(%d e %d)\n", vertice_v, vertice_i, vertice_j);

				//Vizinhos de I+1
				NP_lista[rota[vertice_i]].vertice = rota[vertice_i];
				NP_lista[rota[vertice_i]].vizinhos = (int *) malloc(
						p * sizeof(int));
				pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C,
						NP_lista, rota[vertice_i], p);

				//Vizinhos de J+1
				NP_lista[rota[vertice_j]].vertice = rota[vertice_j];
				NP_lista[rota[vertice_j]].vizinhos = (int *) malloc(
						p * sizeof(int));
				pegar_vertices_proximos_v2(vertices_solucao, tamanho_da_rota, C,
						NP_lista, rota[vertice_j], p);

				copy_int_vetor(rota, rota_obtida, qtd_n);
				//printf("IT 1-Custo antes:-->(%.1f) \n", custo( C, rota_obtida));
				//Inserir TIPO I (usando todos os k pertencente a NP(i+1) e ao caminho ji )
				retorno_insercao = inserir_tipo_1(vertices, rota, vertice_v,
						vertice_i, vertice_j, NP_lista, qtd_n, p, C,
						rota_obtida);
				/*
				 if (rota_obtida == NULL){
				 free(rota_obtida);
				 rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
				 }else
				 */
				if (retorno_insercao == 0)
					if (custo(C, rota_obtida) < melhor_custo) {

						//printf("IT 1-Custo melhor!! / Melhor Rota Atual:-->(%.1f) \n", custo( C, rota_obtida));
						melhor_custo = custo(C, rota_obtida);
						copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
						////imprimir_solucao(melhor_rota);
					}

				//Inserir TIPO II (Usandos todos os k pertencentes a NP(i + 1) no caminho ji)
				//e todos os vertices L no caminho ij, pertencentes a NP(j+1)
				copy_int_vetor(rota, rota_obtida, qtd_n);
				//printf("IT 2-Custo antes:-->(%.1f) \n", custo( C, rota_obtida));
				//printf("\n\n");
				//imprimir_solucao(rota_obtida, vertices);

				retorno_insercao = inserir_tipo_2(rota, vertice_v, vertice_i,
						vertice_j, NP_lista, qtd_n, p, C, rota_obtida);
				////printf("ponto 1\n");
				if (retorno_insercao == 0)
					if (custo(C, rota_obtida) < melhor_custo) {
						////printf("ponto 3\n");
						//printf("IT2-Custo melhor!! / Melhor Rota Atual:-->(%.1f) \n", custo( C, rota_obtida));
						melhor_custo = custo(C, rota_obtida);
						copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
						////imprimir_solucao(melhor_rota);
					}
				//printf("ponto 4\n");

				free(NP_lista[rota[vertice_i]].vizinhos);
				free(NP_lista[rota[vertice_j]].vizinhos);

				//Mudando sentido da rota
				inverter_sentido_rota(rota, vertices_solucao[0],
						vertices_solucao[0]);

			}

		}

	}
	//Atualizar o Rota com Melhor Rota Obtida
	//copy_int_vetor(melhor_rota, rota, qtd_n);

	// Tentar inserir t, na melhor posicao entre os pares de vertices adjacentes
	int primeiro_vertice = vertices_solucao[0];
	int vertice_x = primeiro_vertice;
	int sucessor_vertice_x = rota[primeiro_vertice];

	while (primeiro_vertice != sucessor_vertice_x) {

		/*//De acordo com artigo deve-se tentar inserir somente entre os vertices adjacentes da e estejam em NP(v)
		 if (in_int_vetor(primeiro_vertice, NP_lista[vertice_v].vizinhos, p) == 0 ||
		 in_int_vetor(sucessor_vertice_x, NP_lista[vertice_v].vizinhos, p) == 0){
		 vertice_x = sucessor_vertice_x;
		 sucessor_vertice_x = rota[vertice_x];
		 continue;
		 }*/

		copy_int_vetor(rota, rota_obtida, qtd_n);

		inserir_k_na_solucao(rota_obtida, vertice_v, vertice_x);

		if (custo(C, rota_obtida) < melhor_custo) {

			//printf("Melhor custo entre vertices Adjacentes (%d e %d)\n", vertice_x, sucessor_vertice_x);
			melhor_custo = custo(C, rota_obtida);
			copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
			////imprimir_solucao(melhor_rota);

		}

		vertice_x = sucessor_vertice_x;
		sucessor_vertice_x = rota[vertice_x];

	}

	//Atualizar a rota, com melhor rota obtida
	copy_int_vetor(melhor_rota, rota, qtd_n);

	//printf("Melhor Custo --> %d\n", melhor_custo);
	////imprimir_solucao(rota, vertices);

	//printf("\n<------------------------------ Fim GENI(%d) ---------------------------->\n", vertice_v);

	//Liberar Memória
	free(NP_lista[vertice_v].vizinhos);
	/*for (int i = 0; i < qtd_n; ++i) {
	 free(NP_lista[i].vizinhos);
	 }*/
	free(NP_lista);
	free(rota_obtida);
	free(melhor_rota);
	free(vertices_solucao);

}

void remover_tipo_1(Vertice* vertices, int* rota, int vertice_i, int p,
		NP* NP_lista, int qtd_n, Matriz *C, int *&rota_obtida, int in_us) {

	int vertice_j = -1;
	int vertice_k = -1;
	int antecessor_vertice_i = -1;
	int sucessor_vertice_i = -1;
	int sucessor_vertice_j = -1;

	int qtd_k = -1;

	antecessor_vertice_i = antecessor_de(rota, vertice_i);
	sucessor_vertice_i = rota[vertice_i];

	int *nova_rota = (int *) malloc(qtd_n * sizeof(int));
	int *melhor_rota = (int *) malloc(qtd_n * sizeof(int));
	int melhor_custo = custo(C, rota);
	int * vertices_k;

	copy_int_vetor(rota, nova_rota, qtd_n);
	copy_int_vetor(rota, melhor_rota, qtd_n);

	//printf("Vertice i = %d \n Vizinhos: ", vertice_i);
	/*for (int i = 0; i < p; ++i) {
	 printf("%d ", NP_lista[sucessor_vertice_i].vizinhos[i]);
	 }
	 */

	for (int j = 0; j < p; ++j) {

		vertice_j = NP_lista[sucessor_vertice_i].vizinhos[j];

		//printf("VERTICE J = %d\n", vertice_j);

		if (vertice_j == vertice_i or vertice_j == antecessor_vertice_i)
			continue;
		//Geralmente o vertice i esta entre os proximos de i+1
		//antecessor_vertice_j = antecessor_vertice_j[s_inicial, vertice_j];
		sucessor_vertice_j = rota[vertice_j];
		////printf("VERTICE J = %d\n", vertice_j);

		//obter vertices k
		qtd_k = obter_vertices_k_us1(rota, vertice_i, vertice_j, NP_lista,
				vertices_k, p);
		//printf("qtd vertices_k = (%d);\n", qtd_k);

		if (qtd_k == 0) {
			//rota_obtida = NULL;
			//return;
			continue;
		}

		for (int k = 0; k < qtd_k; ++k) {
			vertice_k = vertices_k[k];
			//printf("Vertice k = %d\n", vertice_k);
			int sucessor_vertice_k = nova_rota[vertice_k];
			//printf("Antes de remover vertice i (%f):\n", custo(C, nova_rota));
			////imprimir_solucao(nova_rota);
			//printf("vertice_i: %d\n", vertice_i);
			//printf("vertice_j: %d\n", vertice_j);
			//printf("vertice_k: %d\n", vertice_k);

			//Retirar vertice_i da solucao e reorganizando solucao
			nova_rota[vertice_i] = -1;
			inverter_sentido_rota(nova_rota, sucessor_vertice_i, vertice_k);
			nova_rota[antecessor_vertice_i] = vertice_k;
			inverter_sentido_rota(nova_rota, sucessor_vertice_k, vertice_j);
			nova_rota[sucessor_vertice_i] = vertice_j;
			nova_rota[sucessor_vertice_k] = sucessor_vertice_j;

			//printf("Depois de remover vertice i (%f):\n", custo(C, nova_rota));
			////imprimir_solucao(nova_rota);

			if (in_us == 0) {
				copy_int_vetor(nova_rota, rota_obtida, qtd_n);
				free(nova_rota);
				free(melhor_rota);
				free(vertices_k);
				return;
			}

			GENI(vertices, nova_rota, C, vertice_i, p, qtd_n);
			//printf("Depois de GENI vertice i (%f):\n", custo(C, nova_rota));
			////imprimir_solucao(nova_rota);

			if (custo(C, nova_rota) < melhor_custo) {
				melhor_custo = custo(C, nova_rota);
				copy_int_vetor(nova_rota, melhor_rota, qtd_n);
			}

			copy_int_vetor(rota, nova_rota, qtd_n);

		}
		//printf("\n");
		//Liberar vertices_k
		if (qtd_k > 0)
			free(vertices_k);

		qtd_k = 0;
	}

	if (in_us == 0) {
		free(rota_obtida);
		rota_obtida = NULL;
		free(nova_rota);
		free(melhor_rota);
		return;
	}

	copy_int_vetor(melhor_rota, rota_obtida, qtd_n);

	free(nova_rota);
	free(melhor_rota);

	return;

}

void remover_tipo_2(Vertice *vertices, int* rota, int vertice_i, int p,
		NP* NP_lista, int qtd_n, Matriz *C, int *&rota_obtida, int in_us) {

	int vertice_j = -1;
	int vertice_k = -1;
	int vertice_l = -1;
	int antecessor_vertice_i = -1;
	int antecessor_antecessor_vertice_i = -1;
	int antecessor_vertice_j = -1;
	int sucessor_vertice_i = -1;
	int sucessor_vertice_l = -1;
	int * vertices_k;			// = (int *) malloc(sizeof(int));
	int * vertices_l;			// = (int *) malloc(sizeof(int));
	int qtd_k = -1;
	int qtd_l = -1;

	antecessor_vertice_i = antecessor_de(rota, vertice_i);
	antecessor_antecessor_vertice_i = antecessor_de(rota, antecessor_vertice_i);
	sucessor_vertice_i = rota[vertice_i];

	int *nova_rota = (int *) malloc(qtd_n * sizeof(int));
	int *melhor_rota = (int *) malloc(qtd_n * sizeof(int));
	int melhor_custo = custo(C, rota);

	copy_int_vetor(rota, nova_rota, qtd_n);
	copy_int_vetor(rota, melhor_rota, qtd_n);

	//printf("Vertices J (NP(i+1): ");
	/*for (int j = 0; j < p; ++j){
	 printf("%d, ", NP_lista[sucessor_vertice_i].vizinhos[j]);
	 }
	 printf("\n");
	 */

	for (int j = 0; j < p; ++j) {

		vertice_j = NP_lista[sucessor_vertice_i].vizinhos[j];

		if (vertice_j == vertice_i or vertice_j == antecessor_vertice_i
				or vertice_j == sucessor_vertice_i
				or vertice_j == antecessor_antecessor_vertice_i)
			continue;

		antecessor_vertice_j = antecessor_de(rota, vertice_j);
		//printf("VERTICE J = %d\n", vertice_j);

		//obter vertices k
		qtd_k = obter_vertices_k_us2(rota, vertice_i, vertice_j, NP_lista,
				vertices_k, p);
		//printf("qtd vertices_k = (%d);\n", qtd_k);

		if (qtd_k == 0) {
			//rota_obtida = NULL;
			continue;
		}

		for (int k = 0; k < qtd_k; ++k) {
			vertice_k = vertices_k[k];
			//printf("Vertice k = %d\n", vertice_k);
			int sucessor_vertice_k = nova_rota[vertice_k];

			if (vertice_j == vertice_k or vertice_k == antecessor_vertice_j
					or sucessor_vertice_k == antecessor_vertice_j)
				continue;

			qtd_l = obter_vertices_l_us2(rota, vertice_k, vertice_j, NP_lista,
					vertices_l, p);
			//printf("Vertices L Obtidos (qtd=%d) --> \n", qtd_l);
			/*
			 for (int i = 0; i < qtd_l; ++i) {
			 printf("%d ", vertices_l[i]);
			 }
			 printf("\n");
			 */

			if (qtd_l == 0) {
				//rota_obtida = NULL;
				continue;
			}

			for (int l = 0; l < qtd_l; ++l) {
				vertice_l = vertices_l[l];
				sucessor_vertice_l = rota[vertice_l];

				if (antecessor_vertice_i == vertice_l
						or sucessor_vertice_i == vertice_l
						or sucessor_vertice_i == sucessor_vertice_l
						or antecessor_vertice_i == sucessor_vertice_l)
					continue;

				//printf("Antes de remover vertice i (%f):\n", custo(C, nova_rota));
				////imprimir_solucao(nova_rota);
				//printf("Ant-Antecessor de vertice_i: %d\n", antecessor_antecessor_vertice_i);
				//printf("Antecessor de vertice_i: %d\n", antecessor_vertice_i);
				//printf("vertice_i: %d\n", vertice_i);
				//printf("vertice_j: %d\n", vertice_j);
				//printf("vertice_k: %d\n", vertice_k);
				//printf("vertice_l: %d\n", vertice_l);

				//APLICAR: Remover tipo II - Retirar vertice_i da solucao e reorganizando solucao = OK
				nova_rota[vertice_i] = -1;
				inverter_sentido_rota(nova_rota, sucessor_vertice_i,
						antecessor_vertice_j);
				nova_rota[vertice_l] = sucessor_vertice_k;
				inverter_sentido_rota(nova_rota, sucessor_vertice_l, vertice_k);
				nova_rota[antecessor_vertice_i] = vertice_k;
				nova_rota[sucessor_vertice_l] = antecessor_vertice_j;
				nova_rota[sucessor_vertice_i] = vertice_j;

				//printf("Depois de remover vertice i (%f):\n", custo(C, nova_rota));
				////imprimir_solucao(nova_rota);

				if (in_us == 0) {
					copy_int_vetor(nova_rota, rota_obtida, qtd_n);
					free(nova_rota);
					free(melhor_rota);
					free(vertices_k);
					free(vertices_l);
					return;
				}

				GENI(vertices, nova_rota, C, vertice_i, p, qtd_n);

				//printf("Depois de GENI vertice i (%f):\n", custo(C, nova_rota));
				////imprimir_solucao(nova_rota);

				if (custo(C, nova_rota) < melhor_custo) {
					melhor_custo = custo(C, nova_rota);
					copy_int_vetor(nova_rota, melhor_rota, qtd_n);
				}

				copy_int_vetor(rota, nova_rota, qtd_n);

			}

			//reset contador de L
			if (qtd_l > 0)
				free(vertices_l);

			qtd_l = 0;

		}
		//printf("\n");

		//reset contador de k
		if (qtd_k > 0)
			free(vertices_k);

		qtd_k = 0;
	}

	if (in_us == 0) {
		free(rota_obtida);
		free(nova_rota);
		free(melhor_rota);
		rota_obtida = NULL;
		return;
	}

	copy_int_vetor(melhor_rota, rota_obtida, qtd_n);

	free(nova_rota);
	free(melhor_rota);

	return;

}

//Fase de Realocacao (Pós-Otimizacao)
void US(Vertice *vertices, int *&rota, Matriz *C, int p, int qtd_n) {

	//printf("\n<--------------------- Inicio US -------------------->\n");

	//Para todos os vertices da solucao, procura um melhor local.
	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

	NP* NP_lista = (NP *) (malloc(qtd_n * sizeof(NP)));
	gerar_NP_lista(rota, C, NP_lista, p, qtd_n, vertices_solucao,
			tamanho_da_rota);
	int *rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
	int *melhor_rota = (int *) (malloc(qtd_n * sizeof(int)));
	int melhor_custo = custo(C, rota);
	int vertice_i = -1;

	for (int i = 0; i < tamanho_da_rota; ++i) {

		vertice_i = vertices_solucao[i];

		//printf("\nInício Vertice i (%d/%d) = %d\n", i+1, tamanho_da_rota, vertice_i);

		for (int sentido = 0; sentido < 2; ++sentido) {

			remover_tipo_1(vertices, rota, vertice_i, p, NP_lista, qtd_n, C,
					rota_obtida, 1);

			if (rota_obtida == NULL) {
				free(rota_obtida);
				rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
			} else if (custo(C, rota_obtida) < melhor_custo) {
				//printf("REMOV 1-Custo melhor!! / Melhor Rota Atual:--> \n");
				melhor_custo = custo(C, rota_obtida);
				copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
				////imprimir_solucao(melhor_rota);
			}

			remover_tipo_2(vertices, rota, vertice_i, p, NP_lista, qtd_n, C,
					rota_obtida, 1);

			if (rota_obtida == NULL) {
				free(rota_obtida);
				rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
			} else if (custo(C, rota_obtida) < melhor_custo) {
				//printf("REMOV 2-Custo melhor!! / Melhor Rota Atual:--> \n");
				melhor_custo = custo(C, rota_obtida);
				copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
				////imprimir_solucao(melhor_rota);

			}

			//Mudando sentido da rota
			inverter_sentido_rota(rota, vertices_solucao[0],
					vertices_solucao[0]);

		}

		//printf("\nFim Vertice i (%d/%d) = %d\n", i+1, tamanho_da_rota, vertice_i);

	}

	if (melhor_custo < custo(C, rota)) {
		copy_int_vetor(melhor_rota, rota, qtd_n);
	}

	//printf("\n<--------------------- Fim US -------------------->\n");
	for (int i = 0; i < qtd_n; ++i) {
		free(NP_lista[i].vizinhos);
	}
	free(NP_lista);
	free(rota_obtida);
	free(melhor_rota);
	free(vertices_solucao);

}

void fixar_componente_conexa(int qtd_de_padroes, DMPadrao& padrao,
		DMPadrao* padroes, Matriz* C, int*& rota, Vertice* vertices) {

	//Fixar três CC na rota e fechar o ciclo
	int usar_3_maiores_cc = 1; // Se 0, fixará somente a maior cc na rota

	printf("\nGENIUS + MINING !!\n");
	int qtd_padroes = qtd_de_padroes;
	padrao = padroes[genrand_int32() % qtd_padroes];
	/*
	 while (padrao.qtd_cc < 10)
	 padrao = padroes[genrand_int32() % qtd_padroes];
	 */

	printf("PADRÕES MINERADOS: \n");
	printf("Qtd de Padrões: %d: \n\n", qtd_padroes);
	printf("Padrão Sorteado = %d - qtd_arestas: %d - suporte: %d \n", padrao.id,
			padrao.tamanho, padrao.suporte);
	printf("\tArestas: ");
	for (int a = 0; a < padrao.tamanho; ++a) {
		printf("[%d]%d->%d ", padrao.arestas[a].id,
				padrao.arestas[a].vertice_origem,
				padrao.arestas[a].vertice_destino);
	}
	printf("\n\n");
	//Imprimir
	for (int c = 0; c < padrao.qtd_cc; ++c) {
		printf("CC %d, qtd_arestas=%d --| %d", c, padrao.ccs[c].tamanho,
				padrao.ccs[c].aresta_id[0]);
		for (int a = 1; a < padrao.ccs[c].tamanho; ++a) {
			printf(" > %d", padrao.ccs[c].aresta_id[a]);
		}
		printf("\n");
	}

	if (padrao.qtd_cc >= 3 && usar_3_maiores_cc == 1) {
		//imprimir_solucao(rota, vertices);
		printf("FiXANDO as tres maiores CC a rota \n");
		int custos_combinacoes[8]; //Guarda o custo das diversas combinações de ligar 3 CC, a indice é a combinacao e valor custo
		ComponenteConexa componentes[3];
		componentes[0] = padrao.ccs[padrao.qtd_cc - 1];
		//primeiro e ultimo vertices da primeira CC
		int v1 = padrao.arestas[componentes[0].aresta_id[0]].vertice_origem;
		int v2 = padrao.arestas[componentes[0].aresta_id[componentes[0].tamanho
				- 1]].vertice_destino;
		componentes[1] = padrao.ccs[padrao.qtd_cc - 2];
		//primeiro e ultimo vertices da segunda CC
		int v3 = padrao.arestas[componentes[1].aresta_id[0]].vertice_origem;
		int v4 = padrao.arestas[componentes[1].aresta_id[componentes[1].tamanho
				- 1]].vertice_destino;
		componentes[2] = padrao.ccs[padrao.qtd_cc - 3];
		//primeiro e ultimo vertices da terceira CC
		int v5 = padrao.arestas[componentes[2].aresta_id[0]].vertice_origem;
		int v6 = padrao.arestas[componentes[2].aresta_id[componentes[2].tamanho
				- 1]].vertice_destino;
		printf("1 2 3 4 5 6 > %d %d %d %d %d %d \n", v1, v2, v3, v4, v5, v6);
		custos_combinacoes[0] = get(C, v2, v3) + get(C, v4, v5)
				+ get(C, v6, v1);
		custos_combinacoes[1] = get(C, v2, v3) + get(C, v4, v6)
				+ get(C, v5, v1);
		custos_combinacoes[2] = get(C, v2, v4) + get(C, v3, v5)
				+ get(C, v6, v1);
		custos_combinacoes[3] = get(C, v2, v4) + get(C, v3, v6)
				+ get(C, v5, v1);
		custos_combinacoes[4] = get(C, v2, v5) + get(C, v6, v3)
				+ get(C, v4, v1);
		custos_combinacoes[5] = get(C, v2, v5) + get(C, v6, v4)
				+ get(C, v3, v1);
		custos_combinacoes[6] = get(C, v2, v6) + get(C, v5, v3)
				+ get(C, v4, v1);
		custos_combinacoes[7] = get(C, v2, v6) + get(C, v5, v4)
				+ get(C, v3, v1);
		int menor_custo = custos_combinacoes[0];
		int indice_menor_custo = 0;
		for (int i = 1; i < 8; ++i) {
			if (custos_combinacoes[i] < menor_custo) {
				menor_custo = custos_combinacoes[i];
				indice_menor_custo = i;
			}
		}
		//Adicionar as tres CC a rota, mesmo desconexas
		for (int j = 0; j < 3; ++j) {
			printf("Adicionando CC..\n");
			for (int i = 0; i < componentes[j].tamanho; ++i) {
				Aresta aresta = padrao.arestas[componentes[j].aresta_id[i]];
				rota[aresta.vertice_origem] = aresta.vertice_destino;
				printf(" [%d](%d -> %d) ", aresta.id, aresta.vertice_origem,
						aresta.vertice_destino);
			}
			printf("\n");
		}

		//Ligar as 3 CC da rota da forma de menor custo
		switch (indice_menor_custo) {
		case 0:
			//2 - 3, 4 - 5, 6 - 1
			rota[v2] = v3;
			rota[v4] = v5;
			rota[v6] = v1;
			break;
		case 1:
			//2 - 3, 4 - 6, 5 - 1 invert(56)
			inverter_sentido_rota(rota, v5, v6);
			rota[v2] = v3;
			rota[v4] = v6;
			rota[v5] = v1;
			break;
		case 2:
			//2 - 4, 3 - 5, 6 - 1 invert(34)
			inverter_sentido_rota(rota, v3, v4);
			rota[v2] = v4;
			rota[v3] = v5;
			rota[v6] = v1;
			break;
		case 3:
			//2 - 4, 3 - 6, 5 - 1 invert(56), invert(34)
			inverter_sentido_rota(rota, v5, v6);
			inverter_sentido_rota(rota, v3, v4);
			rota[v2] = v4;
			rota[v3] = v6;
			rota[v5] = v1;
			break;
		case 4:
			//2 - 5, 6 - 3, 4 - 1
			rota[v2] = v5;
			rota[v6] = v3;
			rota[v4] = v1;
			break;
		case 5:
			//2 - 5, 6 - 4, 3 - 1 invert(34)
			inverter_sentido_rota(rota, v3, v4);
			rota[v2] = v5;
			rota[v6] = v4;
			rota[v3] = v1;
			break;
		case 6:
			//2 - 6, 5 - 3, 4 - 1 invert(56)
			inverter_sentido_rota(rota, v5, v6);
			rota[v2] = v6;
			rota[v5] = v3;
			rota[v4] = v1;
			break;
		case 7:
			//2 - 6, 5 - 4, 3 - 1 invert(34), invert(56)
			inverter_sentido_rota(rota, v3, v4);
			inverter_sentido_rota(rota, v5, v6);
			rota[v2] = v6;
			rota[v5] = v4;
			rota[v3] = v1;
			break;
		}
		imprimir_solucao(rota, vertices);
		printf("CC Fixadas... melhor combinacao = %d \n", indice_menor_custo);
		//scanf("%d",&indice_menor_custo);
	} else {
		//Ficar só a maior CC corrente
		//Adicionar a maior CC Corrente a rota
		//Pegar CC corrente para adicionar na rota
		int primeiro_vertice_cc = -1;
		int ultimo_vertice_cc = -1;
		int indice_cc = padrao.cc_atual;
		printf("Indice cc atual = %d\n", indice_cc);
		if (indice_cc == 0) {
			padroes[padrao.id].cc_atual = padrao.qtd_cc - 1;
		} else {
			padroes[padrao.id].cc_atual--;
		}
		ComponenteConexa maior_cc = padrao.ccs[indice_cc];
		for (int a = 0; a < maior_cc.tamanho; a++) {
			Aresta aresta = padrao.arestas[maior_cc.aresta_id[a]];
			//pegar primeiro e ultimo vértice da cc
			if (a == 0)
				primeiro_vertice_cc = aresta.vertice_origem;

			if (a == maior_cc.tamanho - 1)
				ultimo_vertice_cc = aresta.vertice_destino;

			rota[aresta.vertice_origem] = aresta.vertice_destino;
			printf(" (%d -> %d) ", aresta.vertice_origem,
					aresta.vertice_destino);
		}
		//printf("%d \n\n", padrao.arestas[maior_cc.aresta_id[maior_cc.tamanho-1]].vertice_destino);
		//Ligar ultimo vértice de destino ao primeiro de origem
		printf("Fechando a rota %d -> %d \n\n", ultimo_vertice_cc,
				primeiro_vertice_cc);
		rota[ultimo_vertice_cc] = primeiro_vertice_cc;
		//imprimir_solucao(rota, vertices);
		//scanf("%d", &primeiro_vertice_cc);
	}
}


void GENI_US(std::vector<Vertex> toFix, Vertice *vertices, Matriz* C,
		float alfa, int *&rota, int prize, float d, int qtd_n, int qtd_t,
		int minerou, DMPadrao *padroes, int qtd_de_padroes, int*& subProbSize) {

	int vertices_t[qtd_t];
	int LC[qtd_n]; //inicialmente qtd_v - 2
	int count_t = 0;
	int count_w = 0;
	int prize_obtido = 0;
	//int count_v = 0;
	int count_vt = 0;
	int nvs;
	int a = 0;
	int b = 0;
	int c = 0;

	//p define a quantidade de proximos (Np(v))
	int p = 5;

	//identificando vertices T para selecionar 2 aleatoriamente
	for (int i = 0; i < qtd_n; ++i) {
		rota[i] = -1; //prepar vetor que representara a rota

		////printf("Vertice: %d, Tipo: %s \n", i, vertices[i].tipo);

		if (strcmp(vertices[i].tipo, "t") == 0) {
			vertices_t[count_t] = i;
			count_t++;
		}
		if (strcmp(vertices[i].tipo, "w") == 0) {
			count_w++;
		}
		if (strcmp(vertices[i].tipo, "v/t") == 0) {
			count_vt++;
		}
	}
	//count_v = count_t + count_vt;
	//printf("\n--> |V| = %d (|T| = %d; |V/T| = %d); |W| = %d\n\n", count_v, count_t, count_vt, count_w);

	DMPadrao padrao;
	//printf("GENIUS in Action!!!\n");
	if (minerou == 1) {
		fixar_componente_conexa(qtd_de_padroes, padrao, padroes, C, rota,
				vertices);
	} else {

		//1) Criar solucao parcial com a, b, c pertence a T, aleatoriamente

		while (a == b || b == c || a == c) {
			a = vertices_t[genrand_int32() % qtd_t];
			b = vertices_t[genrand_int32() % qtd_t];
			c = vertices_t[genrand_int32() % qtd_t];
		}

		//Criar rota inicial com 3 vertices
		rota[a] = b;
		rota[b] = c;
		rota[c] = a;

	}

	//imprimir_solucao(rota, vertices);

	//2) Inicializar lista de Candidados LC
	criar_LC_4_genius(vertices, rota, LC, qtd_n, 1);

	prize_obtido = prize_acumulado(vertices, rota, qtd_n);
	int prize_teste = prize_obtido;

	//Vertices Candidatos e qtd
	int * vertices_lc = (int *) malloc(sizeof(int));
	;
	int qtd_lc;
	int flag = 0; // se 1 indica que LC ja foi atualizado com V/Ts

	//Vertice a ser inserir na solucao
	int vertice_v;

	//3) Enquanto nao atender as condicoes adicionar vertices a solucao parcial
//	bool solConstruida = false;
	bool wDescoberto = true;
	bool tForaSolucao = true;
	bool prizeIncompleto = true; // [false] - para criar uma solucao inviavel

	int contw = 0, contt = 0, contp = 0;

	while (wDescoberto || tForaSolucao || prizeIncompleto) {

//			imprimir_solucao(rota, vertices);
//			getchar();

		//4) Obter vertice v, aleatoriamente em LC
		qtd_lc = vertices_de_LC(LC, vertices_lc, qtd_n);
//		printf("qtd_lc = %d \n", qtd_lc);

		if (qtd_lc == 0 and flag == 0) { // Ou seja, ja inseriu todos os obrigatorios ( T )
		//Criar LC com os demais vertices (V/T)
		//criar_LC_4_genius(vertices, rota, a, b, c, LC, qtd_n, 0);
		//qtd_lc = vertices_de_LC(LC, vertices_lc, qtd_n);
		//Marcar flag para indicar que LC ja foi reconstruida
			flag = 1;

			//US(vertices, rota, C, p, qtd_n);
			/*
			 int *melhor_vizinho_rota = (int *)malloc(qtd_n*sizeof(int));
			 busca_local(vertices, C, rota, melhor_vizinho_rota, d, prize, qtd_n);
			 copy_int_vetor(melhor_vizinho_rota, rota, qtd_n);
			 free(melhor_vizinho_rota);
			 */
		}

		if (flag == 1) {
			criar_LC_4_genius(vertices, rota, LC, qtd_n, 0);
			qtd_lc = vertices_de_LC(LC, vertices_lc, qtd_n);
		}

		//scanf("%d", &qtd_lc);

		//vertice_v = vertices_lc[genrand_int32() % qtd_lc];

		vertice_v = criar_LRC_4_genius(vertices, vertices_lc, qtd_lc, rota,
				qtd_n, C, alfa);

		//printf("Vertice Sorteado: %d\n", vertice_v);

		//5) Adicionar v a s_inicial usando GENI
		GENI(vertices, rota, C, vertice_v, p, qtd_n);

		/*
		 * Usando a Mineração:
		 * 1) Verficar se o vértice recém inserido via GENI está no padrão sorteado
		 * 2) Se estiver e o vértice adjacente no padrão ainda não estiver na solução
		 *    adicionar de forma simples na rota parcial
		 */
		//Remover v de LC
		LC[vertice_v] = 0;

		if (minerou == 2) {
			int vertice = vertice_v;
			int qtd_vertices_inseridos = 0;
			ComponenteConexa cc;
			int indice_in_cc = in_componente_conexa(vertice, padrao, cc);
			if (indice_in_cc != -1) {
				//imprimir_solucao(rota,vertices);
				//printf("Vértice %d esta na sub-CC \n", vertice_v);
				for (int c = indice_in_cc; c < cc.tamanho; ++c) {
					Aresta aresta = padrao.arestas[cc.aresta_id[c]];

					//inserir a sub-cc a partir desta aresta até a final
					if (rota[aresta.vertice_destino] != -1)
						break; //já está na rota

					inserir_k_na_solucao(rota, aresta.vertice_destino, vertice);
					qtd_vertices_inseridos++;
					vertice = aresta.vertice_destino;
					LC[vertice] = 0;

					//printf(" (%d -> %d) ", aresta.vertice_origem, aresta.vertice_destino);
				}

				if (qtd_vertices_inseridos > 0) { //Tentar reconectar a solução a lah GENIUS

					//inserir GENIUS-Like
					int vertice_i = antecessor_de(rota, vertice);

					remover_k_da_solucao(rota, vertice);

					inserir_k_na_rota_genius_like(vertices, rota, C, vertice,
							vertice_i, p, qtd_n);

					//US(vertices,rota, C, p, qtd_n);
					if (rota[vertice] == -1) { //Se nao conseguiu inserir via GENI, inserir normal e rodar um US na rota
						inserir_k_na_solucao(rota, vertice, vertice_i);
						US(vertices, rota, C, p, qtd_n);
					}

					//scanf("%d", &vertice);
				}

			}

		}

//		prize_obtido = prize_acumulado(vertices, rota, qtd_n);

		prize_teste += vertices[vertice_v].premio;


		if (!existe_w_descoberto(vertices, C, d, rota, qtd_n)) {
			wDescoberto = false;
		}

		else {
			contw++;
		}

		if (!existe_t_fora_da_solucao(vertices, rota, qtd_n)) {
			tForaSolucao = false;
		} else {
			contt++;
		}

		if (!(prize_teste < prize)) {
			prizeIncompleto = false;
		} else {
			contp++;
		}

//		cout << "prize_obtido: " << prize_obtido << endl;
//		cout << "prize_teste:" << prize_teste << endl;
//
//		getchar();

	}

	//7) Procurar uma melhor posicao para os vertices de s_inicial
//		imprimir_solucao(rota, vertices);
	//printf("Custo After GENI                 : %0.f\n", custo(C, rota));

	//for analisys
	subProbSize[0] = contw;
	subProbSize[1] = contt;
	subProbSize[2] = contp;

//	printf("its para W %d \t", subProbSize[0]);
//	printf("its para T %d \t", subProbSize[1]);
//	printf("its para prize %d \n", subProbSize[2]);
//
//	getchar();

	//

	//Fase de Melhoramento
	US(vertices, rota, C, p, qtd_n);

	//imprimir_solucao(rota, vertices);
	//printf("Custo Pós-Intermediate_US: %0.f\n", custo(C, rota));
	//scanf("%d", &flag);

	free(vertices_lc);
	//printf("FIM GENIUS!!!\n");
//
//	printf("numero de whiles:\n");

//	printf("%d \t", contw);
//	printf("%d \t", contt);
//	printf("%d \n", contp);

//	getchar();
}

int remover_genius(Vertice *vertices, int *&rota, int vertice_v, Matriz *C,
		int qtd_n, int tamanho_da_rota, int *vertices_solucao) {
	//printf("IN remover GENIUS\n");
	int *rota_obtida1 = (int *) (malloc(qtd_n * sizeof(int)));
	int *rota_obtida2 = (int *) (malloc(qtd_n * sizeof(int)));

	NP* NP_lista = (NP *) (malloc(qtd_n * sizeof(NP)));

	gerar_NP_lista(rota, C, NP_lista, 5, qtd_n, vertices_solucao,
			tamanho_da_rota);

	//printf("Remover Tipo 1\n");
	remover_tipo_1(vertices, rota, vertice_v, 5, NP_lista, qtd_n, C,
			rota_obtida1, 0);
	//printf("Remover Tipo 2\n");
	remover_tipo_2(vertices, rota, vertice_v, 5, NP_lista, qtd_n, C,
			rota_obtida2, 0);

	//printf("Dois Tipos de Remocao ja testadas!\n");
	if (rota_obtida1 == NULL and rota_obtida2 == NULL) {
		//printf("As duas retornaram NULL \n");
		//free(rota_obtida1);
		//free(rota_obtida2);
		for (int i = 0; i < qtd_n; ++i) {
			free(NP_lista[i].vizinhos);
		}
		free(NP_lista);

		return 1;
	} else if (rota_obtida1 == NULL) {

		//printf("Rev-Tipo 1 retornou null \n");

		copy_int_vetor(rota_obtida2, rota, qtd_n);

		//free(rota_obtida1);
		free(rota_obtida2);
		for (int i = 0; i < qtd_n; ++i) {
			free(NP_lista[i].vizinhos);
		}
		free(NP_lista);

		return 0;
	} else if (rota_obtida2 == NULL) {
		//printf("Rev-Tipo 2 retornou null \n");
		copy_int_vetor(rota_obtida1, rota, qtd_n);

		free(rota_obtida1);
		//free(rota_obtida2);
		for (int i = 0; i < qtd_n; ++i) {
			free(NP_lista[i].vizinhos);
		}
		free(NP_lista);

		return 0;
	} else {
		//printf("As duas conseguiram remover V \n");
		if (custo(C, rota_obtida1) <= custo(C, rota_obtida2)) {
			copy_int_vetor(rota_obtida1, rota, qtd_n);
		} else {
			copy_int_vetor(rota_obtida2, rota, qtd_n);
		}
	}

	////imprimir_solucao(rota_obtida1, vertices);
	////imprimir_solucao(rota_obtida2, vertices);
	////imprimir_solucao(rota, vertices);

	for (int i = 0; i < qtd_n; ++i) {
		free(NP_lista[i].vizinhos);
	}
	free(NP_lista);

	free(rota_obtida1);
	free(rota_obtida2);

	return 0;
}

void inserir_k_na_rota_genius_like(Vertice* vertices, int *&rota, Matriz *C,
		int vertice_v, int vertice_i, int p, int qtd_n) {
	//int vertice_i = vertice_i;
	int sucessor_vertice_i = rota[vertice_i];
	int vertice_j = -1;	//pertence a NP(vertice)

	//Gerar lista NP
	int* vertices_solucao;
	int tamanho_da_rota = vertices_da_solucao(vertices_solucao, rota);

	NP* NP_lista = (NP *) (malloc(qtd_n * sizeof(NP)));

	for (int i = 0; i < qtd_n; i++) {
		NP_lista[i].vertice = -1;
		if (i == vertice_v || i == sucessor_vertice_i) {
			NP_lista[i].vertice = i;
			NP_lista[i].vizinhos = (int *) malloc(p * sizeof(int));
			pegar_vertices_proximos(vertices_solucao, tamanho_da_rota, C,
					NP_lista, i, p);
		}
	}
	//Fim gerar lista NP

	//Inicio: inserir1 e inserir2
	int *rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
	int *melhor_rota = (int *) (malloc(qtd_n * sizeof(int)));
	int melhor_custo = 1000000000;
	int retorno_insercao;
	if (p > tamanho_da_rota) {
		p = tamanho_da_rota;
	}
	for (int j = 0; j < p; ++j) {
		vertice_j = NP_lista[vertice_v].vizinhos[j];
		if (vertice_j == vertice_i)
			continue;

		//printf("I e J = %d e %d \n", i, j);
		//vertice_j = NP_lista[vertice_v].vizinhos[j];
		int sucessor_vertice_j = rota[vertice_j];

		NP_lista[sucessor_vertice_j].vertice = sucessor_vertice_j;
		NP_lista[sucessor_vertice_j].vizinhos = (int *) malloc(p * sizeof(int));

		pegar_vertices_proximos(vertices_solucao, tamanho_da_rota, C, NP_lista,
				sucessor_vertice_j, p);
		//printf("Tentar insercao de V (%d), entre os vertices(%d e %d)\n", vertice_v, vertice_i, vertice_j);

		copy_int_vetor(rota, rota_obtida, qtd_n);
		//printf("IT 1-Custo antes:-->(%.1f) \n", custo( C, rota_obtida));
		//Inserir TIPO I (usando todos os k pertencente a NP(i+1) e ao caminho ji )
		retorno_insercao = inserir_tipo_1(vertices, rota, vertice_v, vertice_i,
				vertice_j, NP_lista, qtd_n, p, C, rota_obtida);
		/*
		 if (rota_obtida == NULL){
		 free(rota_obtida);
		 rota_obtida = (int *) (malloc(qtd_n * sizeof(int)));
		 }else
		 */
		if (retorno_insercao == 0)
			if (custo(C, rota_obtida) < melhor_custo) {

				//printf("IT 1-Custo melhor!! / Melhor Rota Atual:-->(%.1f) \n", custo( C, rota_obtida));
				melhor_custo = custo(C, rota_obtida);
				copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
				////imprimir_solucao(melhor_rota);
			}

		//Inserir TIPO II (Usandos todos os k pertencentes a NP(i + 1) no caminho ji)
		//e todos os vertices L no caminho ij, pertencentes a NP(j+1)
		copy_int_vetor(rota, rota_obtida, qtd_n);
		//printf("IT 2-Custo antes:-->(%.1f) \n", custo( C, rota_obtida));
		//printf("\n\n");
		//imprimir_solucao(rota_obtida, vertices);

		retorno_insercao = inserir_tipo_2(rota, vertice_v, vertice_i, vertice_j,
				NP_lista, qtd_n, p, C, rota_obtida);
		////printf("ponto 1\n");
		if (retorno_insercao == 0)
			if (custo(C, rota_obtida) < melhor_custo) {
				////printf("ponto 3\n");
				//printf("IT2-Custo melhor!! / Melhor Rota Atual:-->(%.1f) \n", custo( C, rota_obtida));
				melhor_custo = custo(C, rota_obtida);
				copy_int_vetor(rota_obtida, melhor_rota, qtd_n);
				////imprimir_solucao(melhor_rota);
			}
		//printf("ponto 4\n");

	}

	//Atualizar a rota, com melhor rota obtida
	copy_int_vetor(melhor_rota, rota, qtd_n);

	//printf("Melhor Custo --> %d\n", melhor_custo);
	////imprimir_solucao(rota, vertices);

	//printf("\n<------------------------------ Fim GENI(%d) ---------------------------->\n", vertice_v);

	//Liberar Memória
	for (int i = 0; i < qtd_n; ++i) {
		if (NP_lista[i].vertice != -1)
			free(NP_lista[i].vizinhos);
	}
	free(NP_lista);
	free(rota_obtida);
	free(melhor_rota);
	free(vertices_solucao);
}
