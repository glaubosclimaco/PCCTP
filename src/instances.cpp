/*
 * instances.cpp
 *
 *  Created on: 23/05/2012
 *      Author: Rogerio
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
//#include <sys/resource.h> //Nao existe no windows :(

//Tipo Estruturado para Matriz e funcoes relacionadas
typedef struct matriz {
	int linha;
	int coluna;
	float** v;
} Matriz;

float get(Matriz* mat, int i, int j);
void setMatriz(Matriz* mat, int i, int j, float valor);

//Tipo Estruturado para vertice
typedef struct vertice {
	int id; //precisei..entao botei., guarda o numero do vertices (0, 1, 2, ... N)
	int x;
	int y;
	char tipo[4];
	int premio;
} Vertice;

int qtd_nos = 0;

int get_indice_do_menor_valor(int *vetor, int size, int *indices_excluidos,
		int size_excluidos);
int get_indice_do_maior_valor(int *vetor, int size, int *indices_excluidos,
		int size_excluidos);
FILE* abrir_arquivo(char* file_name);
int in_int_vetor(int a, int *v, int tam);

void definir_vertices_T(Matriz* C, float d, int *vertices_t, int qtd_t) {
	//Vetor que guarda qtd vertices a proximos de cada vertice
	int *qtd_proximos = (int *) malloc(qtd_nos * sizeof(int));

	//Calculando qtd vertices proximos cada vertice
//	FILE* out = fopen("saida_definir_vertices_T.txt", "wt");
	for (int i = 0; i < qtd_nos; i++) {
		int count = 0;
		for (int j = 0; j < qtd_nos; j++) {
			if (get(C, i, j) < d && i != j)
				count++;
		}
		qtd_proximos[i] = count;
	}

//	fprintf(out,
//			"***** | Quantidade Vertices Proximo a uma distancia = %.3f | ***** \n",
//			d);
//	for (int i = 0; i < qtd_nos; i++) {
//		fprintf(out, "Qtd de Vertices proximos a v( %d )= %d \n", i,
//				qtd_proximos[i]);
//	}

//	fclose(out);

	//Definindo T..

	int qtd_t_obtidos = 0;
	//printf("Quantidade de T: %d \n", qtd_t);
	//printf("Quantidade de T-Obtidos: %d \n", qtd_t_obtidos);
	int vertices_ja_selecionados[qtd_t];
	//inicilizar com -1
	for (int z = 0; z < qtd_t; ++z) {
		vertices_ja_selecionados[z] = -1;
	}

	int indice_com_menor_qtd_vizinhos = get_indice_do_menor_valor(qtd_proximos,
			qtd_nos, vertices_ja_selecionados, qtd_t);

	vertices_ja_selecionados[qtd_t_obtidos] = indice_com_menor_qtd_vizinhos;
	vertices_t[qtd_t_obtidos] = indice_com_menor_qtd_vizinhos;
	//printf("Adicionando vertice = %d ao conjunto T\n", indice_com_menor_qtd_vizinhos);

	//printf("ZERO ta no vetor = %d\n",in_int_vetor(0,vertices_ja_selecionados, qtd_t));
	//printf("1 ta no vetor = %d\n",in_int_vetor(1,vertices_ja_selecionados, qtd_t));

	qtd_t_obtidos++;

	while (qtd_t_obtidos < qtd_t) {

		//Recalcular distancias
		for (int i = 0; i < qtd_nos; i++) {
			int count = 0;
			for (int j = 0; j < qtd_nos; j++) {
				if ((in_int_vetor(i, vertices_t, qtd_t_obtidos) == 0)
						&& get(C, i, j) < d && i != j)
					count++;
			}
			qtd_proximos[i] = count / 2;
		}

		int indice_com_menor_qtd_vizinhos = get_indice_do_menor_valor(
				qtd_proximos, qtd_nos, vertices_ja_selecionados, qtd_t);
		//printf("Vertice com menor qtd vizinhos: %d \n", indice_com_menor_qtd_vizinhos);

		vertices_ja_selecionados[qtd_t_obtidos] = indice_com_menor_qtd_vizinhos;
		vertices_t[qtd_t_obtidos] = indice_com_menor_qtd_vizinhos;
		//printf("Adicionando vertice = %d ao conjunto T\n", indice_com_menor_qtd_vizinhos);

		qtd_t_obtidos++;

	}

	free(qtd_proximos);

}

void definir_vertices_W(Matriz* C, float d, int *vertices_w, int qtd_w,
		int *vertices_t, int qtd_t) {
	//Vetor que guarda qtd vertices a proximos de cada vertice
	int *qtd_proximos = (int *) malloc(qtd_nos * sizeof(int));

	//Calculando qtd vertices proximos cada vertice
//	FILE* out = fopen("saida_definir_vertices_W.txt", "wt");
	for (int i = 0; i < qtd_nos; i++) {
		int count = 0;
		for (int j = 0; j < qtd_nos; j++) {
			if (get(C, i, j) < d && i != j)
				count++;
		}
		qtd_proximos[i] = count / 2;
	}

//	fprintf(out,
//			"***** | Quantidade Vertices Proximo a uma distancia = %.3f | ***** \n",
//			d);
//	for (int i = 0; i < qtd_nos; i++) {
//		fprintf(out, "Qtd de Vertices proximos a v( %d )= %d \n", i,
//				qtd_proximos[i]);
//	}

//	fclose(out);

	//Definindo W..

	int qtd_w_obtidos = 0;
	//printf("Quantidade de T: %d \n", qtd_t);
	//printf("Quantidade de T-Obtidos: %d \n", qtd_t_obtidos);
	int vertices_ja_selecionados[qtd_w];
	//inicilizar com -1
	for (int z = 0; z < qtd_w; ++z) {
		vertices_ja_selecionados[z] = -1;
	}

	int indice_com_maior_qtd_vizinhos = get_indice_do_maior_valor(qtd_proximos,
			qtd_nos, vertices_ja_selecionados, qtd_w);
	//printf("Vertice com menor qtd vizinhos: %d \n", indice_com_menor_qtd_vizinhos);

	vertices_ja_selecionados[qtd_w_obtidos] = indice_com_maior_qtd_vizinhos;
	vertices_w[qtd_w_obtidos] = indice_com_maior_qtd_vizinhos;
	//printf("Adicionando vertice = %d ao conjunto W\n", indice_com_maior_qtd_vizinhos);

	//printf("ZERO ta no vetor = %d\n",in_int_vetor(0,vertices_ja_selecionados, qtd_t));
	//printf("1 ta no vetor = %d\n",in_int_vetor(1,vertices_ja_selecionados, qtd_t));

	qtd_w_obtidos++;

	while (qtd_w_obtidos < qtd_w) {

		//Recalcular distancias
		for (int i = 0; i < qtd_nos; i++) {
			int count = 0;
			for (int j = 0; j < qtd_nos; j++) {
				if ((in_int_vetor(i, vertices_t, qtd_t) == 0)
						&& (in_int_vetor(j, vertices_t, qtd_t) == 0)
						&& (in_int_vetor(j, vertices_w, qtd_w_obtidos) == 0)
						&& (in_int_vetor(i, vertices_w, qtd_w_obtidos) == 0)
						&& (get(C, i, j) < d) && (i != j))
					count++;
			}
			qtd_proximos[i] = count;
		}

		int indice_do_v_com_maior_qtd_vizinhos = get_indice_do_maior_valor(
				qtd_proximos, qtd_nos, vertices_ja_selecionados, qtd_w);
		//printf("Vertice com menor qtd vizinhos: %d \n", indice_com_menor_qtd_vizinhos);

		vertices_ja_selecionados[qtd_w_obtidos] =
				indice_do_v_com_maior_qtd_vizinhos;
		vertices_w[qtd_w_obtidos] = indice_do_v_com_maior_qtd_vizinhos;
		//printf("Adicionando vertice = %d ao conjunto W\n", indice_do_v_com_maior_qtd_vizinhos);

		qtd_w_obtidos++;
		//printf("Qtd Obtidos = %d\n", qtd_w_obtidos);

	}

	free(qtd_proximos);
}

void def_conjunto_W(Matriz* C, float d, int *vertices_w, int qtd_w) {
	//Vetor que guarda qtd vertices a proximos de cada vertice
	int prox[qtd_nos];

	for (int i = 0; i < qtd_nos; i++) {
		int count = 0;
		for (int j = 0; j < qtd_nos; j++) {
			if (get(C, i, j) < d && i != j)
				count++;
		}
		prox[i] = count / 2;
	}

	//Definindo W..
	int qtd_w_obtidos = 0;
	//printf("Quantidade de W: %d \n", qtd_w);
	//printf("Quantidade de W-Obtidos: %d \n", qtd_w_obtidos);

	while (qtd_w_obtidos < qtd_w) {
		int vertice = 0;
		//Corrigir..
		int maior = prox[47];

		//printf("Pegando= %d \n", qtd+1);
		for (int i = 0; i < qtd_nos; i++) {
			//Que ja nao esteja em T, e tenha a menor qtd de vizinhos a d distancia
			if ((in_int_vetor(i, vertices_w, qtd_w_obtidos) == 0)
					&& (prox[i] > maior)) {
				maior = prox[i];
				vertice = i;
			}
		}
		//printf("Adicionando vertice = %d ao conjunto W\n", vertice);
		vertices_w[qtd_w_obtidos] = vertice;
		qtd_w_obtidos++;
		//Recaular distancias desconsirando os T  e W ja identificados
		//Recalcular distancias
		for (int i = 0; i < qtd_nos; i++) {
			int count = 0;
			for (int j = 0; j < qtd_nos; j++) {
				if (in_int_vetor(prox[i], vertices_w, qtd_w_obtidos) == 0
						&& get(C, i, j) < d && (i != j))
					count++;
			}
			prox[i] = count / 2;
		}

	}

}

void def_premio(Vertice *vertices_tsp, Matriz* C, int qtd_vertice_v) {
	int omega_definido = 0;

	for (int i = 0; i < qtd_nos; i++) {
		unsigned int soma_distacia_ij = 0;
		unsigned int soma_quadrados_distacia_ij = 0;
		if (strcmp(vertices_tsp[i].tipo, "w") != 0) {
			for (int j = 0; j < qtd_nos; j++) {
				if ((strcmp(vertices_tsp[i].tipo, "w") != 0) && (i != j)) {
					soma_distacia_ij += get(C, i, j);
					soma_quadrados_distacia_ij += pow(get(C, i, j), 2);
				}
			}
			unsigned int fator_omega;
			if (omega_definido == 0) {
				fator_omega = 10;
				//Definir valor do Omega
				while (soma_distacia_ij > fator_omega * 10) {
					fator_omega = fator_omega * 10;
				}
			}

			unsigned int primeira_parcela = soma_quadrados_distacia_ij
					/ soma_distacia_ij;

			unsigned int segunda_parcela = soma_distacia_ij
					/ (qtd_vertice_v - 1);

			float terceira_parcela = (float) 1 / fator_omega;

			unsigned int premio = ((primeira_parcela) * (segunda_parcela)
					* (terceira_parcela));

			vertices_tsp[i].premio = premio;

			/*
			 printf("Soma_Distancia: %u \n", soma_distacia_ij);

			 printf("Fator Omega: %u \n", fator_omega);

			 printf("Primeira Parcela: %u \n", primeira_parcela);

			 printf("Primeira Parcela: %u \n", segunda_parcela);

			 printf("Terceira Parcela: %f \n", terceira_parcela);

			 printf("Premio = %d \n", premio);
			 */

			if (vertices_tsp[i].premio < 1) {
				printf("PREMIO INVÁLIDO!");
				exit(1);
			}
		} else {
			vertices_tsp[i].premio = -1;
		}
	}
}

void gerar_matriz_de_distancias(Vertice *V, Matriz* C) {

	// Calcular distancia euclidiana entre os vertices
	double distancia_heuclidiana;
	//gravar matriz de distancias em txt file
	FILE* out = fopen("matriz_distancia.prrcp", "wt");
	fprintf(out,
			"***** | Matriz de Distancia entre todos os vertices | ***** \n");
	fprintf(out, "V");
	for (int i = 0; i < qtd_nos; ++i) {
		if (i == 0)
			fprintf(out, "\t%d", i);
		else
			fprintf(out, "\t\t%d", i);
		//fprintf(out,"\t%d",i);
	}
	unsigned int menor_distancia;
	unsigned int maior_distancia;

	for (int i = 0; i < qtd_nos; i++) {
		fprintf(out, "\n%d", i);
		for (int j = 0; j < qtd_nos; j++) {

			distancia_heuclidiana = sqrt(
					pow((V[j].x - V[i].x), 2) + pow((V[j].y - V[i].y), 2));
			setMatriz(C, i, j, distancia_heuclidiana);

			if (i == 0 && j == 1) {
				menor_distancia = distancia_heuclidiana;
				maior_distancia = distancia_heuclidiana;
			} else if (distancia_heuclidiana < menor_distancia
					&& distancia_heuclidiana > 0) {
				menor_distancia = distancia_heuclidiana;
			} else if (distancia_heuclidiana > maior_distancia) {
				maior_distancia = distancia_heuclidiana;
			}

			if (distancia_heuclidiana < 0) {
				printf("Distancia negativa Localizada! \n");
				exit(1);
			}
			if (distancia_heuclidiana == 0.0)
				fprintf(out, "\t%.2f\t", distancia_heuclidiana);
			else
				fprintf(out, "\t%.2f", distancia_heuclidiana);
		}
	}
	fprintf(out, "\n--FIM--");
	fprintf(out, "\n Menor Distancia: %d", menor_distancia);
	fprintf(out, "\n Maior Distancia: %d", maior_distancia);
	fclose(out);
}

void load_vertices_from_file(char* file_name, Vertice *V) {
	//Carregar arquivo c/ problemas teste
	FILE* file = abrir_arquivo(file_name);
	//Pegar os Vertices (X, Y) e guardar no vetor de Vertices
	char linha[121];
	int v, i;
	float x, y;
	i = 0;
	while (fgets(linha, 121, file) != NULL) {
		int n = sscanf(linha, "%d %f %f", &v, &x, &y);
		if (n > 0) {
			V[i].id = i;
			V[i].x = (int) x;
			V[i].y = (int) y;
			i++;
		}
	}
	//Fechar o arquivo
	fclose(file);
}

int get_prrcp_instance(char* file_name, int qtd_t, int qtd_w, float d,
		int fator_omega, Vertice* vertices_tsp, Matriz* C, int qtd_vertices) {

	//Quantidade de vertices
	qtd_nos = qtd_vertices;

	//Carregar arquivo c/ problemas teste
	load_vertices_from_file(file_name, vertices_tsp);

	//Calcular distancia eucliana entre todos os Pontos
	gerar_matriz_de_distancias(vertices_tsp, C);
	printf("Matriz de Distancia (C) gerada com sucesso!\n");

	//Definir Vertice Obrigatorios (T)
	int* vertices_t;
	printf("Definir %d vertices tipo T.\n", qtd_t);
	vertices_t = (int*) (malloc(qtd_t * sizeof(int)));
	definir_vertices_T(C, d, vertices_t, qtd_t);
	printf("Vertices T defindos com sucesso!.\n");

	//Definir Vertices W
	int* vertices_w;
	vertices_w = (int*) (malloc(qtd_w * sizeof(int)));

	printf("Definir %d vertices tipo W.\n", qtd_w);
	//def_conjunto_W(C, d, vertices_w, qtd_w);
	definir_vertices_W(C, d, vertices_w, qtd_w, vertices_t, qtd_t);

	printf("Vertices W defindos com sucesso!.\n");

	//Definir V/T e atualizar tipo dos verticesTSP(todos os vertices inclusive W)
	int qtd_vt = 0;

	for (int i = 0; i < qtd_nos; i++) {
		if (in_int_vetor(i, vertices_t, qtd_t) == 1) {
			printf("%d = T\n", i);
			strcpy(vertices_tsp[i].tipo, "t");
		} else if (in_int_vetor(i, vertices_w, qtd_w) == 1) {
			strcpy(vertices_tsp[i].tipo, "w");
			printf("%d = W\n", i);
		} else {
			qtd_vt++;
			printf("%d = VT\n", i);
			strcpy(vertices_tsp[i].tipo, "v/t");
		}
	}
	printf("Vertices V/T defindos com sucesso(Qtd = %d)!.\n", qtd_vt);

	int qtd_vertice_v = qtd_nos - qtd_w;

	//Definir premio do Vertices V = (T e V/T)
	def_premio(vertices_tsp, C, qtd_vertice_v);
	printf("Premios do verticos V definidos com sucesso!.\n");

	//gravar instancia em txt
	FILE* out = fopen("instancia.prrcp", "wt");
	fprintf(out, "***** | Instancia Gerada a partir de: %s (TSPLIB)| ***** \n",
			file_name);
	fprintf(out, "Distancia D = %.2f\n", d);
	//Premios da Instancia
	int premio_t = 0;
	int premio_vt = 0;
	for (int i = 0; i < qtd_nos; i++) {
		if (strcmp(vertices_tsp[i].tipo, "t") == 0) {
			premio_t += vertices_tsp[i].premio;
		} else if (strcmp(vertices_tsp[i].tipo, "v/t") == 0) {
			premio_vt += vertices_tsp[i].premio;
		}
	}
	int premio_total = premio_t + premio_vt;
	fprintf(out, "Premio Total = %d (T = %d e V/T = %d)\n", premio_total,
			premio_t, premio_vt);
	printf("Premio Total = %d (T = %d e V/T = %d)\n", premio_total, premio_t,
			premio_vt);
	fprintf(out, "Todos Vertices,\nV\tx\ty\ttipo\tpremio\n");
	for (int i = 0; i < qtd_nos; i++) {
		fprintf(out, "%d\t%d\t%d\t%s\t\%d\n", i, vertices_tsp[i].x,
				vertices_tsp[i].y, vertices_tsp[i].tipo,
				vertices_tsp[i].premio);
	}
	fclose(out);

	return premio_total;

}
