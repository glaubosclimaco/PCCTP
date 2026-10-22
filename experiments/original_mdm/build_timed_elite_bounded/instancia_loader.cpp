/*
 * instancia.cpp
 *
 *  Created on: 16/07/2012
 *      Author: rogerio
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <string>
#include <iostream>

//Tipo Estruturado para vertice
typedef struct vertice{
	int id; //precisei..entao botei., guarda o numero do vertices (0, 1, 2, ... N)
	int x;
	int y;
	char tipo[4];
	int premio;
	int *cobre; //conjunto de vertices W, que este vertice V cobre
	int *coberto_por; //conjunto de vertices V, que cobre o este W
	int qtd_cobre;
	int qtd_coberto_por;
	int reduzido;
}Vertice;

//Tipo Estruturado para Matriz e funcoes relacionadas
typedef struct matriz{
	int linha;
	int coluna;
	float** v;
}Matriz;

float get(Matriz* mat, int i, int j);
void set(Matriz* mat, int i, int j, float valor);

void get_instancia_prrcp(char *file_name, Vertice *vertices_tsp, Matriz* C, int *qtd_nos,int *qtd_T, float *distancia, int *premio_solicitado){

	//Variaveis
	int itoken = 0;
	float ftoken = 0.0;
	int qtd_n = 0;
	int premio = 0;
	int count_t = 0;
	float d = 0.0;
	int valor_otimo = 0;

	std::fstream fin(file_name);

	//Qtd de Vertices
	fin >> itoken;
	qtd_n = itoken;
	*qtd_nos = qtd_n;

	//PRIZE
	fin >> itoken;
	premio = itoken;

	//Distancia
	fin >> ftoken;
	d = ftoken;

	//Valor Otimo
	fin >> itoken;
	valor_otimo = itoken;

	// printf("Qtd de Vertices: %d\n", qtd_n);
	// printf("Premio Solicitado: %d\n", premio);
	// printf("Distancia D: %.3f\n", d);
	// printf("Valor Otimo: %d\n", valor_otimo);

	if (!fin){
		printf("Erro ao abrir o arquivo!\n");
	}

	*premio_solicitado = premio;
	*distancia = d;

	int count_w = 0;
	//Tipo dos Vertices
	for (int i = 0; i < qtd_n; ++i) {
		fin >> itoken;
		vertices_tsp[i].id = i;
		switch (itoken) {
			case 0:
				//std::cout << "Vertice: "<< i << " eh Obrigatorio."<< std::endl;
				strcpy(vertices_tsp[i].tipo,"v/t");
				break;
			case 1:
				//std::cout << "Vertice: "<< i << " eh Opcional."<< std::endl;
				strcpy(vertices_tsp[i].tipo,"t");
				count_t++;
				break;
			case 2:
				//std::cout << "Vertice: "<< i << " eh de Cobertura."<< std::endl;
				strcpy(vertices_tsp[i].tipo,"w");
				count_w++;
				break;
			default:
				break;
		}

	}
	*qtd_T = count_t;
	//settando variável global com qtd de vertice V

	//Premio dos Vertices e
	for (int i = 0; i < qtd_n; ++i) {
			fin >> itoken;
			vertices_tsp[i].premio = itoken;
	}


	//gravar instancia em txt
	FILE* out = fopen("instancia.prrcp", "wt");
	fprintf(out,"***** | Instancia Gerada a partir de: %s (TSPLIB)| ***** \n", file_name);
	fprintf(out, "Distancia D = %.2f\n", d);
	//Premios da Instancia
	int premio_t = 0;
	int premio_vt = 0;
	for (int i = 0; i < qtd_n; i++) {
		if ( strcmp(vertices_tsp[i].tipo, "t")==0) {
			premio_t += vertices_tsp[i].premio;
		}else if (strcmp(vertices_tsp[i].tipo, "v/t") == 0) {
			premio_vt += vertices_tsp[i].premio;
		}
	}
	int premio_total = premio_t + premio_vt;

	fprintf(out, "Premio Total = %d (T = %d e V/T = %d)\n", premio_total, premio_t, premio_vt);

	fprintf(out,"Todos Vertices,\nV\tx\ty\ttipo\tpremio\n");
		for (int i = 0; i < qtd_n; i++) {
			fprintf(out,"%d\t%d\t%d\t%s\t\%d\n",i,0,0, vertices_tsp[i].tipo, vertices_tsp[i].premio );
		}
	fclose(out);

	//Matriz de Distancias
	for (int i = 0; i < qtd_n; ++i) {
		for (int j = 0; j < qtd_n; ++j) {
			fin >> ftoken;
			set(C, i, j, ftoken);
		}
	}

	//Grava Matriz de Distancias
	FILE* out2 = fopen("matriz_distancia_.prrcp", "wt");
	fprintf(out2,"***** | Matriz de Distancia entre todos os vertices | ***** \n");
	fprintf(out2,"V");

	for (int i = 0; i < qtd_n; ++i) {
		fprintf(out2,"\t%d",i);
	}
	unsigned int menor_distancia;
	unsigned int maior_distancia;

	for (int i = 0; i < qtd_n; i++) {
		fprintf(out2,"\n%d",i);
		for (int j = 0; j < qtd_n; j++) {

			if (i == 0 && j ==1 ){
				menor_distancia = get(C, i, j);
				maior_distancia = get(C, i, j);
			}else if (get(C, i, j) < menor_distancia && get(C, i, j) > 0) {
				menor_distancia = get(C, i, j);
			}else if (get(C, i, j) > maior_distancia){
				maior_distancia = get(C, i, j);
			}

			if (get(C, i, j) < 0){
				printf("Distancia negativa Localizada! \n");
				exit(1);
			}
			fprintf(out2,"\t%.2f",get(C, i, j));
		}
	}
	fprintf(out2,"\n--FIM--");
	fprintf(out2, "\n Menor Distancia: %d", menor_distancia);
	fprintf(out2, "\n Maior Distancia: %d", maior_distancia);
	fclose(out2);

	//definir conjunto de vertices de cobertura

		int aux_cobre[count_w]; //aux para evitar realloc
		int aux_coberto[qtd_n - count_w];
		int count_cobre = 0;
		int count_coberto_por = 0;
		for (int i = 0; i < qtd_n; ++i) {

			vertices_tsp[i].qtd_coberto_por = 0;
			vertices_tsp[i].qtd_cobre = 0;

				//Definir quais W cada V cobre
			if ( strcmp(vertices_tsp[i].tipo, "w") != 0 ){
				//montar conjunto de w q cobre
				count_cobre = 0;
				for (int j = 0; j < qtd_n; ++j) {
					if (strcmp(vertices_tsp[j].tipo, "w") == 0 and get(C, i, j) <= d){
						aux_cobre[count_cobre] = j;
						count_cobre++;
					}
				}
				if (count_cobre > 0)
					vertices_tsp[i].cobre = (int *) malloc(count_cobre * sizeof(int));

				vertices_tsp[i].qtd_cobre = count_cobre;
				for (int k = 0; k < count_cobre; ++k) {
					vertices_tsp[i].cobre[k] = aux_cobre[k];
				}

				//Verbose
				/*
				printf("\nVertice V/T: %d cobre %d vertices de W. sao eles: ", i, vertices_tsp[i].qtd_cobre);

				for (int k = 0; k < vertices_tsp[i].qtd_cobre; ++k) {
					printf("%d (%.1f), ", vertices_tsp[i].cobre[k], get(C, i, vertices_tsp[i].cobre[k]));
				}*/

			}else {
				//montar conjunto de V q cobrem o W atual
				count_coberto_por = 0;
				for (int j = 0; j < qtd_n; ++j) {
					if (strcmp(vertices_tsp[j].tipo, "w") != 0 and get(C, i, j) <= d){
						aux_coberto[count_coberto_por] = j;
						count_coberto_por++;
					}
				}
				if (count_coberto_por >0)
					vertices_tsp[i].coberto_por = (int *) malloc(count_coberto_por * sizeof(int));

				vertices_tsp[i].qtd_coberto_por = count_coberto_por;
				for (int k = 0; k < count_coberto_por; ++k) {
					vertices_tsp[i].coberto_por[k] = aux_coberto[k];
				}

				//Verbose
				/*
				printf("\nVertice W: %d e coberto %d vertices de V. sao eles: ", i, vertices_tsp[i].qtd_coberto_por);

				for (int k = 0; k < vertices_tsp[i].qtd_coberto_por; ++k) {
					printf("%d (%.1f), ", vertices_tsp[i].coberto_por[k], get(C, i, vertices_tsp[i].coberto_por[k]));
				}
				*/

			}

		}

}

int premio_total(Vertice* vertices_tsp, int qtd_n){
	int premio = 0;
	for (int i = 0; i < qtd_n; ++i) {
		if (strcmp(vertices_tsp[i].tipo,"w") != 0){
			premio += vertices_tsp[i].premio;
		}
	}

	return premio;

}
