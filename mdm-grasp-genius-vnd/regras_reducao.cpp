/*
 * regras_reducao.cpp
 *
 *  Created on: 18/11/2012
 *      Author: rogerio
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

//Tipo Estruturado para vertice
typedef struct vertice{
	int id; //precisei..entao botei., guarda o numero do vertices (0, 1, 2, ... N)
	int x;
	int y;
	char tipo [3];
	int premio;
	int *cobre; //conjunto de vertices W, que este vertice V cobre
	int *coberto_por; //conjunto de vertices V, que cobre o este W
	int qtd_cobre;
	int qtd_coberto_por;
	int reduzido;
}Vertice;

int regra_01(Vertice *vertices, int qtd_n, char *nome_instancia){
	/* R6 do Matheus
	 * Remover todos os W que sao coberto por algum T
	 * Obs.: Nao reduz em nada a complexidade do problema!
	 */

	int qtd_w = 0;
	int qtd_w_reduzidos = 0;

	for (int i = 0; i < qtd_n; ++i) {

		if (strcmp(vertices[i].tipo,"w")==0){
			qtd_w++;
			//printf("W = %d (", i);
			vertices[i].reduzido = 0;
			for (int j = 0; j < vertices[i].qtd_coberto_por; ++j) {
				int vertice_cobertura = vertices[i].coberto_por[j];
				//printf(" %d,",vertice_cobertura);
				if ( strcmp(vertices[vertice_cobertura].tipo, "t") == 0 ){
					//printf(" (REDUZIR ESTE W!) ");
					qtd_w_reduzidos++;
					vertices[i].reduzido = 1;
					break;
				}
			}
			//printf(")\n");
		}

	}

	// printf("Quantidade W Reduzidos %d/%d \n", qtd_w_reduzidos, qtd_w);

	return qtd_w_reduzidos;
}

int regra_02(Vertice *vertices, int qtd_n, char *nome_instancia, int prize){

	/* R2 do Matheus
	 * Transforma em obrigatorio os vertices V/T o somatorio dos demais nao atinja o prize
	 * OPS!: Transformar um vertice opcional em obrigatorio é reducao?? SIM, reduz a complexidade do problema
	 */

	int premio = 0;
	int qtd_vt = 0;
	int qtd_vt_transformados = 0;

	for (int i = 0; i < qtd_n; ++i) {

		if (strcmp(vertices[i].tipo, "v/t") == 0){
			//printf("V/T = %d ", i);
			qtd_vt++;
			for (int j = 0; j < qtd_n; ++j) {
				if (strcmp(vertices[i].tipo, "w") != 0 and i != j )
					premio += vertices[j].premio;
			}

			if (premio < prize){
				//printf( "--> Transformar em Obrigatorio T");
				qtd_vt_transformados++;
				strcpy(vertices[i].tipo, "t");
			}

			//printf("\n");
			premio = 0;
		}


	}

	// printf("Quantidade V/T Transformados em T: %d/%d \n", qtd_vt_transformados, qtd_vt);

	return qtd_vt_transformados;
}

int aplicar_regras_reducao(Vertice *vertices, int qtd_n, char *nome_instancia, int prize, FILE *file){

	int qtd_w_reduzidos = regra_01(vertices, qtd_n, nome_instancia);
	int qtd_vt_transformados_em_t = regra_02(vertices, qtd_n, nome_instancia, prize);
	fprintf(file, "%s\t\t\t%d\t\t%d\n", nome_instancia, qtd_vt_transformados_em_t, qtd_w_reduzidos );

	return qtd_vt_transformados_em_t;
}






