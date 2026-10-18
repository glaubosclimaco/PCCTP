/*
 * principal.cpp
 *
 * Este e modulo principal, responsavel por chamar todos os auxiliares
 * para:
 * 1) Carregar os problemas-teste
 * 2) Registrar os Tempos
 * 3) Dispara o GRASP
 *
 *  Created on: 20/05/2012
 *      Author: Rogerio
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
//#include <sys/resource.h> //Nao existe no windows :(
#include <fstream>
#include <ctime>
#include<iostream>

using namespace std;

//Tipo Estruturado para Matriz e funcoes relacionadas
typedef struct matriz {
	int linha;
	int coluna;
	float* v;
} Matriz;

Matriz* create_matriz(int linhas, int colunas);
void free_matriz(Matriz* mat);

//Tipo Estruturado para vertice
typedef struct vertice {
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
} Vertice;

char path_instancias[50] = "instancias/";

//Headers
void get_instancia_prrcp(char *file_name, Vertice *vertices_tsp, Matriz* C,
		int *qtd_N, int *qtd_T, float *d, int *premio_solicitado);
float grasp_basico(Vertice *vertices, Matriz *C, float alfa, int grasp_max,
		float tempo_max, int qtd_t, int prize, int seed, float d, int qtd_n,
		FILE* out, char *nome_instancia, int premio_solicitado,
		int *&s_melhor_rota, int valor_otimo_instancia, int aplicar_mineracao,
		int &ultima_iteracao, int &iteracao_melhor_custo, int sup_min, int percentual_estabilidade);
//void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total);
int prize_acumulado(Vertice *vertices, int *s_inicial, int qtd_n);
int aplicar_regras_reducao(Vertice *vertices, int qtd_n, char *nome_instancia,
		int prize, FILE *file);
int qtd_vertices_em_s(int *s, int qtd_n);
float pegar_tempo_matheus(char nome_instancia[30], int seed);

//Fim-Headers

int main(int argc, char **argv) {

	int qtd_n = 0;

	int qtd_T = 0;

	float d = 0.0;

	//Nome do arquivo de instancia PRRCP
	char file_name[100];
	char nome_instancia[30];

	int premio_solicitado = 0;
	int valor_otimo_instancia = 0;

	//Paramentos pra Solucao
	float custo_obtido = 0.0;
	float alfa = 0.6;
	int grasp_max = 100; // 18;
	int criterio_tempo = 1;
	float tempo_max = 0; //usado quando o critério de parada for tempo..deixar zerado
	int seed = 1;
	int seed_inicio = 1;
	int seed_fim = 4;
	int qtd_problemas_teste = 144;
	int aplicar_mineracao = 0;
	int indice_instancia;
	int numIt = 70;
	int sup_min = 0;
	int percentual_estabilidade = 0;
	if (argc != 5) {
		printf("./executavel <data/nome_instancia> <sed> <sup_min> <intero para percentual de estabilidade. Ex: 5 para 5%%>\n");
		return -1;
	}

	if (argc == 5) {
		// printf("\nParametros ok!\n");
		//indice_instancia = atof(argv[1]);
		strcpy(file_name, argv[1]);
		seed = atoi(argv[2]);
		sup_min = atoi(argv[3]);
		percentual_estabilidade = atoi(argv[4]);
		//seed_fim = atoi(argv[3]);
		//file_name = argv[1];
		//alfa = atof(argv[2]);
		//seed = atoi(argv[3]);
	} else {
		printf("\nParametros incorretos.\n");
		exit(1);

	}

	//Carregar nome das instancias a partir de um arquivo de indexacao
	std::fstream get("nomes_instancias2.txt");
	FILE* out1 = fopen("solucoes_log.txt", "at");
	FILE* out2 = fopen("solucoes_finais_log.txt", "at");
	FILE* out3 = fopen("reducoes_log.txt", "at");
	fclose(out1);
	fclose(out2);
	fclose(out3);

	FILE* file_log = fopen("solucoes_log.txt", "at");
	FILE* file_solucoes_finais = fopen("solucoes_finais_log.txt", "at");
	FILE* file_reducoes = fopen("reducoes_log.txt", "at");

	//strcpy(file_name, path_instancias);
	//formar nome do arquivo
	//strcat(file_name, nome_instancia);

	//Pegar qtd de vertices para inicializar estruturas auxiliares
	std::fstream fin(file_name);

	fin >> qtd_n;
	fin >> premio_solicitado;
	fin >> d;
	fin >> valor_otimo_instancia;

	//Matriz das distancias(Euclidiana) entre os vertices
	Matriz *C = create_matriz(qtd_n, qtd_n);

	Vertice* vertices_instancia = (Vertice*) malloc(qtd_n * sizeof(Vertice));

	//Carregar instancia do arquivo
	get_instancia_prrcp(file_name, vertices_instancia, C, &qtd_n, &qtd_T, &d,
			&premio_solicitado);

//	qtd_T += aplicar_regras_reducao(vertices_instancia, qtd_n, nome_instancia,
//			premio_solicitado, file_reducoes);

	//Melhor Solucao
	int * s_melhor_rota = (int *) malloc(qtd_n * sizeof(int));
	int ultima_iteracao_executada = -1;
	int iteracao_melhor_custo = -1;

	//Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

	//Se o critério for tempo pegar o tempo de referencia..
	if (criterio_tempo == 1)
		tempo_max = 1200;	//pegar_tempo_matheus(nome_instancia, seed);

	// printf("\nTempo max = %.2f\n", tempo_max);

	// printf("%d\t", valor_otimo_instancia);

	clock_t start = clock();
	custo_obtido = grasp_basico(vertices_instancia, C, alfa, numIt, tempo_max,
			qtd_T, premio_solicitado, seed, d, qtd_n, file_log, nome_instancia,
			premio_solicitado, s_melhor_rota, valor_otimo_instancia,
			aplicar_mineracao, ultima_iteracao_executada,
			iteracao_melhor_custo, sup_min, percentual_estabilidade);
	clock_t end = clock();

//	cout << custo_obtido << endl;

	double tempo_total = (end - start) / (double) CLOCKS_PER_SEC;
	
	cout << endl;
	cout << "FIM"<<endl;
	cout << file_name << "\t ótimo:" << valor_otimo_instancia << "\t custo_obtido:" << custo_obtido << "\t tempo:" << tempo_total << endl;


	free(s_melhor_rota);

	//}

	//}
	//Executar o Algoritmo Heuristico 10x para cada seed
	/*
	 for (int i = 1; i <= 10; ++i) {

	 alfa = ((float)i/10);

	 for (seed = 1; seed <= 10; seed++) {

	 custo = grasp_basico(vertices_tsp, &C, alfa, grasp_max, qtd_T, premio_solicitado, seed, d, qtd_n, out, nome_instancia, premio_solicitado);

	 }

	 }
	 */
	for (int i = 0; i < qtd_n; ++i) {

		if (vertices_instancia[i].qtd_cobre > 0)
			free(vertices_instancia[i].cobre);

		if (vertices_instancia[i].qtd_coberto_por > 0)
			free(vertices_instancia[i].coberto_por);
	}

	free_matriz(C);

	free(vertices_instancia);
	fclose(file_log);
	fclose(file_reducoes);
	fclose(file_solucoes_finais);
	//}

	/*
	 int premio_total_ = premio_total(vertices_tsp, qtd_n);

	 printf("\nPRIZE Total: %d\n", premio_total_);

	 printf("\nPRIZE Solicitado: %.2d \n ", premio_solicitado);

	 printf("\nMelhor custo: %.2f\n", custo);

	 printf("\nGRASP Max: %d\n", grasp_max);

	 printf("\nAlfa: %.2f\n", alfa);

	 printf("\nDistancia D: %.0f\n", d);

	 */

	//Finalizar a contagem do tempo	 */
	//TODO
	return 0;

}
