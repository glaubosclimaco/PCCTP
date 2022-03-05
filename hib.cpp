//============================================================================
// Name        : hib.cpp
// Author      : Glaubos Clímaco
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <string.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "PR.h"
#include "model/Mycallback.h"
#include "model/Regras.h"

using namespace std;

#define infinito 99999999999.0
#define db 0
#define lb 0
#define hib 1
#define fix_one 0

Matriz* create_matriz(int linhas, int colunas);
void free_matriz(Matriz* mat);

//Headers

double HM(double, int &cont_parada, vector<edgeSol>& sol, int &itExecutadas,
          double gamma, double &pE, int **& matriz_BL, Vertice *vertices_grasp,
          Matriz *C, float alfa, int grasp_max, float tempo_max, int qtd_t,
          int prize, int seed, float d, char *nome_instancia,
          int premio_solicitado, int *&s_melhor_rota, int numIt,
          vector<aresta>& L, int delta, int fazLB, bool removeE0, int root,
          int qnd_vt, int qtd_n, vector<Vertex> vertices, double **matriz_custo,
          double prize_solver, bool exigir_r, double prize_V, double otimo,
          double limPrimal, int*);

void get_instancia_prrcp(char *file_name, Vertice *vertices_tsp, Matriz* C,
                         int &qtd_N, int *qtd_T, double &d, int &premio_solicitado, int &root,
                         int &qnd_w, int &qnd_vt, int &qnd_t, double &prize, double &D,
                         double &valor_otimo, vector<Vertex> &vertices, double **&matriz_custo,
                         double& prize_V, vector<int> &R_set);

int aplicarRegras(int*&, int ordem, bool debug, double** matriz_custo,
                  vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t,
                  int &qnd_w, double limPrimal, bool &exigir_r, vector<int> R_set, int**,
                  int);
bool ordenaCustoFreq(edgeSol a, edgeSol b);
void liberar_arestas(vector<int>& S, int *solucao, Vertice* vertices,
                     int**& matriz_BL, double** matriz_custo);
bool verifDesigTriangular(vector<Vertex> vertices, double **matriz_custo);
void imprimir_solucao(int *solucao, Vertice* vertices);
void get_solucao(int *solucao, Vertice* vertices, vector<int> *S);
float grasp_basico(vector<Vertex> vsol, int &, vector<edgeSol>&, int&, double,
                   double &pE, int**& matriz_BL, Vertice *vertices, Matriz *C, float alfa,
                   int grasp_max, float tempo_max, int qtd_t, int prize, int seed, float d,
                   int qtd_n, char *nome_instancia, int premio_solicitado,
                   int *&s_melhor_rota, int numIt, int*&);
void liberar_arestas(vector<aresta>& L, int *solucao, Vertice* vertices,
                     int**& matriz_BL, double** matriz_custo);
//void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total);
int prize_acumulado(Vertice *vertices, int *s_inicial, int qtd_n);
//int aplicar_regras_reducao(Vertice *vertices, int qtd_n, char *nome_instancia,
//		int prize, FILE *file);
int qtd_vertices_em_s(int *s, int qtd_n);
void printSolVector(vector<edgeSol> sol);
float pegar_tempo_matheus(char nome_instancia[30], int seed);
void printPython(double** m, vector<Vertex> v, int n, double limPrimal);
double solver(vector<Vertex>& vsol, vector<aresta>&, int, int, vector<edgeSol>&,
              int**, bool, int root, int qtd_t, int qnd_vt, int qtd_n,
              vector<Vertex> vertices, double **matriz_custo, double prize,
              bool exigir_r, double prize_V, double otimo, double limPrimal);
void contRemovidosV(vector<Vertex> vertices, int qtd_n, int** matriz_BL,
                    int qnd_vt, int qnd_t);
void reduzGrafo(vector<Vertex>& vertices, int& qnd_vt, int& qtd_n,
                double limPrimal, double** matriz_custo, int**);
bool myfunctionEdgeSol(edgeSol a, edgeSol b);

bool biggerPrize(Vertex a, Vertex b) {
	return a.premio > b.premio;
}

//void analisaE0(vector<Vertex> vertices, vector<aresta>& arestas, int qtd_n,
//		int** matriz_BL, double** matriz_custo);

//Fim-Headers

int main(int argc, char **argv) {

	int arestas_reais = 0, delta, fazLB = 0, aplicaE = 0, contE0 = 0,
	    aplicar_regras = 0, root = 0, qnd_w = 0, qnd_t = 0, qnd_vt = 0,
	    itExecutadas = 0, numIt = 70, **matriz_BL, *s_melhor_rota, seed,
	    grasp_max = 100, qtd_n, qtd_T = 0, valor_otimo_instancia = 0,
	    premio_solicitado = 0, ordem = -1;
	double gamma = 0.0, prize_V = 0.0, **matriz_custo, valor_otimo, prize, D,
	       tempo_grasp = 0.0, pE = 0.0, tempo_max = 0, d, limPrimal =
	                                   std::numeric_limits<double>::max(), custo_obtido = 0.0,
	                                                                       alfa = 0.6, tempo_regras_esparso = 0.0, tempo_regras_original = 0.0;
//	Nome do arquivo de instancia PRRCP
	char *file_name, nome_instancia[30];
	file_name = new char[100];
	vector<Vertex> vertices;
	vector<int> R_set;
	int *z, *nReducoes, *subProbSize;
	subProbSize = new int[4];
	z = new int[4];
	nReducoes = new int[15];
	bool liberar = true, exigir_r = false;
	vector<aresta> arestas, L;
	vector<edgeSol> solucao;
	vector<Vertex> vsol;
	double perfix;
	int rReduzidos, tReduzidos, wReduzidos, totalReducoes;


#if lb

	if (argc != 10) {
		printf(
		    "./executavel <data/nome_instancia> <sed> <aplica_cortes> <remove Eº> <gamma: [0,1;1,0]> <aplica regras: {1;0}> <ordem: [1;5]> <LB 0/1> <delta>\n");
		return (-1);
	}

	strcpy(file_name, argv[1]);
	seed = atoi(argv[2]);
	do_cuts = atoi(argv[3]);
	aplicaE = atoi(argv[4]);
	gamma = strtod(argv[5], NULL);
	string s;
	aplicar_regras = atoi(argv[6]);
	ordem = atoi(argv[7]);
	fazLB = atoi(argv[8]);
	delta = atoi(argv[9]);

#endif

#if hib

	if (argc != 9) {
		printf(
		    "./executavel <data/nome_instancia> <sed> <aplica_cortes> <remove Eº> <gamma: [0,1;1,0]> <aplica regras: {1;0}> <ordem: [1;5]> perfix\n");
		return (-1);
	}

	strcpy(file_name, argv[1]);
	seed = atoi(argv[2]);
	do_cuts = atoi(argv[3]);
	aplicaE = atoi(argv[4]);
	gamma = strtod(argv[5], NULL);
	string s;
	aplicar_regras = atoi(argv[6]);
	ordem = atoi(argv[7]);
	perfix = strtod(argv[8], NULL);
	//fazLB = atoi(argv[8]);
	//delta = atoi(argv[9]);

#endif

//	Pegar qtd de vertices para inicializar estruturas auxiliares
	std::fstream fin(file_name);

	fin >> qtd_n;
	fin >> premio_solicitado;
	fin >> d;
	fin >> valor_otimo_instancia;

//	Matriz das distancias entre os vertices
	Matriz *C = create_matriz(qtd_n, qtd_n);

//	alocação da solucao do grasp
	s_melhor_rota = new int[qtd_n];

//	alocação da matriz BL
	matriz_BL = new int*[qtd_n];

	for (int i = 0; i < qtd_n; i++) {
		matriz_BL[i] = new int[qtd_n];
	}

	Vertice* vertices_instancia = (Vertice*) malloc(qtd_n * sizeof(Vertice));

//	Carregar instancia do arquivo
	get_instancia_prrcp(file_name, vertices_instancia, C, qtd_n, &qtd_T, d,
	                    premio_solicitado, root, qnd_w, qnd_vt, qnd_t, prize, D,
	                    valor_otimo, vertices, matriz_custo, prize_V, R_set);

	int qnd_vt_i = qnd_vt, qnd_t_i = qnd_t, qnd_w_i = qnd_w;

	fin.close();

	/* HM */
	int cont_parada = 0;
	itExecutadas = 0;

	// clock_t start = clock();

	// double custoHib = HM(perfix, cont_parada, solucao, itExecutadas, gamma, pE,
	//                      matriz_BL, vertices_instancia, C, alfa, grasp_max, tempo_max, qtd_T,
	//                      premio_solicitado, seed, d, nome_instancia, premio_solicitado,
	//                      s_melhor_rota, numIt, L, delta, fazLB, aplicaE, root, qnd_vt, qtd_n,
	//                      vertices, matriz_custo, prize, exigir_r, prize_V, valor_otimo,
	//                      limPrimal, subProbSize);

	// clock_t end = clock();
	// double tempoHibrido = (end - start) / (double) CLOCKS_PER_SEC;

	// printf("%s \t %.0f \t \%.2f \n", file_name, custoHib, tempoHibrido);

	// return 1;

//	printSolVector(solucao);

	/*
	 * Regras de redução em GRAFO ORIGINAL
	 */

	if (aplicar_regras) {

//		printf("%s \t", file_name);

		clock_t start = clock();
		aplicarRegras(nReducoes, ordem, false, matriz_custo, vertices, prize,
		              qnd_vt, qnd_t, qnd_w, limPrimal, exigir_r, R_set, matriz_BL, 1);

		//		printf("Depois |R| = %d \t |T| = %d \t |W| = %d \n\n", qnd_vt, qnd_t,
		//				qnd_w);

//		printf("%d \t %d \t %d \t", qnd_vt, qnd_t, qnd_w);

		//
		//		printf("Total de reduções: %d \n\n",
		//						abs(qnd_vt - qnd_vt_i) + abs(qnd_t - qnd_t_i)
		//								+ abs(qnd_w - qnd_w_i));
		clock_t end = clock();
		tempo_regras_original = (end - start) / (double) CLOCKS_PER_SEC;

		rReduzidos = abs(qnd_vt - qnd_vt_i);
		tReduzidos = abs(qnd_t - qnd_t_i);
		wReduzidos = abs(qnd_w - qnd_w_i);
		totalReducoes = rReduzidos + tReduzidos + wReduzidos;

		// printf("%d \t %d \t %d \t %d \t", rReduzidos, tReduzidos, wReduzidos,
		//        totalReducoes);

		//	n reducoes
//		for (int i = 1; i < 15; i++)
//			printf("%d \t", nReducoes[i]);
//		printf("%.2f \n", tempo_regras_original);
	}

//	return 0;

	/*
	 * GRASP
	 */

	if (aplicaE) {

		clock_t start = clock();

		custo_obtido = grasp_basico(vsol, cont_parada, solucao, itExecutadas,
		                            gamma, pE, matriz_BL, vertices_instancia, C, alfa, grasp_max,
		                            tempo_max, qtd_T, premio_solicitado, seed, d, qtd_n,
		                            nome_instancia, premio_solicitado, s_melhor_rota, numIt,
		                            subProbSize);

		clock_t end = clock();

		// printf("GRASP: %s \t %.0f \t %.0f \n", file_name, valor_otimo,
		//        custo_obtido);

//		return 1;

		if (liberar) {
//			liberando arestas de Eº
			liberar_arestas(L, s_melhor_rota, vertices_instancia, matriz_BL,
			                matriz_custo);
		}

		// printf("arestas liberadas: %d \n", L.size());

		tempo_grasp = (end - start) / (double) CLOCKS_PER_SEC;

		limPrimal = custo_obtido; //limite primal para o solver

		for (int i = 0; i < qtd_n; ++i) {

			if (vertices_instancia[i].qtd_cobre > 0) {
				free(vertices_instancia[i].cobre);
			}

			if (vertices_instancia[i].qtd_coberto_por > 0) {
				free(vertices_instancia[i].coberto_por);
			}
		}

		free_matriz(C);
		free(vertices_instancia);

//		contabilização para analise de Eº
//		contRemovidosV(vertices, qtd_n, matriz_BL, qnd_vt, qnd_t);

		/*
		 * Reduzindo o grafo
		 */
		reduzGrafo(vertices, qnd_vt, qtd_n, limPrimal, matriz_custo, matriz_BL);

//		printPython(solucao, vertices, qtd_n, limPrimal);

	}

	/*info do grafo reduzido*/
//	printf("|R| = %d, |T| = %d, |W| = %d \n", qnd_vt, qnd_t, qnd_w);
//	,: ");
	/*
	 * Regras de redução em GRASP ESPARSO
	 */

	if (aplicar_regras) {
		//		printf("Antes |R| = %d \t |T| = %d \t |W| = %d \n", qnd_vt, qnd_t,
		//				qnd_w);

		//		printf("%s \t %d \t %d \t %d \t", file_name, qnd_vt, qnd_t, qnd_w);
//		printf("%s \t", file_name);

		clock_t start = clock();
		aplicarRegras(nReducoes, ordem, false, matriz_custo, vertices, prize,
		              qnd_vt, qnd_t, qnd_w, limPrimal, exigir_r, R_set, matriz_BL, 0);

		//		printf("Depois |R| = %d \t |T| = %d \t |W| = %d \n\n", qnd_vt, qnd_t,
		//				qnd_w);

		//		printf("%d \t %d \t %d \t", qnd_vt, qnd_t, qnd_w);

		//		printf("%d \t",
		//				abs(qnd_vt - qnd_vt_i) + abs(qnd_t - qnd_t_i)
		//						+ abs(qnd_w - qnd_w_i));
		//
		//		printf("Total de reduções: %d \n\n",
		//						abs(qnd_vt - qnd_vt_i) + abs(qnd_t - qnd_t_i)
		//								+ abs(qnd_w - qnd_w_i));
		clock_t end = clock();
		tempo_regras_esparso = (end - start) / (double) CLOCKS_PER_SEC;

		//	n reducoes
//		for (int i = 1; i < 15; i++)
//			printf("%d \t", nReducoes[i]);

		double tempo_total = tempo_regras_esparso + tempo_regras_original;

		rReduzidos = abs(qnd_vt - qnd_vt_i);
		tReduzidos = abs(qnd_t - qnd_t_i);
		wReduzidos = abs(qnd_w - qnd_w_i);
		totalReducoes = rReduzidos + tReduzidos + wReduzidos;
//
		printf("%s \t %d \t %d \t %d \t %d \t", file_name, rReduzidos,
		       tReduzidos, wReduzidos, totalReducoes);
//
		printf("%.2f \n", tempo_total);

	}

	return 1;

	/*
	 * SOLVER
	 */

	clock_t startR = clock();

	double custo_solver = solver(vsol, L, delta, fazLB, solucao, matriz_BL,
	                             aplicaE, root, qnd_t, qnd_vt, qtd_n, vertices, matriz_custo, prize,
	                             exigir_r, prize_V, valor_otimo, limPrimal);

	clock_t endR = clock();

	double tempo_solver = (endR - startR) / (double) CLOCKS_PER_SEC;
	double tempoMPO = 0.0;
	tempoMPO = tempo_solver + tempo_grasp;

#if db
	cout << "instancia: " << file_name << endl;
	cout << "otimo:" << valor_otimo << endl;
	cout << "custo GRASP: " << custo_obtido << endl;
	cout << "pE: " << pE << endl;
	cout << "it executadas: " << itExecutadas << endl;
	cout << "custo SOLVER: " << custo_solver << endl;
	cout << "tempo solver: " << tempo_solver << endl;
#endif

//	printf("%s \t %.0f \t  %.0f \t  %.2f \t %.0f \t %.2f \n", file_name,
//			valor_otimo, custo_obtido, tempo_grasp, custo_solver, tempoMPO);

	// print do local branching
	printf("MPO: %s \t %.0f \t %.0f \t %.2f \n", file_name, valor_otimo,
	       custo_solver, tempoMPO);

//	analisaE0(vertices, arestas, qtd_n, matriz_BL, matriz_custo);

	return (0);

}

double solver(double perfix, vector<Vertex>& vSol, vector<aresta>& L, int delta, int fazLB,
              vector<edgeSol>& solucao, int** matriz_BL, bool removeE0, int root,
              int qnd_t, int qnd_vt, int qtd_n, vector<Vertex> vertices,
              double **matriz_custo, double prize, bool exigir_r, double prize_V,
              double otimo, double limPrimal) {

//	printf("CHAMADA DO SOLVER \n\n");

	double custo_primal = limPrimal + 0.000001;
	bool fixOne = false;

	try {
//		int **solucao_matriz;
//
//		solucao_matriz = new int*[qtd_n];
//
//		for (int i = 0; i < qtd_n; i++) {
//			solucao_matriz[i] = new int[qtd_n];
//		}
//
//		for (int i = 0; i < qtd_n; i++) {
//			for (int j = 0; j < qtd_n; j++) {
//				solucao_matriz[i][j] = 0;
//			}
//		}
//
//		for (int i = 0; i < solucao.size(); i++) {
//			int a = solucao[i].o;
//			int b = solucao[i].d;
//			solucao_matriz[a][b] = 1;
//			solucao_matriz[b][a] = 1;
//		}

		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		env.set(GRB_IntParam_Presolve, 0);
		env.set(GRB_IntParam_Cuts, 0);
		model.getEnv().set(GRB_IntParam_PreCrush, 1); //necessario quando se utiliza seus próprios cortes
		model.getEnv().set(GRB_IntParam_Threads, 1);
		env.set(GRB_DoubleParam_Heuristics, 0.0);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 1800.0);
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.getEnv().set(GRB_IntParam_SolutionLimit, 1);
//		model.getEnv().set(GRB_DoubleParam_Cutoff, custo_primal);

		// Must set LazyConstraints parameter when using lazy constraints
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		GRBVar *y = 0; //representa os vertices
		GRBVar **x = 0; //representa as arestas

		int bound = qtd_n * (qtd_n + 1) / 2;

//		GRBVar delta1, delta2; //variavel delta do lb
//		delta1 = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
//				"delta_1");
//		delta2 = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
//				"delta_2");

		// GRBVar** z = 0;
		GRBLinExpr lhs, rhs, obj;
//		obj = delta1 + delta2;

		char varName[100];

		//cadastrando variavies
		y = new GRBVar[qtd_n];

		for (unsigned i = 0; i < vertices.size(); i++) {
			if (vertices[i].reduzido || vertices[i].tipo == 2) {
				continue;
			}

			// cout << "indice Y: " << indiceY << endl;
			sprintf(varName, "y_%d", i);
			// if (vertices[indiceY].reduzido)
			// 	y[i] = model.addVar(0.0, 0.0, 0.0, GRB_BINARY, varName);
			// else
			y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, varName);
//			cout << "y_" << i << ":" << vertices[i].tipo << "  ";
//			getchar();
		}

		// printf("\n");

		x = new GRBVar*[qtd_n];

		for (int i = 0; i < qtd_n; i++) {
			x[i] = new GRBVar[qtd_n];
		}

		for (int i = 0; i < qtd_n; i++) {
			if (vertices[i].reduzido || vertices[i].tipo == 2) {
				continue;
			}

			for (int j = i + 1; j < qtd_n; j++) {
				if (vertices[j].reduzido || vertices[j].tipo == 2) {
					continue;
				}

				sprintf(varName, "x_%d_%d", i, j);
//				cout<<varName<<endl;
//				cin.get();
				x[i][j] = model.addVar(0.0, 1.0, matriz_custo[i][j],
				                       GRB_BINARY, varName);
			}
		}

		model.update();

		/* ***************************************************
		 *
		 * 						Restrições
		 *
		 * ***************************************************/

		/*
		 * Local Branching
		 */
		if (0) {
			lhs = 0;

			/* Local branching com as arestas liberadas
			 *
			 * for (long unsigned j = 0; j < L.size(); j++) {
			 int a, b;
			 a = L[j].o.id;
			 b = L[j].d.id;
			 if (a < b)
			 lhs += x[a][b];
			 else
			 lhs += x[b][a];
			 }*/

//			model.addConstr(lhs <= delta2, "local_branching_2");
			lhs = 0;

		}

		/**
		 //		  Fixação de variáveis: conjuntos Eº e E*
		 //		 */
//		cout << "remove E0: " << removeE0 << endl;
		if (removeE0) {

			for (int i = 0; i < qtd_n; i++) {
				if (vertices[i].reduzido || vertices[i].tipo == 2) {
					continue;
				}

				for (int j = i + 1; j < qtd_n; j++) {
					if (vertices[j].reduzido || vertices[j].tipo == 2) {
						continue;
					}

					if (matriz_BL[i][j] == 0) {

						lhs = x[i][j];
						model.addConstr(lhs == 0);
						lhs = 0;
					}
				}
			}

		}

		/*
		 * Fixacao variaveis em um
		 * */

		if (perfix > 0) {
			//ordena a solucao por frequencia
			sort(solucao.begin(), solucao.end(), ordenaCustoFreq);

			//		printSolVector(solucao);

			lhs = 0;
			int nFixed = solucao.size() * perfix;
//			printf("arestas fixadas: %d \n", nFixed);
			int a, b, i;

			for (i = 0; i < nFixed; i++) {
				a = solucao[i].o;
				b = solucao[i].d;

//				printf("%d - %d: %0.f\n", a,b,matriz_custo[a][b]);
				if (a < b) {
					lhs = x[a][b];
				} else {
					lhs = x[b][a];
				}

				model.addConstr(lhs == 1);
				lhs = 0;
			}

//			printf("i final: %d\n", i);

		}

		/**
		 * Fim da fixação
		 */

		// restricao obtida pela regra R7:
		// \sum_{i \in VT} yi >= 1
		if (exigir_r) {
			for (unsigned i = 0; i < vertices.size(); i++) {
				if (vertices[i].reduzido || vertices[i].tipo != 0) {
					continue;
				}

				lhs += y[i];
			}

			model.addConstr(lhs >= 1.0);
			lhs = 0;
		}

		//(0)
		lhs = 0;

		for (int i = 0; i < qtd_n; i++) {
			if (vertices[i].tipo != 1) {
				continue;
			}

			lhs += y[i];
		}

		model.addConstr(lhs == qnd_t);
		lhs = 0;

		//(2)
		for (int i = 0; i < qtd_n; i++) {
			if (vertices[i].reduzido || vertices[i].tipo == 2) {
				continue;
			}

			for (int j = i + 1; j < qtd_n; j++) {
				if (vertices[j].reduzido || vertices[j].tipo == 2) {
					continue;
				}

				lhs += x[i][j];
			}
		}

		for (int r = 0; r < qtd_n; r++) {
			if (vertices[r].reduzido || vertices[r].tipo != 0) {
				continue;
			}

			rhs += y[r];
		}

		rhs += qnd_t;
		model.addConstr(lhs == rhs, "2");
		lhs = rhs = 0;

		//(3)
		for (unsigned k = 0; k < vertices.size(); ++k) {
			// cout << "indice: " << indice << endl;
			if (vertices[k].reduzido || vertices[k].tipo == 2) {
				continue;
			}

			lhs += vertices[k].premio * y[k];
		}

		model.addConstr(lhs >= prize, "3");
		lhs = 0;

		// (4) ∑_{k ∈ S_l}yk ≥ 1 (∀l ∈ W )
		for (unsigned l = 0; l < vertices.size(); l++) {
			if (vertices[l].tipo != 2 || vertices[l].reduzido) {
				continue;
			}

			// if (vertices[indiceL].reduzido) continue;
			// cout << "w_" << w << endl;
			if (vertices[l].coberto_por.size() > 0) {
				for (unsigned j = 0; j < vertices[l].coberto_por.size(); j++) {
					int indiceK = vertices[l].coberto_por[j];

					if (vertices[indiceK].reduzido
					        || vertices[indiceK].tipo == 2) {
						continue;
					}

					// cout << "indice K: " << indiceK << endl;

					lhs += y[indiceK];
				}

				// cout << lhs << endl;
				sprintf(varName, "4.%d", l);
				model.addConstr(lhs >= 1.0, varName);
				lhs = 0;
			}
		}

		// (5)
		for (int t = 0; t < qtd_n; t++) {
			if (vertices[t].reduzido || vertices[t].tipo != 1) {
				continue;
			}

			for (int i = 0; i < qtd_n; i++) {
				if (vertices[i].reduzido || vertices[i].tipo == 2) {
					continue;
				}

				for (int j = i + 1; j < qtd_n; j++) {
					if (vertices[j].reduzido || vertices[j].tipo == 2) {
						continue;
					}

					if (vertices[i].id == t || vertices[j].id == t) {
						lhs += x[i][j];
					}
				}
			}

			sprintf(varName, "5.%d", t);
			model.addConstr(lhs == 2.0, varName);
			lhs = 0;
		}

		// (6)
		for (int r = 0; r < qtd_n; r++) {
			if (vertices[r].reduzido || vertices[r].tipo != 0) {
				continue;
			}

			for (int i = 0; i < qtd_n; i++) {
				if (vertices[i].reduzido || vertices[i].tipo == 2) {
					continue;
				}

				for (int j = i + 1; j < qtd_n; j++) {
					if (vertices[j].reduzido || vertices[j].tipo == 2) {
						continue;
					}

					if (vertices[i].id == r || vertices[j].id == r) {
						lhs += x[i][j];
					}
				}
			}

			rhs = 2 * y[r];
			sprintf(varName, "6.%d", r);
			model.addConstr(lhs == rhs, varName);
			lhs = rhs = 0;
		}

		// (10)
		for (int i = 0; i < qtd_n; i++) {
			if (vertices[i].reduzido || vertices[i].tipo == 2) {
				continue;
			}

			for (int j = i + 1; j < qtd_n; j++) {
				if (vertices[j].reduzido || vertices[j].tipo == 2) {
					continue;
				}

				lhs = x[i][j];
				rhs = y[i];
				sprintf(varName, "10.%d.%d.%d", i, j, i);
				model.addConstr(lhs <= rhs, varName);

				rhs = y[j];
				sprintf(varName, "10.%d.%d.%d", i, j, j);
				model.addConstr(lhs <= rhs, varName);

				lhs = rhs = 0;
			}
		}

		/* funcao objetivo diferente para o lb com delta adaptativo
		 * */

//		model.setObjective(obj, GRB_MINIMIZE);
		// Callback
		// Mycallback cb = Mycallback(x, root, qtd_n, matriz_custo, vertices, y,
		//                            &model, prize, prize_V);
		// model.setCallback(&cb);

		model.update();
		// model.write("nao_dir.lp");
		model.optimize();

//		if (model.get(GRB_DoubleAttr_ObjVal) > limPrimal) {
//			printf("maior que o limPrimal \n");
//			return (limPrimal);
//		} else {
		double c = 0.0;
//		solucao.clear();
//		for (int i = 0; i < qtd_n; i++) {
//			if (vertices[i].reduzido || vertices[i].tipo == 2)
//				continue;
////				vSol.push_back(i);
//			for (int j = i + 1; j < qtd_n; j++) {
//				if (vertices[j].reduzido || vertices[j].tipo == 2)
//					continue;
//
//				if (x[i][j].get(GRB_DoubleAttr_X) <= 0.9)
//					continue;
//
////				printf("x_%d_%d=%0.2f \t", i, j, x[i][j].get(GRB_DoubleAttr_X));
//
////				if (x[i][j].get(GRB_DoubleAttr_X) >= 0.60) {
////					GRBLinExpr exp = x[i][j];
////					model.addConstr(exp == 1.0);
////				}
//
//				edgeSol e;
//				e.o = i;
//				e.d = j;
//				e.custo = matriz_custo[i][j];
//				c += matriz_custo[i][j];
//				solucao.push_back(e);
//			}
//		}

		return (model.get(GRB_DoubleAttr_ObjVal));
//		}

//		return (model.get(GRB_DoubleAttr_ObjVal));

	} catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		return (limPrimal);
	} catch (...) {
		cout << "Exception during optimization" << endl;
		return (limPrimal);
	}

//	if (model.get(GRB_DoubleAttr_ObjVal) > limPrimal)
//		return (limPrimal);

}

bool verifDesigTriangular(vector<Vertex> vertices, double **matriz_custo) {
	bool desigualdade_triangular = true;

	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].reduzido || vertices[i].tipo == 2) {
			continue;
		}

		for (int j = 0; j < vertices.size(); j++) {
			if (vertices[j].reduzido || vertices[j].tipo == 2) {
				continue;
			}

			if (i == j) {
				continue;
			}

			for (int k = 0; k < vertices.size(); k++) {
				if (vertices[k].reduzido || vertices[k].tipo == 2) {
					continue;
				}

				if (k == j) {
					continue;
				}

				if (matriz_custo[i][j]
				        > matriz_custo[i][k] + matriz_custo[k][j]) {
					desigualdade_triangular = false;
//					printf("\n(%d_%d)=%0.f > (%d_%d,%d_%d)=%0.f \n\n", i, j,
//							matriz_custo[i][j], i, k, k, j,
//							matriz_custo[i][k] + matriz_custo[k][j]);
					break;
				}

				if (!desigualdade_triangular) {
					break;
				}
			}

			if (!desigualdade_triangular) {
				break;
			}
		}

		if (!desigualdade_triangular) {
			break;
		}
	}

	return desigualdade_triangular;

}

void printPython(double** m, vector<Vertex> v, int n, double limPrimal) {
	printf("\n");

	for (int i = 0; i < n; i++) {
		if (v[i].reduzido || v[i].tipo == 2) {
			continue;
		}

		for (int j = i + 1; j < n; j++) {
			if (v[j].reduzido || v[j].tipo == 2) {
				continue;
			}

			if (m[i][j] < limPrimal) {
				printf("(\"%d(%d)\", \"%d(%d)\", \"1\"),\t", i, v[i].tipo, j,
				       v[j].tipo, m[i][j]);
			}

		}
	}
}

void get_solucao(int *solucao, Vertice* vertices, vector<int> *S) {
//printf("\n*** Print da Solucao ***\n");
	int p = 0;

	while (solucao[p] == -1) {
		p++;
	}

//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = solucao[antecessor];
//	printf("\n\n Arestas: %d [%s] --> %d[%s]", antecessor,
//			vertices[antecessor].tipo, sucessor, vertices[sucessor].tipo);
	int count = 1;
	S->push_back(antecessor);
	S->push_back(sucessor);

	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = solucao[sucessor];
//		printf(" --> %d [%s]", sucessor, vertices[sucessor].tipo);
		S->push_back(sucessor);
		count++;
//		if (count > 100)
//			scanf("%d", &count);

	}

//	cout << "numero de arestas na sol: " << S->size() << endl;

//	S->push_back(solucao[1]);
}

void liberar_arestas(vector<aresta>& L, int *solucao, Vertice* vertices,
                     int**& matriz_BL, double** matriz_custo) {
	bool debug = false;
	vector<int> S;
	get_solucao(solucao, vertices, &S);

//	printf("\n\n");

	if (debug) {
		printf("\nLiberar arestas: \n\n");
	}

	unsigned tam = S.size();
	int a, b, c;
	double custoS = 0.0;

	for (unsigned i = 1; i < tam - 1; i++) {
		custoS = 0.0;
		a = b = c = 0;

		a = S[i - 1];
		b = S[i];
		c = S[i + 1];

		if (matriz_BL[a][c] == 1 or matriz_BL[c][a] == 1) {
			continue;
		}

		//			printf("%d -> %d -> %d \n", a, b, c);
		custoS = matriz_custo[a][b] + matriz_custo[b][c];

		if (custoS > matriz_custo[a][c]) {
			if (debug) {
				printf(
				    "%d-%d: %.0f é mais barata que %d-%d-%d: %0.f, e será liberada.\n",
				    a, c, matriz_custo[a][c], a, b, c, custoS);
			}

			matriz_BL[a][c] = 1;
			matriz_BL[c][a] = 1;
			aresta e;
			no origem, destino;
			origem.id = a;
			destino.id = c;
			e.o = origem;
			e.d = destino;
			e.custo = matriz_custo[a][c];
			L.push_back(e);
		} else if (custoS > matriz_custo[c][a]) {
			if (debug) {
				printf(
				    "%d-%d: %.0f é mais barata que %d-%d-%d: %0.f, e será liberada.\n",
				    a, c, matriz_custo[a][c], a, b, c, custoS);
			}

			matriz_BL[a][c] = 1;
			matriz_BL[c][a] = 1;
			aresta e;
			no origem, destino;
			origem.id = a;
			destino.id = c;
			e.o = origem;
			e.d = destino;
			e.custo = matriz_custo[a][c];
			L.push_back(e);
		}
	}

	a = b = c = 0;
	custoS = 0.0;

	a = S[tam - 2];
	b = S[0];
	c = S[1];
	custoS = matriz_custo[a][b] + matriz_custo[b][c];

	if (matriz_BL[a][c] == 0 and matriz_BL[c][a] == 0) {

		if (custoS > matriz_custo[a][c]) {
			if (debug) {
				printf(
				    "%d-%d: %.0f é mais barata que %d-%d-%d: %0.f, e será liberada.\n",
				    a, c, matriz_custo[a][c], a, b, c, custoS);
			}

			matriz_BL[a][c] = 1;
			matriz_BL[c][a] = 1;
			aresta e;
			no origem, destino;
			origem.id = a;
			destino.id = c;
			e.o = origem;
			e.d = destino;
			e.custo = matriz_custo[a][c];
			L.push_back(e);
		} else if (custoS > matriz_custo[c][a]) {
			if (debug) {
				printf(
				    "%d-%d: %.0f é mais barata que %d-%d-%d: %0.f, e será liberada.\n",
				    a, c, matriz_custo[a][c], a, b, c, custoS);
			}

			matriz_BL[a][c] = 1;
			matriz_BL[c][a] = 1;
			aresta e;
			no origem, destino;
			origem.id = a;
			destino.id = c;
			e.o = origem;
			e.d = destino;
			e.custo = matriz_custo[a][c];
			L.push_back(e);
		}

	}
}

int aplicarRegras(int*& nReducoes, int ordem, bool debug, double** matriz_custo,
                  vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t,
                  int &qnd_w, double limPrimal, bool &exigir_r, vector<int> R_set,
                  int ** matriz_BL, int grafoOriginal) {
//	printf("\nimpacto das regras DEPOIS do grasp: \n");
	bool dt = false;
	dt = verifDesigTriangular(vertices, matriz_custo);
//	int nReducoes[15];

	if (dt and debug) {
		printf("\nApplying rules r3 to r5 \n");
	}

	if (debug) {
		printf("Reduction rules applied\n");
		cout << "\nANTES:" << endl;
		printf("\nAntes |R| = %d \t |T| = %d \t |W| = %d \n", qnd_vt, qnd_t,
		       qnd_w);
	}

	switch (ordem) {
	case 1:
		//Order 1: r6, r2, r13, r8, r5, r9, r7, r3, r11, r10, r4, r1, r14, r12.
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);
//		printf("regra 6 aplicada\n");
		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);
//		printf("regra 2 aplicada\n");
		nReducoes[13] = regra_13(vertices, qnd_w);
//		printf("regra 13 aplicada\n");
		nReducoes[14] = regra_R8(vertices, prize); //dado que sempre |T|>1 entao R8 é sempre ativada

//		printf("regra 14 aplicada\n");
		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
//			printf("regra 5 aplicada\n");
		}

		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);
//		printf("regra 9 aplicada\n");
		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
		nReducoes[7] = exigir_r;

//		printf("regra 7 aplicada\n");
		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
//			printf("regra 3 aplicada\n");
		}

		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);
//		printf("regra 11 aplicada\n");
		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal); //r10

//		printf("regra 10 aplicada\n");
		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);

//			printf("regra 4 aplicada\n");
		}

		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);
//		printf("regra 1 aplicada\n");
		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
//		printf("regra 14 aplicada\n");
		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
//		printf("regra 12 aplicada\n");
		break;

	case 2:
		//Order 2: r6, r5, r2, r13, r8, r7, r9, r4, r14, r1, r10, r11, r3, r12.
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);

		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
		}

		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);
		nReducoes[13] = regra_13(vertices, qnd_w);
		nReducoes[8] = regra_R8(vertices, prize);
		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
		nReducoes[7] = exigir_r;
		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);

		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);
		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal);
		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);

		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
		break;

	case 3:
		//sequência 3: r10, r11,  r12, r1, r3, r9, r2, r4, r6,r7,r13, r5, r14, r8.
		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal);
//		printf("regra 10 aplicada\n");
		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);
//		printf("regra 11 aplicada\n");
		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
//		printf("regra 12 aplicada\n");
		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);

//		printf("regra 1 aplicada\n");
		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
//			printf("regra 3 aplicada\n");
		}

		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);
//		printf("regra 9 aplicada\n");
		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);

//		printf("regra 2 aplicada\n");
		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
//			printf("regra 4 aplicada\n");
		}

		nReducoes[8] = regra_R8(vertices, prize);
//		printf("regra 8 aplicada\n");
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);
//		printf("regra 6 aplicada\n");
		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
//		printf("regra 7 aplicada\n");
		nReducoes[7] = exigir_r;
		nReducoes[13] = regra_13(vertices, qnd_w);

//		printf("regra 13 aplicada\n");
		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
//			printf("regra 5 aplicada\n");
		}

		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
//		printf("regra 14 aplicada\n");
		break;

	case 4:
		//sequência 4: r11, r3, r10, r12, r1, r4, r9, r6, r2, r8, r7, r5, r14, r13
		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);

		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal);
		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);

		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);
		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);
		nReducoes[8] = regra_R8(vertices, prize);
		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
		nReducoes[7] = exigir_r;

		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
		}

		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
		nReducoes[13] = regra_13(vertices, qnd_w);
		break;

	case 5:
		//sequência 5: r2, r11, r3, r12, r1, r4, r10, r9, r6,  r8, r7, r5, r14, r13
		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);

		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);

		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
		nReducoes[7] = exigir_r;

		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);

		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal);
		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);

		nReducoes[8] = regra_R8(vertices, prize);

		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
		}

		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
		nReducoes[13] = regra_13(vertices, qnd_w);
		break;

	case 6: //r6,r10, r2, r11, r13, r12, r8, r1, r5, r3, r9, r4, r7, r14 (combinação de r1 e r3).
		nReducoes[6] = regra_R6(vertices, prize, qnd_vt, qnd_t, qnd_w);
		nReducoes[10] = regra_10(vertices, prize, qnd_vt, qnd_t, qnd_w,
		                         matriz_custo, limPrimal, grafoOriginal);
		nReducoes[2] = regra_R2(vertices, prize, qnd_vt, qnd_t);
		nReducoes[11] = regra_11(vertices, qnd_vt, qnd_t, qnd_w, matriz_BL,
		                         grafoOriginal);
		nReducoes[13] = regra_13(vertices, qnd_w);
		nReducoes[12] = regra_12(vertices, qnd_vt, qnd_t, qnd_w, grafoOriginal);
		nReducoes[8] = regra_R8(vertices, prize);
		nReducoes[1] = regra_R1(vertices, prize, qnd_vt, qnd_t);

		if (dt) {
			nReducoes[5] = regra_R5(vertices, prize, qnd_vt, qnd_t,
			                        matriz_custo);
		}

		if (dt) {
			nReducoes[3] = regra_R3(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		nReducoes[9] = regra_R9(vertices, prize, matriz_custo, qnd_vt);

		if (dt) {
			nReducoes[4] = regra_R4(vertices, prize, qnd_vt, qnd_t, qnd_w,
			                        grafoOriginal);
		}

		exigir_r = regra_R7(vertices, prize, qnd_vt, qnd_t, qnd_w, exigir_r,
		                    R_set);
		nReducoes[7] = exigir_r;
		nReducoes[14] = regra_14(vertices, qnd_vt, matriz_custo, prize, qnd_t,
		                         qnd_w, grafoOriginal);
		break;

	default:
		printf("Opção inválida de regra\n");
		return -1;

	}

	/*
	 * FIM   Regras de redução
	 */
//	if (debug) {
//		printf("|R| = %d \t |T| = %d \t |W| = %d \n\n", qnd_vt, qnd_t, qnd_w);
//	}
//	if (debug) {
//		cout << "\nDEPOIS: \n";
////		printf("|E| after = %d \n", arestas_reais - contE0);
//		printf("After |R| = %d \t |T| = %d \t |W| = %d \n\n", qnd_vt, qnd_t,
//				qnd_w);
//
//		printf("Ativação: \n");
//		for (int i = 1; i < 15; i++) {
//			cout << i << " : " << ativadas[i] << endl;
//		}
//		cout << endl;
//	}
//	printf("\n-%s \t %d \t %d \t %d \n", file_name, qnd_vt, qnd_t, qnd_w);
	return 0;
}

void contRemovidosV(vector<Vertex> vertices, int qtd_n, int** matriz_BL,
                    int qnd_vt, int qnd_t) {
	vector<int> verticesRemovidos;

	for (int i = 0; i < qtd_n; i++) {
		if (vertices[i].tipo == 2 or vertices[i].reduzido) {
			continue;
		}

		int somaIncid = 0;

		for (int j = i + 1; j < qtd_n; j++) {
			if (vertices[j].tipo == 2 or vertices[j].reduzido or i == j) {
				continue;
			}

			somaIncid += matriz_BL[i][j];
		}

		if (somaIncid == 0) {
			verticesRemovidos.push_back(i);
		}
	}

//	printf("%d vertices removidos de %d (VT + T)\n", verticesRemovidos.size(),
//			(qnd_vt + qnd_t));
//	printf("v_removido[0] = %d \n", verticesRemovidos[0]);
}

void reduzGrafo(vector<Vertex>& vertices, int& qnd_vt, int& qtd_n,
                double limPrimal, double** matriz_custo, int** matriz_BL) {
	for (int i = 0; i < qtd_n; i++) {
		if (vertices[i].tipo != 0) {
			continue;
		}

		double somaAdja = 0.0;
		int cont = 0;

		for (int j = 0; j < qtd_n; j++) {
			if (vertices[j].tipo == 2 or i == j) {
				continue;
			}

			somaAdja += matriz_custo[i][j];
			cont++;
		}

		//			printf("somaAdja de %d = %f \n", i, somaAdja);
		//			printf("comp = %f \n", limPrimal * cont );
		//			getchar();
		if (somaAdja == limPrimal * cont) {
			vertices[i].reduzido = true;
			qnd_vt--;
		}
	}

//definindo custos infinitos para arestas que nao pertencem mais a E
	for (int i = 0; i < qtd_n; ++i) {
		for (int j = i + 1; j < qtd_n; ++j) {
			if (matriz_BL[i][j] == 0) {
				matriz_custo[i][j] = matriz_custo[j][i] = limPrimal;
			}

		}
	}

//setting vertices that does not belong to the graph
	for (int i = 0; i < qtd_n; i++) {
		if (vertices[i].tipo != 0) {
			continue;
		}

		double somaAdja = 0.0;
		int cont = 0;

		for (int j = 0; j < qtd_n; j++) {
			if (vertices[j].tipo == 2 or i == j) {
				continue;
			}

			somaAdja += matriz_custo[i][j];
			cont++;
		}

		//			printf("somaAdja de %d = %f \n", i, somaAdja);
		//			printf("comp = %f \n", limPrimal * cont );
		//			getchar();
		if (somaAdja == limPrimal * cont) {
			vertices[i].reduzido = true;
			qnd_vt--;
		}
	}

}
//
//void analisaE0(vector<Vertex> vertices, vector<aresta>& arestas, int qtd_n,
//		int** matriz_BL, double** matriz_custo) {
////colocando as arestas de Eº num vector
//	for (int i = 0; i < qtd_n; i++) {
//		for (int j = i + 1; j < qtd_n; j++) {
//			if (vertices[i].tipo == 2 or vertices[j].tipo == 2)
//				continue;
//			if (matriz_BL[i][j] > 0)
//				continue;
//			//					if (i == 2 and j == 113) {
//			//						cin.get();
//			//						cout << i << endl;
//			//						cout << j << endl;
//			//					}
//			aresta e;
//			no origem, destino;
//			origem.id = i;
//			destino.id = j;
//			e.o = origem;
//			e.d = destino;
//			e.custo = matriz_custo[i][j];
//			arestas.push_back(e);
//		}
//	}
//
//	printf("tamanho de Eº = %d. \n", arestas.size());
//
////ordenas arestas de Eº de acordo com o custo
//	sort(arestas.begin(), arestas.end(), myfunctionAresta);
//
//	cout << "aresta mais cara de S*: \n";
//	int pos = 99999999, e_i, e_j;
//	while (cin >> e_i >> e_j) {
//		pos = 99999999;
//		cout << " arestas buscadas: " << e_i << "-" << e_j << endl;
//		for (int i = 0; i < arestas.size(); i++) {
//			printf("%d. %d-%d : %.0f \n", i, arestas[i].o.id, arestas[i].d.id,
//					arestas[i].custo);
//			if ((arestas[i].o.id == e_i and arestas[i].d.id == e_j)
//					or (arestas[i].o.id == e_j and arestas[i].d.id == e_i)) {
//				printf("%d. %d-%d : %.0f \n", i, arestas[i].o.id,
//						arestas[i].d.id, arestas[i].custo);
//				pos = i;
//				break;
//			}
//		}
//
//		if (pos == 99999999) {
//			printf("Aresta nao pertence a Eº \n");
//		} else {
//			cout << "pos: " << pos << endl;
//			cout << "tam de Eº: " << arestas.size() << endl;
//
//			printf("%f das arestas mais custosas de Eº nao prestam \n",
//					((double) arestas.size() - (double) pos)
//							/ (double) arestas.size());
//		}
//	}
//
//}

//the hybrid method

double HM(double perfix, int &cont_parada, vector<edgeSol>& sol, int &itExecutadas,
          double gamma, double &pE, int **& matriz_BL, Vertice *vertices_grasp,
          Matriz *C, float alfa, int grasp_max, float tempo_max, int qtd_t,
          int prize, int seed, float d, char *nome_instancia,
          int premio_solicitado, int *&s_melhor_rota, int numIt,
          vector<aresta>& L, int delta, int fazLB, bool removeE0, int root,
          int qnd_vt, int qtd_n, vector<Vertex> vertices, double **matriz_custo,
          double prize_solver, bool exigir_r, double prize_V, double otimo,
          double limPrimal, int* cont) {

	int aplicaE = 1;
	int liberar = 1;
	double custo_grasp_real = 0.0, custo_grasp = 0.0, custo_solver = 0.0;
	bool solverChamado = false;
	vector<edgeSol> arestasFreq;

	vector<Vertex> vsol;

//	while (itExecutadas < numIt) {

	/*chamando o grasp*/
	custo_grasp = grasp_basico(vsol, cont_parada, arestasFreq, itExecutadas,
	                           gamma, pE, matriz_BL, vertices_grasp, C, alfa, grasp_max, tempo_max,
	                           qtd_t, premio_solicitado, seed, d, qtd_n, nome_instancia,
	                           premio_solicitado, s_melhor_rota, numIt, cont);

	custo_grasp_real = custo_grasp;

#ifdef db
//	printSolVector(arestasFreq);
//	printf("otimo: %0.f \n", otimo);
//	printf("solucao do grasp = %.0f \n", custo_grasp);
//	printf("tam da sol: %d \n", arestasFreq.size());

//	checarSolucao(vertices_grasp, C, s_melhor_rota, prize, d, qtd_n);

//	printf("num it = %d \n", itExecutadas);

#endif

	/*liberando arestas de Eº*/
	if (liberar) {
		//			liberando arestas de Eº
		liberar_arestas(L, s_melhor_rota, vertices_grasp, matriz_BL,
		                matriz_custo);
#ifdef db
//		printf("liberando arestas \n");
#endif
	}

	/*chamando o solver*/
	limPrimal = custo_grasp_real + 0.01;
	custo_solver = solver(perfix, vsol, L, delta, fazLB, arestasFreq, matriz_BL,
	                      aplicaE, root, qtd_t, qnd_vt, qtd_n, vertices, matriz_custo,
	                      prize_solver, exigir_r, prize_V, otimo, limPrimal);

//	std::sort(vsol.begin(), vsol.end(), biggerPrize);

//		for (int i = 0; i < vsol.size(); i++) {
//			cout << vsol[i].id << ": " << vertices[vsol[i].id].premio << endl;
//		}
//
//		cin.get();

#ifdef db
//	printf("solucao do solver = %0.f \n", custo_solver);
//	printSolVector(arestasFreq);
#endif

	/*chamando o grasp  novamente*/

//	vsol.clear(); //desativando o bias
//	gamma = 1.0;
//	custo_grasp = grasp_basico(vsol, cont_parada, sol, itExecutadas, gamma, pE,
//			matriz_BL, vertices_grasp, C, alfa, grasp_max, tempo_max, qtd_t,
//			premio_solicitado, seed, d, qtd_n, nome_instancia,
//			premio_solicitado, s_melhor_rota, numIt, cont);
//
//	custo_grasp_real = custo_grasp;
//	printf("solucao do grasp = %.0f \n", custo_grasp);
//	printf("num it = %d \n", itExecutadas);
//	}
	if (custo_solver < custo_grasp) {
		return custo_solver;
	}
	else {
		return custo_grasp;
	}

}

bool myfunctionEdgeSol(edgeSol a, edgeSol b) {
	return (a.custo < b.custo);
}

bool ordenaCustoFreq(edgeSol a, edgeSol b) {
//	if (a.custo == b.custo) {
//		return a.frequencia > b.frequencia;
//	}
	return (a.custo < b.custo);
}

//void solutionCost(vector<edgeSol> sol) {
//	double totalCost=0.0;
//	for(int i=0;i<sol.size();i++){
//		totalCost += sol[i].custo;
//	}
//	printf("custo: ")
//}

