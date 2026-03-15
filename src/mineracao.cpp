/*
 * mineracao.cpp
 *
 *  Created on: 09/03/2013
 *      Author: rogerio
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "mineracao.h"

using namespace std;

/* função de comparação de Componentes Conexas */
int compara_cc(const void* p1, const void* p2) {
	/* converte ponteiros genéricos para ponteiros de ComponenteConexa */
	ComponenteConexa *cc1 = (ComponenteConexa*) p1;
	ComponenteConexa *cc2 = (ComponenteConexa*) p2;
	/* dados os ponteiros de float, faz a comparação */
	if (cc1->tamanho < cc2->tamanho)
		return -1;
	else if (cc1->tamanho > cc2->tamanho)
		return 1;
	else
		return 0;
}

//Mapear solucoes elite praa arquivo bd.txt
void encode(Matriz* conjunto_elite, int elite_size, int qtd_n) {

	int primeiro = -1;
	int vertice = -1;
	int sucessor_vertice;

	std::ofstream in("bd.txt");
	if (!in) {
		printf("Erro ao abrir arquivo.");
		return;
	}

	for (int i = 0; i < elite_size; ++i) {

		//pegar primeiro vértice da solucao
		primeiro = 0;
		while ((get(conjunto_elite, i, primeiro)) == -1)
			primeiro++;

		//gravar areastar codificadas para o arquivo
		vertice = primeiro;
		sucessor_vertice = (int) get(conjunto_elite, i, vertice);
		std::cout << "Elite " << i;
		do {
			//in << "(" << vertice <<"-"<< sucessor_vertice << ") ";
			std::cout << "(" << vertice << "-" << sucessor_vertice << ") ";
			in << vertice * qtd_n + sucessor_vertice << " ";
			vertice = sucessor_vertice;
			sucessor_vertice = (int) get(conjunto_elite, i, vertice);
		} while (primeiro != vertice);
		in << "\n";

	}
	in.close();
	printf("\nEncode arestas ok.\n");
}

int decode(int qtd_n, DMPadrao *&padroes) {

	FILE *fp = fopen("padroes.txt", "r");
	if (!fp)
		exit(1);

	int r, supp, tam;
	r = fscanf(fp, "%d;%d;", &tam, &supp);

	//cout << "suporte:" <<supp << " ";
	//cout << "tam:" << tam << " ";
	//cout << endl;
	int qtd_padroes = 10;
	int p = 0;
	//Criar lista de Padroes
	//DMPadrao *
	padroes = (DMPadrao *) malloc(qtd_padroes * sizeof(DMPadrao));

	//Enquanto houver padrões
	while (r == 2) {
		int elem;
		padroes[p].id = p;
		padroes[p].suporte = supp;
		padroes[p].tamanho = tam;
		padroes[p].arestas = (Aresta*) malloc(
				padroes[p].tamanho * sizeof(Aresta));
		padroes[p].qtd_cc = 0;

		for (int i = 0; i < tam; i++) {
			int l = fscanf(fp, "%d", &elem);
			if (l < 1) {
				cout << "Erro fscanf" << endl;
				exit(1);
			}
			int origem = elem / qtd_n;
			int destino = elem - (qtd_n * origem);
			padroes[p].arestas[i].id = i;
			padroes[p].arestas[i].vertice_origem = origem;
			padroes[p].arestas[i].vertice_destino = destino;
			padroes[p].arestas[i].analisada = 0;
			//cout << origem << "->" << destino << " , ";
		}

		r = fscanf(fp, "%d;%d;", &tam, &supp);
		p++;

	}

	fclose(fp);
	qtd_padroes = p;
	printf("Decode das arestas ok! \n");
	printf("Qtd de Padrões: %d: \n\n", qtd_padroes);

	for (int i = 0; i < qtd_padroes; ++i) {
		printf("Padrão %d - qtd_arestas: %d - suporte: %d \n", i + 1,
				padroes[i].tamanho, padroes[i].suporte);
		printf("\tArestas: ");
		for (int a = 0; a < padroes[i].tamanho; ++a) {
			printf(" [%d]%d->%d", padroes[i].arestas[a].id,
					padroes[i].arestas[a].vertice_origem,
					padroes[i].arestas[a].vertice_destino);
		}
		printf("\n\n");

		//Pegar componentes conexas
		int aresta_inicial_atual_id = 0; //inicio de cc
		int cc_atual;
		int qtd_aresta_analisadas = 0;
		int qtd_aresta_na_cc = 1;
		int novo_cc = 0;

		while (qtd_aresta_analisadas < padroes[i].tamanho) {
			aresta_inicial_atual_id = 0;
			while (padroes[i].arestas[aresta_inicial_atual_id].analisada == 1) {
				aresta_inicial_atual_id++;
				continue;
			}

			Aresta aresta_atual = padroes[i].arestas[aresta_inicial_atual_id];
			Aresta aresta_antecessora = in_padrao2(aresta_atual.vertice_origem,
					padroes[i]);

			int contador = 0;

			/*Dada um aresta verificar qual a primeira aresta da CC que ela pertence*/
			while (aresta_antecessora.vertice_destino != -1) {
				printf("Aresta Atual: %d\n", aresta_atual.id);
				//printf("Aresta antecessora: %d\n",aresta_antecessora.id);
				aresta_inicial_atual_id = aresta_antecessora.id;
				aresta_atual = aresta_antecessora;
				aresta_antecessora = in_padrao2(aresta_atual.vertice_origem,
						padroes[i]);
				contador++;
				if (contador >= padroes[i].tamanho + 5) {
					printf(
							"Provavelmente este padrão é uma solução elite...!!!!\n");
					//scanf("%d", &contador);
					//exit(0);
					break;
				}
			}

			padroes[i].arestas[aresta_inicial_atual_id].analisada = 1;
			qtd_aresta_analisadas++;

			Aresta proxima_aresta = in_padrao(aresta_atual.vertice_destino,
					padroes[i]);

			while (aresta_atual.vertice_destino == proxima_aresta.vertice_origem) {

				//Primeira CC
				if (padroes[i].qtd_cc == 0) {

					cc_atual = 0;

					//Inicializado o conjunto de componentes conexas com 1und
					padroes[i].ccs = (ComponenteConexa *) malloc(
							sizeof(ComponenteConexa));
					padroes[i].qtd_cc = 1;

					//Inicializando a conjunto de areastas da recem criada CC, com 2und
					padroes[i].ccs[cc_atual].aresta_id = (int *) malloc(
							2 * sizeof(int));
					padroes[i].ccs[cc_atual].tamanho = 2;

					//Adicionar primeiras duas arestas na CC recém criada.
					padroes[i].ccs[cc_atual].aresta_id[0] = aresta_atual.id;
					padroes[i].ccs[cc_atual].aresta_id[1] = proxima_aresta.id;

				}

				if (novo_cc == 1) { //significa nova CC

					novo_cc = 0;
					cc_atual++;

					padroes[i].qtd_cc++;
					//Adicionando espaço para mais 1cc
					padroes[i].ccs = (ComponenteConexa *) realloc(
							padroes[i].ccs,
							padroes[i].qtd_cc * sizeof(ComponenteConexa));

					//Inicializando a conjunto de areastas da recem criada CC, com 2und
					padroes[i].ccs[cc_atual].aresta_id = (int *) malloc(
							2 * sizeof(int));
					padroes[i].ccs[cc_atual].tamanho = 2;

					//Adicionar primeiras duas arestas na CC recém criada.
					padroes[i].ccs[cc_atual].aresta_id[0] = aresta_atual.id;
					padroes[i].ccs[cc_atual].aresta_id[1] = proxima_aresta.id;

				} else if (qtd_aresta_na_cc != 1) { // Significa nova aresta na CC Atual

					//Adicionando espaço para uma aresta na cc atual
					padroes[i].ccs[cc_atual].aresta_id = (int *) realloc(
							padroes[i].ccs[cc_atual].aresta_id,
							(padroes[i].ccs[cc_atual].tamanho + 1)
									* sizeof(int));

					//Adicionar nova aresta na CC atual
					padroes[i].ccs[cc_atual].aresta_id[padroes[i].ccs[cc_atual].tamanho] =
							proxima_aresta.id;

					//Atualiza a quantidade de arestas na CC atual
					padroes[i].ccs[cc_atual].tamanho++;

				}

				//atualiza o contador de aresta na cc atual
				qtd_aresta_na_cc++;

				padroes[i].arestas[proxima_aresta.id].analisada = 1;
				aresta_atual = proxima_aresta;
				proxima_aresta = in_padrao(aresta_atual.vertice_destino,
						padroes[i]);
				qtd_aresta_analisadas++;

				///printf("Atual: %d (ana=%d) e Proxima: %d (ana=%d) \n", aresta_atual.id, aresta_atual.analisada, proxima_aresta.id, proxima_aresta.analisada);

				//printf("Qtd de Arestas Analisadas: %d \n", qtd_aresta_analisadas);

			} //fim while

			if (qtd_aresta_na_cc > 1) {
				qtd_aresta_na_cc = 1;
				novo_cc = 1;
			}

		}

		//Settar a ultima CC como a atual a ser usada
		padroes[i].cc_atual = padroes[i].qtd_cc - 1;

		//Ordenar CC
		qsort(padroes[i].ccs, padroes[i].qtd_cc, sizeof(ComponenteConexa),
				compara_cc);

		//Marcar Excluidas
		/*
		 int qtd_excluidas = 0;
		 for (int c = 0; c < padroes[i].qtd_cc; ++c) {
		 printf("CC %d, qtd_arestas=%d --| %d",c, padroes[i].ccs[c].tamanho, padroes[i].ccs[c].aresta_id[0]);
		 for (int a = 1; a < padroes[i].ccs[c].tamanho; ++a) {
		 printf(" > %d", padroes[i].ccs[c].aresta_id[a]);
		 }

		 Aresta primeira_aresta = padroes[i].arestas[padroes[i].ccs[c].aresta_id[0]];
		 int indice_outra_cc = in_componente_conexa2(primeira_aresta.vertice_origem, padroes[i], c);

		 if ( indice_outra_cc != -1){
		 //Ou seja nao adicionar a temp
		 printf(" APAGAR ESTA!!, pois é sub-segmento da CC = %d \n ", indice_outra_cc);
		 padroes[i].ccs[c].excluida = 1;
		 qtd_excluidas++;
		 }else {
		 padroes[i].ccs[c].excluida = 0;
		 }

		 printf("\n");
		 }
		 */

		//Imprimir
		for (int c = 0; c < padroes[i].qtd_cc; ++c) {
			printf("CC %d, qtd_arestas=%d --| %d", c, padroes[i].ccs[c].tamanho,
					padroes[i].ccs[c].aresta_id[0]);
			for (int a = 1; a < padroes[i].ccs[c].tamanho; ++a) {
				printf(" > %d", padroes[i].ccs[c].aresta_id[a]);
			}
			printf("\n");

		}
	}

	return qtd_padroes;
}

int mining(Matriz* conjunto_elite, int elite_size, int qtd_n,
		DMPadrao *&padroes) {

	//Codificar as aresta da solucao elite
	encode(conjunto_elite, elite_size, qtd_n);

	//LINUX
	//Executar FPMAX (./fpmax seed temp bd_file bd_size suporte qtd_padroes padroes_file)
	system("./fpmax_hnmp 1 100 bd.txt 10 3 10 padroes.txt");

	//WINDOWS
	//system("eclat.exe -s20 -tm bd.txt padroes_20.txt");
	//system("eclat.exe -s30 -tm bd.txt padroes_30.txt");
	//system("python arquivo.py");

	//Decoficar arestas mais frequente mineradas
	return decode(qtd_n, padroes);
}

Aresta in_padrao(int vertice_v, DMPadrao padrao) {
	//Procura no padrao uma aresta com vertice_v seja origem

	//Implementar método de busca binaria
	for (int i = 0; i < padrao.tamanho; i++) {
		if (padrao.arestas[i].vertice_origem == vertice_v)
			return padrao.arestas[i];
	}
	Aresta a;
	a.vertice_origem = -1;
	a.vertice_destino = -1;
	return a;

}

Aresta in_padrao2(int vertice_v, DMPadrao padrao) {
	//Procura no padrao uma aresta com vertice_v seja destino

	//Implementar método de busca binaria
	for (int i = 0; i < padrao.tamanho; i++) {
		if (padrao.arestas[i].vertice_destino == vertice_v)
			return padrao.arestas[i];
	}
	Aresta a;
	a.vertice_origem = -1;
	a.vertice_destino = -1;
	return a;

}

int in_componente_conexa(int vertice_v, DMPadrao padrao, ComponenteConexa &cc) {
	/* Recebu um vértice e verifica se ele faz parte de um componente conexa
	 * retorna -1 caso não esteja em algum CC ou seja o último vertice de alguma
	 * retorna por referencia a CC a qual o vértice_v está contido
	 * retorna a valor inteiro representando o indice da aresta na CC
	 */
	for (int i = 0; i < padrao.qtd_cc; ++i) {
		for (int j = 0; j < padrao.ccs[i].tamanho; ++j) {
			int aresta_id = padrao.ccs[i].aresta_id[j];
			if (padrao.arestas[aresta_id].vertice_origem == vertice_v) {
				cc = padrao.ccs[i];
				return j;
			}
		}
	}
	return -1;
}

int in_componente_conexa2(int vertice_v, DMPadrao padrao, int cc_index) {
	/* Recebu um vértice e verifica se ele faz parte de outra componente conexa
	 * retorna -1 ou o indice da CC ao qual ele tb pertence
	 */
	for (int i = 0; i < padrao.qtd_cc; ++i) {
		for (int j = 0; j < padrao.ccs[i].tamanho; ++j) {
			int aresta_id = padrao.ccs[i].aresta_id[j];
			if ((padrao.arestas[aresta_id].vertice_origem == vertice_v)
					and (i != cc_index)) {
				return i;
			}
		}
	}
	return -1;
}

