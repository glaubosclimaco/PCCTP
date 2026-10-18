#include <stdio.h>
#include <stdlib.h>

typedef struct matriz{
	int linha;
	int coluna;
	float** v;
}Matriz;

Matriz* create_matriz(int linhas, int colunas);
float get(Matriz* mat, int i, int j);
void set(Matriz* mat, int i, int j, float valor);


void teste() {


	Matriz matriz = *create_matriz(5, 10);

	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 10; ++j) {
			set(&matriz, i,j,i+j);
		}
	}

	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 10; ++j) {
			printf("Valor %d x %d = %f \n",i,j,get(&matriz, i,j));
		}
	}
}
