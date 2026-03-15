/*
 * Regras.h - Declaracoes das regras de reducao (implementacoes ausentes no repo)
 */

#ifndef REGRAS_H_
#define REGRAS_H_

#include "Graph.h"
#include <vector>

/* Regras R1-R9 (retornam int, exceto R7 que retorna bool) */
int regra_R1(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t);
int regra_R2(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t);
int regra_R3(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, int &qnd_w, int grafoOriginal);
int regra_R4(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, int &qnd_w, int grafoOriginal);
int regra_R5(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, double **matriz_custo);
int regra_R6(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, int &qnd_w);
bool regra_R7(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, int &qnd_w, bool &exigir_r, std::vector<int> R_set);
int regra_R8(std::vector<Vertex> &vertices, double &prize);
int regra_R9(std::vector<Vertex> &vertices, double &prize, double **matriz_custo, int &qnd_vt);

/* Regras 10-14 */
int regra_10(std::vector<Vertex> &vertices, double &prize, int &qnd_vt, int &qnd_t, int &qnd_w, double **matriz_custo, double limPrimal, int grafoOriginal);
int regra_11(std::vector<Vertex> &vertices, int &qnd_vt, int &qnd_t, int &qnd_w, int **matriz_BL, int grafoOriginal);
int regra_12(std::vector<Vertex> &vertices, int &qnd_vt, int &qnd_t, int &qnd_w, int grafoOriginal);
int regra_13(std::vector<Vertex> &vertices, int &qnd_w);
int regra_14(std::vector<Vertex> &vertices, int &qnd_vt, double **matriz_custo, double &prize, int &qnd_t, int &qnd_w, int grafoOriginal);

#endif /* REGRAS_H_ */
