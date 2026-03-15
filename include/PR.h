/*
 * pr.h
 *
 *  Created on: 24/05/2018
 *      Author: Glaubos
 */

#ifndef PR_H_
#define PR_H_

//#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>    // std::find
#include "Utils.h"



namespace std {

#define MAX 501
/*
 * Classe do path-relinking
 */
class PR {

private:
	int** eliteSet;
//	vector<pair <int*,double> > es;
	unsigned eliteMaxSize;
	long * hashTable;
	Vertice* v;
	int n; //solution size
	unsigned poolSize;
	Matriz* C;
	double* costPool;

//	vector<edgeSol> ordenadasFreq;

public:
	unsigned frequencias[MAX][MAX];
	PR();
	PR(Vertice*, Matriz*, unsigned, int);
	virtual ~PR();
	void extrairFreq(vector<edgeSol>& ordenadasFreq);
	void print(int *solucao);
	bool run(int*& si, int* sg, int n, int*& sm, Matriz* C, double custoSi,
			int&, double, double);
//	string longestCommonSubstring(const string& str1, const string& str2);
	int antecessor(int *s, int vertice);
	bool fullEliteSet();
	int sucessor(int *s, int vertice);
//	int PR::sucessor(int *s, int vertice);
	void solToVector(int *s, vector<int> &v);
	int tamSolucao(int* s, int*& vs);
	void remover_k_da_solucao(int *&rota, int k);
	void inserir_k_na_solucao(int *&solucao, int k, int pos);
	void copySol(int *vetor_from, int* vetor_to, int size);
	int calculateDiff(int *si, int *sg, int n);
	bool areDifferent(int* a, int* b);
	long hash(int* s);
	bool updatePool(int *, double cost);
	bool InPool(long);
	void printEliteSet();
	static bool costComparable(pair<int*, double> pair1,
			pair<int*, double> pair2) {
		return pair1.second < pair2.second;
	}

	unsigned getCurrentSize() const;
	void setCurrentSize(unsigned currentSize);
	unsigned getEliteMaxSize() const;
	void setEliteMaxSize(unsigned eliteMaxSize);
	const vector<pair<int*, double> >& getEs() const;
	void setEs(const vector<pair<int*, double> >& es);
	const vector<long>& getHashTable() const;
	void setHashTable(const vector<long>& hashTable);
	int findBiggerCostIndex(double cost);
	void replaceHash(long oldHash, long newHash);
	const Matriz*& getC() const;
	void setC(const Matriz*& c);
	double* getCostPool() const;
	void setCostPool(double* costPool);
	int getN() const;
	void setN(int n);
	unsigned getPoolSize() const;
	void setPoolSize(unsigned poolSize);
	const Vertice*& getV() const;
	void setV(const Vertice*& v);
	int** getEliteSet() const;
	void setEliteSet(int** eliteSet);


};

} /* namespace std */

#endif /* PR_H_ */

