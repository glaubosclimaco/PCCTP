/*
 * pr.cpp
 *
 *  Created on: 24/05/2018
 *      Author: Glaubos
 */

#include "PR.h"
#define debug 0
namespace std {

PR::PR() {
	// TODO Auto-generated constructor stub

}

PR::~PR() {
	// TODO Auto-generated destructor stub
}

void PR::extrairFreq(vector<edgeSol>& ordenadasFreq) {
	for (unsigned i = 0; i < poolSize; i++) {
		int p = 0;
		while (eliteSet[i][p] == -1) {
			p++;
		}

		//p e primeiro vertice pertencente a solucao
		int antecessor = p;
		int sucessor = eliteSet[i][antecessor];
//		printf("\n%d  --> %d", antecessor, sucessor);
		++frequencias[antecessor][sucessor];
		int count = 1;

//		edgeSol e;
//		e.o = antecessor;
//		e.d = sucessor;
//		ordenadasFreq.frequencia = frequencias[antecessor][sucessor];
//		ordenadasFreq.push_back(e);

		while (sucessor != p) {
			antecessor = sucessor;
			sucessor = eliteSet[i][sucessor];
//			printf(" --> %d", sucessor);
			++frequencias[antecessor][sucessor];

//			edgeSol e;
//			e.o = antecessor;
//			e.d = sucessor;
//			e.frequencia = frequencias[antecessor][sucessor];
//			ordenadasFreq.push_back(e);

		}
	}

	//ordenando de forma decrescente de acordo com a frequência
//	sort(ordenadasFreq.begin(), ordenadasFreq.end(), ordenaFreq);

	printf("extraiu as frequencias do pool \n");

}


bool PR::run(int*& si, int*sg, int n, int*& sm, Matriz* C, double custoSi,
		int& indexSg, double prize, double d) {
//	printf("custo de si: %f \n", custoSi);
#if debug
	printf("--> PR:: run \n");
#endif

	int tsi, tsg, index = indexSg;
	int* vsi;
	int* vsg;
	long hashSg = 1, hashSi = 1;
	int* tmp;
	int nMoves = 0, totalMoves = 0;
	double custoSg = 0.0;
	double custoSm = custo(C, si);
	double bestCost = custoSi;

	vsi = new int[n];
	vsg = new int[n];
	tmp = new int[n];

//initializing arrays
	for (int i = 0; i < n; i++) {
		vsi[i] = -1;
		vsg[i] = -1;
	}

	bool improvedSi = false;

//choosing guide solution
//	custoSg = custo(C, si);
//	while (custo(C, si) == custoSg) {
//		index = rand() % (poolSize);
//		indexSg = index;
//
//		custoSg = custo(C, eliteSet[index]);
//	}
////
	printf("index: %d\n", index);
	printf("\n\nsolucao guia:");
	print(sg);
//	printf("  hash: ");
//	cout << hashSg << endl;
//
//	tsi = tamSolucao(si, vsi);
//	tsg = tamSolucao(eliteSet[index], vsg);

	totalMoves = calculateDiff(si, sg, n);

//verificando o que deve ser modificado em si
	for (int i = 0; i < n; ++i) {
		if (nMoves >= totalMoves) {
			break;
		}
//		printf("i: %d\n", i);
		if (si[i] != -1 and sg[i] == -1) {
//			printf("aqui 1 \n");
			//drop
//			int pos = eliteSet[index][i];

			copy_int_vetor(si, tmp, n);
			remover_k_da_solucao(si, i);
			int viavel = checarSolucao(v, C, si, prize, d, n);
			if (!viavel) {

				printf("sol nao viavel apos drop. Returning false...\n");
				copy_int_vetor(tmp, si, n);
//				return false;
				printf("sm voltou para: \n");
				print(si);
			}
//				inserir_k_na_solucao(si, i, pos);
//				break;
//			} else {
			custoSm = custo(C, si);
			if (custoSm < bestCost) {
#if debug
				printf("\n\ncusto caiu de %f para %f \n", bestCost, custoSm);
				printf("drop %d\n", i);
#endif
				if (hash(si) == hashSg) {
#if debug
					printf(
							"sm and sg are exactly the same :(. Not copying...\n");
#endif

				} else {

					improvedSi = true;
#if debug
					cout << "hash si: " << hash(si) << "\nhash sm: " << hashSg
					<< endl;
					printf("sm and sg are different :) \n");
#endif
					bestCost = custoSm;
					copySol(si, sm, n);
					++nMoves;
					//					print(sm);
//				printf("\ncusto Sm: %f \n", custoSm);

//				this->eliteSet.push_back(middle);
				}
			}
		} else {
			if (si[i] == -1 and sg[i] != -1) {
				//add
//				int k = eliteSet[index][i];
				printf("nao tem %d na sol inicial\n", i);
				int pos = antecessor(sg, i);
//				print(eliteSet[index]);
//				printf("k: %d \n", k);
//				printf("vsg[i]: %d \n", vsg[i]);
//				printf("i: %d \n", i);
//				printf("index: %d \n", index);
//				printf("pool size: %d \n", poolSize);
//				getchar();
				inserir_k_na_solucao(si, i, pos);
#if debug
				printf("add %d com antecesor %d\n", i, pos);
#endif
				custoSm = custo(C, si);
				printf("si: \n");
				print(si);
//				getchar();
				if (custoSm < bestCost) {
#if debug
					printf("\n\ncusto caiu de %f para %f \n", bestCost,
							custoSm);
//					printf("add %d\n", i);
#endif
					if (hash(si) == hashSg) {
#if debug
						printf(
								"sm and sg are exactly the same :(. Not copying...\n");
#endif
					} else {
						improvedSi = true;
#if debug
						cout << "hash si: " << hash(si) << "\nhash sm: "
						<< hashSg << endl;
						printf("sm and sg are different :) \n");
#endif
						bestCost = custoSm;
						copySol(si, sm, n);
						++nMoves;
//						return true;
						//					print(sm);
						//				printf("\ncusto Sm: %f \n", custoSm);

						//				this->eliteSet.push_back(middle);
					}
				}
			}
		}

	}

	printf("\n\nentering swap... \n");
	printf("bestCost: %f \n", bestCost);
	for (int i = 0; i < n; ++i) {
		if (nMoves >= totalMoves) {
			break;
		}
//swap
//		printf("loop %d \n", i);
		if ((si[i] != -1 and sg[i] != -1) and (si[i] != sg[i])) {

			//vertices para fazer swap
			int a = si[i];
			int b = sg[i];

			int antecessor_a = antecessor(si, a);
			int sucessor_a = si[a];

			int antecessor_b = antecessor(sg, b);
			int sucessor_b = sg[b];

			remover_k_da_solucao(si, a);
			inserir_k_na_solucao(si, a, antecessor_b);

//			//Vertices adjancentes
//			if (sucessor_a == b) { //b é sucessor de a
//				si[antecessor_a] = b;
//				si[b] = a;
//				si[a] = sucessor_b;
//			} else if (sucessor_b == a) { //a é sucessor de b
//				si[antecessor_b] = a;
//				si[a] = b;
//				si[b] = sucessor_a;
//			} else { //Vertices nao ajdacentes
//				//Trocar Posicao dos vertices
//				si[antecessor_a] = b;
//				si[b] = sucessor_a;
//				si[antecessor_b] = a;
//				si[a] = sucessor_b;
//			}

			custoSm = custo(C, si);

//			int ssi = si[i]; //sucessor de i em si
//			int ssg = eliteSet[index][i]; //sucessor de i em sg

//			printf("ssi:%d e ssg:%d \n", ssi, ssg);

//			if (ssi != ssg) {
//				remover_k_da_solucao(si, i);
//				inserir_k_na_solucao(si, i, ssg);
//				custoSm = custo(C, si);
#if debug
			printf("custoSm: %f \t bestCost: %f \n", custoSm, bestCost);
			printf("swap %d-%d\n", a, b);
			print(si);
//			getchar();
#endif
			if (custoSm < bestCost) {
#if debug
				printf("\n\ncusto caiu de %f para %f \n", bestCost, custoSm);
//					printf("drop %d\n", i);
#endif
				if (hash(si) == hashSg) {
#if debug
					printf(
							"sm and sg are exactly the same :(. Not copying...\n");
#endif
				} else {
					improvedSi = true;
#if debug
					cout << "hash si: " << hash(si) << "\nhash sm: " << hashSg
					<< endl;
					printf("sm and sg are different :) \n");
#endif
					bestCost = custoSm;
					copySol(si, sm, n);
//						return true;
					++nMoves;
					//					print(sm);
					//				printf("\ncusto Sm: %f \n", custoSm);

					//				this->eliteSet.push_back(middle);
				}
			}
//			}
		}
	}

//	printf("fim loop swap \n");

//	delete vsi;
//	delete vsg;

//	printf("-- fim run \n");

//	printf("sg:");
//	print(eliteSet[index].first);
//	printf("custo de sg: %f \n", custo(C, eliteSet[index].first));

	if (!improvedSi)
		return false;

//	//imprimindo middle
//
//	printf("----------- imprimindo Sm ----------- \n");
//	int p = 0;
//	int cont = 0;
//	while (sm[p] == -1) {
//		p++;
//	}
//
//	//p e primeiro vertice pertencente a solucao
//	int antecessor = p;
//	int sucessor = sm[antecessor];
//	cont = 1;
//	printf("%d -> %d ", antecessor, sucessor);
//	while (sucessor != p) {
//		antecessor = sucessor;
//		sucessor = sm[sucessor];
//		++cont;
//		printf("-> %d ", sucessor);
//	}
//
//	printf("\n");
//	getchar();

	return true;
}

void PR::remover_k_da_solucao(int *&rota, int k) {
	int ant = antecessor(rota, k);
	rota[ant] = rota[k];
	rota[k] = -1;
}

void PR::inserir_k_na_solucao(int *&solucao, int k, int pos) {

	int sucessor_pos = solucao[pos];

	solucao[pos] = k;
	solucao[k] = sucessor_pos;
}

int PR::tamSolucao(int *s, int*& vs) {
	int p = 0;
	int cont = 0;
	while (s[p] == -1) {
		p++;
	}

//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = s[antecessor];
	cont = 1;
	vs[antecessor] = 1;
	vs[sucessor] = 1;
	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = s[sucessor];
		++cont;
		vs[sucessor] = 1;
	}

//	for (int i = 0; i < 48; i++) {
//		printf("v[%d]=%d \t", i, vs[i]);
//	}

	return cont;
}

void PR::print(int *solucao) {
	int p = 0;
	while (solucao[p] == -1) {
		p++;
	}

//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = solucao[antecessor];
	printf("\n%d  --> %d", antecessor, sucessor);
	int count = 1;
	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = solucao[sucessor];
		printf(" --> %d", sucessor);

		count++;
//		if (count > 100)
//			scanf("%d", &count);

	}

//	printf("\n%d vertices\n", count);
	printf("\ncost: %f \n", custo(C, solucao));
}

//string PR::longestCommonSubstring(const string& str1, const string& str2) {
//	if (str1.empty() || str2.empty()) {
//		return 0;
//	}
//	int *curr = new int[str2.size()];
//	int *prev = new int[str2.size()];
//	int *swap = NULL;
//	int maxSubstr = 0;
//	string longest;
//	for (unsigned int i = 0; i < str1.size(); ++i) {
//		for (unsigned int j = 0; j < str2.size(); ++j) {
//			if (str1[i] != str2[j]) {
//				curr[j] = 0;
//			} else {
//				if (i == 0 || j == 0) {
//					curr[j] = 1;
//				} else {
//					curr[j] = 1 + prev[j - 1];
//				}
//				if (maxSubstr < curr[j]) {
//					maxSubstr = curr[j];
//					longest.clear();
//				}
//				if (maxSubstr == curr[j]) {
//					longest += str1.substr(i - maxSubstr + 1, i + 1);
//				}
//			}
//		}
//		swap = curr;
//		curr = prev;
//		prev = swap;
//	}
//	delete[] curr;
//	delete[] prev;
//	return longest.substr(0, maxSubstr);
//}

int PR::antecessor(int *s, int vertice) {
	int i = 0;
	int achou = 0;

	while (achou == 0) {
		if (s[i] == vertice) {
			achou = 1;
		} else {
			i++;
		}
	}

	return i;
}

int PR::sucessor(int *s, int vertice) {
	return s[vertice];
}

void PR::copySol(int *vetor_from, int* vetor_to, int size) {
//	printf("copying ... \n");
	for (int i = 0; i < size; ++i) {
		vetor_to[i] = vetor_from[i];
//		printf("%d ", vetor_to[i]);
	}
//	getchar();
}

long PR::hash(int* s) {
	long hash;
	int p = 0;
	while (s[p] == -1) {
		p++;
	}

//p e primeiro vertice pertencente a solucao
	int antecessor = p;
	int sucessor = s[antecessor];

	hash += antecessor + v[antecessor].premio;
	hash += sucessor + v[sucessor].premio;

	while (sucessor != p) {
		antecessor = sucessor;
		sucessor = s[sucessor];
		hash += sucessor + v[sucessor].premio;
	}

	return hash;

}

/*
 * try to update the elite set
 */

bool PR::updatePool(int* s, double cost) {
#if debug
	printf("-> PR::updatePool\n");
#endif
	long hashOfs;
	int indexPool = -1;
	hashOfs = hash(s);

	if (poolSize == 0) { //if eliteset is empty
#if debug
			printf("elite set empty \n");
			cout << "hash of S: " << hashOfs << endl;
#endif
		copy_int_vetor(s, eliteSet[poolSize], n);
		costPool[poolSize] = cost;
		hashTable[poolSize] = hashOfs; //add hash id to hash table
		++poolSize;

		return true;
	} else {
		bool inPool = InPool(hashOfs);
		if (inPool) {
#if debug
			cout << " S ja esta no pool \n";
			cout << "hash of S: " << hashOfs << endl;
#endif
			return false;
		} else {
			if (poolSize < eliteMaxSize) {
#if debug
				cout << "eliteSet.size() < eliteMaxSize || hash of S: "
				<< hashOfs << endl;
				printf("%d < %d", poolSize, eliteMaxSize);
#endif
				copy_int_vetor(s, eliteSet[poolSize], n);
				costPool[poolSize] = cost;
				hashTable[poolSize] = hashOfs;
				++poolSize;
				return true;
			} else {

				indexPool = findBiggerCostIndex(cost);
#if debug
				printf("index pool: %d \n", indexPool);
#endif
				if (indexPool > -1) { //there is solution inside the pool which is worse the one to be inserted
#if debug
						printf("uma sol de pior custo vai sair do pool\n");
#endif
					copy_int_vetor(s, eliteSet[indexPool], n);
					costPool[indexPool] = cost;
					hashTable[indexPool] = hashOfs;
					indexPool = -1;
				}
			}

		}

	}
	return false;
}
PR::PR(Vertice* _v, Matriz* _C, unsigned maxSizeES, int _n) :
		v(_v), C(_C) {

	poolSize = 0;
	this->eliteMaxSize = maxSizeES;
	this->n = _n;

	this->eliteSet = new int*[maxSizeES];
	for (unsigned i = 0; i < maxSizeES; ++i) {
		eliteSet[i] = new int[n];
	}

	hashTable = new long[maxSizeES];
	costPool = new double[n];

}

bool PR::InPool(long hash) {
	for (unsigned i = 0; i < poolSize; ++i) {
		if (hashTable[i] == hash)
			return true;
	}
	return false;
}

double* PR::getCostPool() const {
	return costPool;
}

void PR::setCostPool(double* costPool) {
	this->costPool = costPool;
}

int PR::getN() const {
	return n;
}

void PR::setN(int n) {
	this->n = n;
}

unsigned PR::getPoolSize() const {
	return poolSize;
}

void PR::setPoolSize(unsigned poolSize) {
	this->poolSize = poolSize;
}

int** PR::getEliteSet() const {
	return eliteSet;
}

int PR::calculateDiff(int* si, int* sg, int n) {
	int p = 0;
	int nMoves = 0;

	for (int i = 0; i < n; ++i) {
		if (si[i] != sg[i]) {
			++nMoves;
		}
	}

	return nMoves;
}



void PR::setEliteSet(int** eliteSet) {
	this->eliteSet = eliteSet;
}

void PR::replaceHash(long oldHash, long newHash) {
	vector<long>::iterator it;
	for (unsigned i = 0; i < poolSize; ++i) {
		if (hashTable[i] == oldHash) {
			hashTable[i] = newHash;
		} else {
//			cout << "\n" << oldHash << "isnt in the pool\n";
		}
	}
}

unsigned PR::getCurrentSize() const {
	return poolSize;
}

void PR::setCurrentSize(unsigned currentSize) {
	this->poolSize = currentSize;
}

unsigned PR::getEliteMaxSize() const {
	return eliteMaxSize;
}

void PR::setEliteMaxSize(unsigned eliteMaxSize) {
	this->eliteMaxSize = eliteMaxSize;
}

void PR::printEliteSet() {
	printf("\n----------------- ELITE SET -----------------------\n");
	printf("size: %d \n", poolSize);
	for (unsigned i = 0; i < poolSize; i++) {
		cout << "\nsolucao " << i << ":";
//		print(es[i].first);
		print(eliteSet[i]);
//		cout << "hash: " << hash(L[i].first) << endl;
		cout << "hash: " << hash(eliteSet[i]) << endl;
		cout << "custo: " << custo(C, eliteSet[i]) << endl;
	}
	printf("\n----------------- ELITE SET -----------------------\n");
}

int PR::findBiggerCostIndex(double cost) {
	unsigned i;
	int ind = -1;
	double biggerCost;
	biggerCost = cost;

	for (i = 0; i < poolSize; ++i) {
		if (biggerCost <= costPool[i]) {
			biggerCost = costPool[i];
			ind = i;
		}
	}

	return ind;
}

} /* namespace std */
