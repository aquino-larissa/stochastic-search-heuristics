#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

int max_agent(int d, int M, int A[M], double fit[d]){
	double max = 0.;
	int i, num;
	int k = 0;
	for (i = 0; i < M; i++){
		num = A[i];
		if (fit[num] > max){
			max = fit[num];
			k =   i;
		}
	}
	return k;}
int random_change(int N, int Aj){
	int V[N], vk, dec, i, jcross;
	for (i = 0; i < N; i++){
		V[N-i-1] = Aj%2;
		Aj = (Aj - V[N-i-1])/2;
		}
	jcross = rand()%N;
	V[jcross] = 1 - V[jcross];
	dec = array_to_dec(N, V);
	return dec;
	}
int select_agent(int d, int M, int A[M], double f[d], int a){
    //sorteia um segundo agente, distinto do primeiro, com probabilidade proporcional ao fitness
	int i, num;
	double soma[M], u, v;
    soma[0] = 0.;
    if(a!= 0){
        num = A[0];
        soma[0]  = f[num];
    }
	for (i = 1; i < M; i++){
		if (i != a){
			num = A[i];
			soma[i] = soma[i-1] + f[num];
		}
	}
	u = rand();
	u = u/RAND_MAX;
    v = 0;
    for(i = 0; i < M; i++){
        if (i != a && a != (M-1)){
            v = soma[i];
            v = v/soma[M-1];
            }
        if(v > u){break;}
        }
	return i;}
int crossover(int N, int Aa, int Ab){
	int sun, Va[N], Vb[N], Vsun[N], i, jcross;
	for(i = 0; i < N; i++){
        Vsun[i] = 0;
		Va[N-i-1] = Aa%2;
		Vb[N-i-1] = Ab%2;
		Aa = Aa/2;
		Ab = Ab/2;
	}
    jcross = 1 + (rand()%(N-3));
	for (i = 0; i < N; i++){
		if(i<jcross){Vsun[i] = Va[i];}
		else{Vsun[i] = Vb[i];}
	}
	sun = array_to_dec(N, Vsun);
	return sun;}
int imitate_state(int N, int Aj, int Amax){ //Amax é o agente com máximo fitness e Aj o agente que irá copiá-lo
	int Vmax[N], V[N], i, num, elements[N], q, dec, k;
	int count = 0;

	for (i = 0; i < N; i++){
    V[N-i-1] = Aj % 2;
    Aj = (Aj - V[N-i-1])/2;
    Vmax[N-i-1] = Amax % 2;
		Amax = (Amax - Vmax[N-i-1])/2;
    if (V[N-i-1] != Vmax[N-i-1]){
      elements[count] = N-i-1;
		  count = count + 1;}
	}
	if (count > 0){
		q = rand()%count;
		num = elements[q];
		V[num] = Vmax[num];
		dec = array_to_dec(N, V);
  }
  else{
    k = rand()%N;
  	V[k] = 1 - V[k];
  	dec = array_to_dec(N, V);}
	return  dec;}
int imitative_learning(int N, int K, int M, int d, double fit[d], int max, double p){
  int l, i, new_A[M], A[M], max_a;
  int flag = 0;
  double q;
  int time_search = 1;
  for (l = 0; l < M; l++){
    A[l] = rand()%d;
    if(A[l] == max){
      flag = 1;
      break;}
  }
    
  while(flag == 0){
    time_search = time_search + 1;
    max_a = max_agent(d, M, A, fit);

    for (i = 0; i < M; i++){
      q = rand();
      q = q/RAND_MAX;
      if (q < p ){new_A[i] = random_change(N, A[i]);}
      else {
        if(A[i] == A[max_a]){new_A[i] = A[i];}
        else{new_A[i] = imitate_state(N, A[i], A[max_a]);}
        }
      if (new_A[i] == max){
        flag = 1;
        break;
        }
      }
    for(i = 0; i < M; i++){A[i] = new_A[i];}
  }
  return time_search;
  }
int sexual_genetic_algorithm(int N, int K, int M, int d, double fit[d], int max, double p){ 
    int l, m, new_A[M], A[M], max_a, parentA, parentB;
    double q;
    int flag = 0;
    int time_search = 1;
    for (l = 0; l < M; l++){
        A[l] = rand()%d;
        if(A[l] == max){
            flag = 1;
            break;}
    }
    
    while(flag == 0 ){
        time_search = time_search + 1;
        for (m = 0; m < M; m++){
            q = rand();
            q = q/RAND_MAX;
            if (q < p ){
                new_A[m] = random_change(N, A[m]);
                if (new_A[m] == max){
                    flag = 1;
                    break;
                }
            }
            else {
                parentA =  select_agent(d, M, A, fit, -1);
                parentB =  select_agent(d, M, A, fit, parentA);
                new_A[m] = crossover(N, A[parentA], A[parentB]);
                if (new_A[m] == max){
                    flag = 1;
                    break;
                }
            }
        }
        for(m = 0; m < M; m++){A[m] = new_A[m];}
        }  
    return time_search;
    }
int asexual_genetic_algorithm(int N, int K, int M, int d, double fit[d], int max, double p){
    int l, m, new_A[M], A[M], max_a, parentA, parentB;
    double q;
    int flag = 0;
    int time_search = 1;
    for (l = 0; l < M; l++){
        A[l] = rand()%d;
        if(A[l] == max){
            flag = 1;
            break;}
    }
    
    while(flag == 0 ){
        time_search = time_search + 1;
        for(m = 0; m < M; m++){
            q = rand();
            q = q/RAND_MAX;
            if (q < p ){
                new_A[m] = random_change(N, A[m]);
                if (new_A[m] == max){
                    flag = 1;
                    break;
                }
            }
            
            else {
                parentA =  select_agent(d, M, A, fit, -1);
                new_A[m] = A[parentA]; // random_change(N, A[parentA]);
                if (new_A[m] == max){
                    flag = 1;
                    break;
                }
            }
        }
        for(m = 0; m < M; m++){A[m] = new_A[m];}
        }  
    return time_search;
}

