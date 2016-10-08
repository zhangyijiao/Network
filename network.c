#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mt19937ar.c"

struct Link
{
	int num;
	struct Link *next;
};

void AddNode(struct Link *net, int i) 
{
	struct Link *p;
	p = (struct Link *)malloc(sizeof(struct Link));
	p->next = net->next;
	p->num = i;
	net->next = p;
}

void AddLink(struct Link *net, int i, int j)
{
	AddNode(&net[i], j);
	AddNode(&net[j], i);
}

void DeleteNode(struct Link *net, int i, int j) 
{
    struct Link *p, *pp;
    p = &net[i];
    while(p->next != NULL) {
        if(p->next->num == j) {
            pp = p->next;
            p->next = pp->next;
            free(pp);
            break;
        }
        p = p->next;
    }
}

void DeleteLink(struct Link *net, int i, int j) 
{
    DeleteNode(net, i, j);
    DeleteNode(net, j, i);
}

void ER_Network(struct Link *net, int n, double probability)
{
	int i, j;
	double ran;
	for(i = 0; i < n; i++) {
		for(j = i + 1; j < n; j++) {
			ran = genrand_real3();  //generates a random number on (0,1)-real-interval
			if(ran < probability) {
				AddLink(net, i, j);
			}
		}
	}
}

//each step add a new node with m links
void BA_Network(struct Link *net, int n, int m)
{
	int i, j, k, total_degree = 0;
	double temp;
	int *degree = (int *)malloc(n*sizeof(int));
	struct Link *p;
	for(i = 0 ; i < n; i++) {
		degree[i] = 0;
	}
	for(i = 0; i < m+1; i++) {
		for(j = i+1; j < m+1; j++) {
			AddLink(net, i, j);
		}
		degree[i] = m;
		total_degree += m;
	}
	for(i = m+1; i < n; i++) {
		j = 0;
		while(j < m) {
			temp = total_degree*genrand_real3();
			k = -1;
			while(temp > 0) {
				k++;
				temp = temp - degree[k];
			}
			p = &net[k];
			while(p->next != NULL) {
				p = p->next;
				if(p->num == i) break; //the chosen node can't be i's neighbor already
			}
			if(p->next == NULL) {
				AddLink(net, i, k);
				degree[i] += 1;
				degree[k] += 1;
				total_degree += 1;
				j++;
 			}
		}
		total_degree += m;
	}
	free(degree);
}

//if j == i or j is i's neighbor, return 1
int YesOrNot(struct Link *net, int i, int j)
{
    struct Link *p;
    int index = 0;
    p = &net[i];
    if(i == j) {
        index = 1;
    }
    while(p->next != NULL) {
        p = p->next;
        if(p->num == j){
            index = 1;
            break;
        }
    }
    return index;
}

//break the link between i and j, then choose a new node connect to i
void Reconnect(struct Link *net, int n, int i, int j) 
{
	int ran;
	DeleteLink(net, i, j);
	ran = (int)(n*genrand_real2()); //[0,1)
	while(YesOrNot(net, i, ran)) {
		ran = (int)(n*genrand_real2());
	}
	AddLink(net, i, ran);
}

//at the beginning, a node connects to 2r neighbors
//average degree is 2r
//prob: reconnect probability 
void SW_network(struct Link *net, int n, int r, double prob) 
{
	int i, j, neighbor;
	struct Link *p;
	for(i = 0; i < n; i++) {
		for(j = 0; j < r; j++) {
			AddLink(net, i, (i + j + 1)%n);
		}
	}
	for(i = 0; i < n; i++) {
		for(j = 0; j < r; j++){
			neighbor = (i+j+1)%n;
			if(genrand_real3() < prob) {
				Reconnect(net, n, i, neighbor);
			}
		}
	}
}

void Lattice1D(struct Link *net, int n)
{
	int i;
	for(i = 0; i < n; i++) {
		AddLink(net, i ,(i+1)%n);
	}
}

//take the square root of n equals l
//periodic boundary condition
void Lattice2D(struct Link *net, int n) 
{
	int i, j, l;
	double size = (double)n;
	l = sqrt(size);
	for(i = 0; i < l; i++) {
		for(j = 0; j < l; j++) {
			AddLink(net, i*l + j, i*l + (j + 1)%l);
			AddLink(net, i*l + j, ((i + 1)%l)*l + j);
		}
	}
}

float AverageDegree(struct Link *net, int n)
{
	int i, temp, total_degree = 0;
	float average_degree;
	struct Link *p;
	for(i = 0; i < n; i++) {
		p = &net[i];
		temp = 0;
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		total_degree += temp;
	}
	average_degree = total_degree*1.0/n;
	return average_degree;
}

void FreeNetwork(struct Link *net, int n)
{
	int i;
	struct Link *p, *pp;
	for(i = 0; i < n; i++) {
		p = &net[i];
		while(p->next != NULL) {
			pp = p->next;
			p->next = pp->next;
			free(pp);
		}
	}
}

int main()
{
	int i, n = 10000, m = 3, r = 2; 
	double probability = 4./(n - 1); // 4 is average degree of a random network, changeable
	double prob = 0.2;
	struct Link *net = (struct Link*)malloc(n*sizeof(struct Link));
	struct Link *p;
	float average_degree;
	//initialize random number generator
	unsigned long idum = (unsigned long)time(NULL);
	init_genrand(idum); 
	//initialize network
	for(i = 0; i < n; i++) {
		net[i].num = i;
		net[i].next = NULL;
	}
	Lattice2D(net, n);

	return 0;
}