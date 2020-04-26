// C++ for C programmers Homework 2
// Implement a Monte Carlo simulation that calculates the average shortest path in a graph
// connectivity matrices is used to represent graphs 
#include <cstdlib> // srand() and rand()
#include <time.h>  // time()
#include <iostream>
#include <float.h>
#include <iomanip>
using namespace std;

// number of vertex
const int vertex = 50;

// function to calculate probability
double prob() {
	return (double)rand() / (double)RAND_MAX;
}

// class to create graph using connectivity matrices representation
class Graph{
	int V;
public:
	double ** graph;
	// array constructor
	Graph(int vertices){
        V = vertices;
        graph = new double*[V];
		for (int i = 0; i < V; ++i) {
			graph[i] = new double[V];  
		}
    }
    // assign edges to the graph
	// number of edges determined by the density 
	void initialize(double density){
		for (int i = 0; i < V; ++i) {
			for (int j = i; j < V; ++j) {
			// if the index is same, do not create edge
				if (i == j) 
					graph[i][j] = false;
				else{
				// if the probability lower than density, do not create edge
				// add random distance to the assigned edges
					double rWeight = (1.0 + 1) + (((double) rand()) / (double) RAND_MAX) * (10.0 - (1.0 + 1));
					graph[i][j] = graph[j][i] = (prob() < density)*rWeight;
				}
			}
		}
	}
	// print graph 
	void printGraph(){
	cout << "This graph has " << V << " vertices" << endl;
	for (int i = 0; i < V; i++) 
       for (int j = 0; j < V; j++)       
          cout << graph[i][j] << " \n"[j == V-1]; 
}
};

// function to find minimal distance
double minDistance(double dist[], bool Cset[]){
	double min = DBL_MAX;
	int index;
	// scan the list of distance and closed set
	for (int v = 0; v < vertex; v++){
		// find the minimal number in open set
		if (dist[v] <= min && Cset[v] == false){
			min = dist[v];
			index = v;
		}
	}
	return index;
}

// function to calculate average path
double avgPathLength(double dist[]){
	double sum = 0;
	for (int d = 0; d < vertex; d++)
		sum += dist[d];
	// as we only calculate the path from node 1, we divided it by vertices-1
	double avg = sum*1.0/(vertex-1);
	return avg;
}

// function to find the shortest path
double dijkstra(double ** graph, int src){
	// create 2 array to:
	// dist = store the minimum distance from src
	// Cset = to define whether the node is on closed set or not 
	double dist[vertex];
	bool Cset[vertex];

	// initialize the array with infinity and all on open set
	for (int i = 0; i < vertex; i++){
		dist[i] = DBL_MAX;
		Cset[i] = false;
	}

	// the distance from the src always zero
	dist[src] = 0;

	// loop through other nodes
	for (int c = 0; c < vertex-1; c++){
		// set the shortest distance from src
		int u = minDistance(dist, Cset);
		// add to closed set
		Cset[u] = true;

		// loop through adjacent nodes
		for (int v = 0; v < vertex; v++){
			if (!Cset[v] && graph[u][v] && dist[u] != DBL_MAX && dist[u] + graph[u][v] < dist[v])
				dist[v] = dist[u] + graph[u][v]; 
		}
	}
	// calculate average path length
	double avg = avgPathLength(dist);
	return avg;
}

int main() {
	// set 2 significant numbers
	std::cout << std::setprecision(2) << std::fixed;
	srand(time(0));
	// set a variable g as the graph
	Graph g(vertex);
	// Monte Carlo simulation
	// calculate the average path length for a random graph N times
	// and then find the average for each density
	double dens[] = {0.2, 0.4};
	int N = 100;
	cout << "Monte Carlo simulation with " << N << " samples" << endl;
	for (int i = 0; i < 2; i++){
		double sum = 0.0;
		double density = dens[i];
		for (int j = 0; j < N; j++){
			// initialize the edges of the graph
			g.initialize(density);
			double res = dijkstra(g.graph, 0);
			sum += res;
		}
		cout << "average path length when the density is " << dens[i]*100 << "% is " <<  sum/N << endl;
	}
}