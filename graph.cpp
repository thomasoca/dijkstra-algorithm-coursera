#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cstdlib> // srand() and rand()
#include <time.h>  // time()
using namespace std;

// create a pair to store node and weight
typedef pair<int, double> Pair;

// prob function to assign edges
double prob() {
	return (double)rand() / (double)RAND_MAX;
}

// class to create definition of graph
class Graph{
    int V;
public:
    // use vector of vectors for adjacency list
    vector<vector<Pair>> adjList;

    // graph constructor
    Graph(int vertices){
        V = vertices;
        adjList.resize(V);
    }

    // edge constructor for undirected graph
    void addEdge(int src, int target, double weight){
        adjList[src].push_back(make_pair(target,weight));
        adjList[target].push_back(make_pair(src,weight));
    }
};

void printGraph(Graph const &graph, int N)
{
	for (int i = 0; i < N; i++)
	{
		// print all neighboring vertices of given vertex
		for (Pair v : graph.adjList[i])
			cout << "(from " << i << " to " << v.first <<
					"; length: " << v.second << ") ";
		cout << endl;
	}
}


int main()
{
    std::cout << std::setprecision(2) << std::fixed;
    // create graph with N nodes
    int N = 5;
    Graph g(N);
    int edges = 0;
    // create edges in the graph (src, target, weight)
    // use monte carlo simulation to create edges and weights
    // number of edges determined by the density 
    srand(time(0));
    double density = 0.5;
    for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
            // if the index is same or the probability lower than density
            // do not create edge
			if (i == j || (prob() < density) == false) 
                continue;
            else{
                // assign random weight between 1.0 and 10.0
                double rWeight = (1.0 + 1) + (((double) rand()) / (double) RAND_MAX) * (10.0 - (1.0 + 1));
                g.addEdge(i, j, rWeight);
                edges += 1;
            }    
        }
    }
    // print the graph
    printGraph(g, N);
    cout << "edges: " << edges << " | " << "nodes: " << N << endl;
    return 0;
}