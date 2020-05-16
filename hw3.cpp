// C++ for C programmers Homework 3
// Compute the minimum spanning tree for an inputted graph
// Kruskal's algorithm is used to find MST
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

// create a pair to store edges
typedef pair<int, int> edge;

// class to create definition of graph
class Graph{
    int V;
public:
    // use vector of pairs for adjacency list
    vector<pair<int,edge>> adjList;
    // graph constructor
    Graph(int vertices){
        V = vertices;
    }

    // edge constructor for undirected graph
    void addEdge(int src, int target, int weight){
        adjList.push_back(make_pair(weight,edge(src,target)));
        adjList.push_back(make_pair(weight,edge(target,src)));
    }
};

// function to find the connection between nodes
int find(int parent[], int i)  
{  
    // if there is a connection, return node
    if (parent[i] == i)  
        return i;  
    // recursion to find parent nodes
    return find(parent, parent[i]);  
}  

// function to connect 2 nodes
void Union(int parent[],int u, int v) {
    parent[u] = parent[v];
}

// implementation of Kruskal's algorithm
void kruskal(Graph &graph) {
    int i, uRep, vRep;

    // size of the graph edges
    int S = graph.adjList.size();
    // sort the edges according to weights
    std::sort(graph.adjList.begin(), graph.adjList.end()); 

    // Define vector Tree to store MST
    // define parent array to store representative nodes
    vector<pair<int,edge>> Tree;
    int* parent;
    parent = new int[S];

    // initialize parent array, each subset contains only single element
    for (int v = 0; v < S; v++)
        parent[v] = v;
    
    // loop through all edges in graph
    for (i = 0; i < S; i++) {
        // use union-find algorithm to detect cycles
        // find connection between two nodes in the edge
        uRep = find(parent, graph.adjList[i].second.first);
        vRep = find(parent, graph.adjList[i].second.second);
        // if not a cycle, add to tree
        if (uRep != vRep) {
            Tree.push_back(graph.adjList[i]); 
            // create union between nodes, update the parent array
            Union(parent, uRep, vRep);
        }
    }
    // print the tree and calculate the MST
    int sum = 0;
    cout << "Edges forming the MST: " << endl;
    for (int j = 0; j < Tree.size(); j++){
        cout << "From " << Tree[j].second.first << " to " << Tree[j].second.second << 
        " with length of " << Tree[j].first << endl;
        sum += Tree[j].first;
    }
    // MST
    cout << "MST length: " << sum << endl;
}


int main() {
    // Read the input file and set the number of vertices
    string sLine;
    ifstream infile;
    infile.open("SampleTestData_mst_data.txt");
    // Read the first line for the number of vertices
    int N;
    getline(infile, sLine);
    N = stoi(sLine);
    // Create the graph
    Graph g(N);
    // Read the rest of the line to set the edges
    int i, j, cost;
    while(infile >> i >> j >> cost){
        g.addEdge(i, j, cost);
    }
    // Kruskal algorithm
    kruskal(g);
	return 0;
}