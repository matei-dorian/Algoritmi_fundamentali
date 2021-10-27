#include <iostream>
#include <vector>
#include <queue>
#include <fstream>

using namespace std;

ifstream fin("dfs.in");
ofstream fout("dfs.out");

class Graph{

protected:
    int n, m; /// n -> number of nodes / m -> number of edges
    vector<vector<int> > l; /// list of adjacency

public:
    Graph(int n = 1, int m = 0);
    Graph(const Graph &g);

    Graph& operator= (const Graph &g);


    int get_numEdges();
    int get_numNodes();

    friend istream& operator>>(istream& in, Graph& obj);
    friend ostream& operator<<(ostream& out, const Graph& obj);

    virtual ostream& print(ostream& out) const;
    virtual istream& read(istream& in);

    void dfs(int start);
    void do_dfs(int start, vector<bool> &viz, bool print);
    int num_of_components();

    void bfs(int start);
    void do_bfs(int start, vector<int> &dist, bool print);
    void dist_from_node(int start);

    virtual ~Graph() = 0;
};


///constructors
Graph :: Graph (int n, int m) : l(n + 1){
    /// the default case is to construct a graph
    /// with one node and 0 edges

    this->n = n;
    this->m = m;

}

Graph :: Graph (const Graph &g) : l(g.n) {
    /// copy constructor

    this->n = g.n;
    this->m = g.m;

    for(int i = 0; i < g.n; ++i)
        l.push_back(g.l[i]);
}


/// getters
int Graph :: get_numEdges() {
    return m;
}

int Graph :: get_numNodes() {
    return n;
}


/// operators
Graph& Graph :: operator= (const Graph &g){
     if(this != &g){
        this->n = g.n;
        this->m = g.m;

        if(!this->l.empty())
            l.clear();

        for(int i = 0; i <= g.n; ++i)
            l.push_back(g.l[i]);
     }

     return *this;
}

istream& operator>>(istream& in, Graph& g){
    g.read(in);
    return in;
}

ostream& operator<<(ostream& out, Graph& g){
    g.print(out);
    return out;
}


/// destructor
Graph :: ~Graph(){
    if(!this->l.empty())
        l.clear();
}


/// methods
istream& Graph :: read(istream& in){
    /// the default format is n m (same line) followed
    /// by m pairs x y, each on one line, representing
    /// the edges

    cin >> n >> m;
    return in;
}

ostream& Graph :: print(ostream& out) const{
    /// the default format is n m on one line followed
    /// by the number of the current node and a list of every
    /// reachable nodes

    cout << n << " " << m << "\n";

    for(int i = 1; i <= n; ++i){
        cout << i << ": ";
            for(unsigned int j = 0; j < l[i].size() ; ++j)
                cout << l[i][j] << " ";
        cout << "\n";
    }

    return out;
}

void Graph :: do_dfs(int start, vector<bool> &viz, bool print){
    /// do DFS from start node
    /// if print is true -> print the node as
    /// you visit them

    if(print)
        fout << start << " ";

    viz[start] = 1;

    for(auto it = l[start].begin(); it != l[start].end(); ++it)
        if(!viz[*it])
            do_dfs(*it, viz, print);

}

void Graph::dfs(int start){
    /// this function initializes the vector viz
    /// used to mark visited node and calls DFS
    /// this function prints the DFS traversal
    /// of the graph with the start node - start

    vector<bool> viz(n + 1, 0);
    do_dfs(start, viz, 1);
}

int Graph :: num_of_components(){
    /// this function returns the number
    /// of connected components of the graph

    int ct = 0; /// variable to count the components
    vector<bool> viz(n + 1, 0);

    for(int i = 1; i <= n; ++i){
        if(!viz[i]){
            ct++;
            do_dfs(i, viz, 0);
        }
    }

    return ct;
}

void Graph :: do_bfs(int start, vector<int>& dist, bool print = 0){
    /// do BFS from start node
    /// mark in dist[i] the minimum distance
    /// between i and start

    queue<int> q;
    vector<bool> viz(n + 1, 0);

    /// initialization
    viz[start] = 1;
    q.push(start);
    dist[start] = 1;

    /// BFS
    while(!q.empty()){
        int k = q.front();

        if(print)
            cout<< k << " ";

        for(auto i: l[k])
            if(viz[i] == 0){
                viz[i] = 1;
                dist[i] = dist[k] + 1;
                q.push(i);
           }
        q.pop();
    }
}

void Graph :: bfs(int start){
    /// this function initializes the vector dist
    /// and calls BFS
    /// this function prints the BFS traversal
    /// of the graph with the start node - start

    vector<int> dist(n + 1, 0);
    do_bfs(start, dist, 1);
}

void Graph :: dist_from_node (int start){
    /// this function calculates the distance between
    /// the start node and every other node in the graph
    /// if the node is unreachable -> the distance is -1

    vector<int>dist(n + 1, 0);
    do_bfs(start, dist);

    for(int i = 1; i <= n; i ++)
        fout << dist[i] - 1<<" ";
}

//----------------------------------------------------------------------

class UnorientedGraph : public Graph{

public:
    UnorientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    UnorientedGraph(int n, int m, vector<pair<int, int> > v);
    UnorientedGraph& operator+=(const pair<int, int> edge);
    istream& read(istream& in) override;

};


/// Constructor
UnorientedGraph :: UnorientedGraph(int n, int m, vector<pair<int, int> > v) : Graph(n, m){
    for(auto it : v){
        l[it.first].push_back(it.second);
        l[it.second].push_back(it.first);
    }
}


/// read method
istream& UnorientedGraph::read(istream& in){
    Graph::read(in); /// call read from base
    UnorientedGraph temp(n, m);

    for(int i = 1; i <= m; i ++){
        int x, y;
        cin >> x >> y;
        temp.l[x].push_back(y);
        temp.l[y].push_back(x);
    }

    operator=(temp);
    return in;
}


/// Operators
UnorientedGraph& UnorientedGraph :: operator+=(const pair<int, int> edge){
    /// add edge to graph

    l[edge.first].push_back(edge.second);
    l[edge.second].push_back(edge.first);
    return *this;
}

//----------------------------------------------------------------------

class OrientedGraph : public Graph{

public:
    OrientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    OrientedGraph(int n, int m, vector<pair<int, int> > v);
    OrientedGraph& operator+=(const pair<int, int> edge);
    istream& read(istream& in) override;

};

/// Constructor
OrientedGraph :: OrientedGraph(int n, int m, vector<pair<int, int> > v) : Graph(n, m){
    for(auto it : v)
        l[it.first].push_back(it.second);
}


/// read method
istream& OrientedGraph :: read(istream& in){
    Graph::read(in); /// call read from base
    OrientedGraph temp(n, m);

    for(int i = 1; i <= m; i ++){
        int x, y;
        cin >> x >> y;
        temp.l[x].push_back(y);
    }

    operator=(temp);
    return in;
}


/// Operators
OrientedGraph& OrientedGraph :: operator+=(const pair<int, int> edge){
    /// add edge to graph

    l[edge.first].push_back(edge.second);
    return *this;
}

//----------------------------------------------------------------------


/*void Solve(){
    int n, m, start;

    fin >> n >> m >> start;

    OrientedGraph g(n, m);

    for(int i = 1; i <= m; i ++){
        int x, y;
        fin >> x >> y;
        g += make_pair(x, y);
    }

    g.dist_from_node(start);
} BFS infoarena 100 pcte */

void Solve(){
    int n, m;
    fin >> n >> m;

    UnorientedGraph g(n, m);
    for(int i = 1; i <= m; i ++){
        int x, y;
        fin >> x >> y;
        g += make_pair(x, y);
    }

    fout << g.num_of_components();
}


int main()
{
    Solve();
    return 0;
}

