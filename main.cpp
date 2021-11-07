#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>
using namespace std;

ofstream fout("file.out");

class Graph{

protected:
    int n, m; /// n -> number of nodes / m -> number of edges
    vector<vector<int> > l; /// list of adjacency

    void do_dfs(const int start, vector<bool> &viz, const bool print) const;
    void do_bfs(const int start, vector<int> &dist, const bool print) const;

public:
    Graph(int n = 1, int m = 0);
    Graph(const Graph &);
    Graph& operator= (const Graph &);

    int get_numEdges() const;
    int get_numNodes() const;

    friend istream& operator>>(istream&, Graph&);
    friend ostream& operator<<(ostream&, const Graph&);

    virtual ostream& print(ostream&) const;
    virtual istream& read(istream&);
    virtual void read_from_file(string) = 0;

    void dfs(const int) const;
    int num_of_components() const;
    void bfs(const int) const;
    void dist_from_node(const int) const;

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
int Graph :: get_numEdges() const {
    return m;
}

int Graph :: get_numNodes() const {
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

ostream& Graph :: print(ostream& out) const {
    /// the default format is n m on one line followed
    /// by the number of the current node and a list of every
    /// reachable nodes

    cout << n << " " << m << "\n";

    for(int i = 1; i <= n; ++i){
        cout << i << ": ";
            for(auto node : l[i])
                cout << node << " ";
        cout << "\n";
    }

    return out;
}

void Graph :: do_dfs(const int start, vector<bool> &viz, const bool print) const {
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

void Graph :: dfs(const int start) const {
    /// this function initializes the vector viz
    /// used to mark visited node and calls DFS
    /// this function prints the DFS traversal
    /// of the graph with the start node - start

    vector<bool> viz(n + 1, 0);
    do_dfs(start, viz, 1);
}

int Graph :: num_of_components() const {
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

void Graph :: do_bfs(const int start, vector<int>& dist, const bool print = 0) const {
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

void Graph :: bfs(const int start) const {
    /// this function initializes the vector dist
    /// and calls BFS
    /// this function prints the BFS traversal
    /// of the graph with the start node - start

    vector<int> dist(n + 1, 0);
    do_bfs(start, dist, 1);
}

void Graph :: dist_from_node (const int start) const {
    /// this function calculates the distance between
    /// the start node and every other node in the graph
    /// if the node is unreachable -> the distance is -1

    vector<int>dist(n + 1, 0);
    do_bfs(start, dist);

    for(int i = 1; i <= n; i ++)
        fout << dist[i] - 1<<" ";
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


class UnorientedGraph : public Graph{

public:
    UnorientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    UnorientedGraph(int , int , vector<pair<int, int> > );
    UnorientedGraph& operator+=(const pair<int, int>);
    istream& read(istream&) override;
    void read_from_file(string);
    void print_bridges() const;
    void biconex() const;
    bool havel_hakimi(const vector<int>);

private:
    void h_merge(vector<pair<int, int> > &, const int) const;
    void do_dfs_with_bridges(const int, int &, vector<bool> &, vector<int> &, vector<int> &, vector<int> &) const;
    void do_biconex(const int, const int, vector<bool> &, vector<int> &, vector<int> &, stack<int> &, string &, int &) const;
};


/// Constructor
UnorientedGraph :: UnorientedGraph(int n, int m, vector<pair<int, int> > v) : Graph(n, m){
    for(auto it : v){
        l[it.first].push_back(it.second);
        l[it.second].push_back(it.first);
    }
}


/// read method
istream& UnorientedGraph :: read(istream& in){
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

void UnorientedGraph :: read_from_file(string filename){
    /// function to read a graph from a file filename

    int n, m;

    ifstream input;
    input.open(filename);

    input >> n >> m;

    UnorientedGraph g(n, m);
    for(int i = 1; i <= m; i ++){
        int x, y;
        input >> x >> y;
        g += make_pair(x, y);
    }
    *this = g;
    input.close();

}


/// Operators
UnorientedGraph& UnorientedGraph :: operator+=(const pair<int, int> edge){
    /// add edge to graph
    m++;
    l[edge.first].push_back(edge.second);
    l[edge.second].push_back(edge.first);
    return *this;
}

/// methods
void UnorientedGraph :: do_dfs_with_bridges(const int start, int &time, vector<bool> &viz, vector<int> &disc, vector<int> &low, vector<int> &parent) const {
    /// DFS that prints bridges
    /// used by print_bridges
    /// disc[i] = discovery time of node i
    /// parent[i] = x if x is parent of i
    /// low[i] = earliest visited node reachable from subtree rooted in i

    viz[start] = 1;
    disc[start] = low[start] = ++time;

    for(auto node : l[start]){
        if(!viz[node]){
            parent[node] = start;
            do_dfs_with_bridges(node, time, viz, disc, low, parent);

            low[start] = min(low[start], low[node]);
            if (low[node] > disc[start])
                cout << start <<" " << node << endl;
        }
        else if(node != parent[start])  /// back-edge
            low[start]  = min(low[start], disc[node]);
    }

}

void UnorientedGraph :: print_bridges() const {
    /// Method that prints all bridges in a graph

    vector<bool> viz(n + 1, 0);
    vector<int> disc(n + 1);
    vector<int> low(n + 1);
    vector<int> parent(n + 1, -1);

    int time = 0;

    for(int i = 1; i <= n; i++)
        if (!viz[i])
            do_dfs_with_bridges(i, time, viz, disc, low, parent);
}

void UnorientedGraph :: h_merge (vector<pair<int, int> > &d, const int k) const {
    /// Utility function used by havel_hakimi to keep d
    /// sorted - merges first k elements with last n - k

    int i = 0;
    unsigned int j = k + 1;
    vector<pair<int, int> > temp;

    while(i <= k && j < d.size())
        if(d[i]. first > d[j]. first){
            temp.push_back(d[i]);
            i++;
        }
        else{
            temp.push_back(d[j]);
            j++;
        }

    while(i <= k){
        temp.push_back(d[i]);
        i++;
    }

    while(j < d.size()){
        temp.push_back(d[j]);
        j++;
    }

    d = temp;
}

bool UnorientedGraph :: havel_hakimi(vector<int> v){
    ///this method receives a vector of integers
    /// as parameters and returns true if there is a graph
    /// with the given integers as degrees or false otherwise
    /// if it returns true, *this will be the graph formed

    int n_nodes = v.size(), s = 0;
    UnorientedGraph temp(n_nodes);
    vector<pair<int, int> > d(n_nodes); /// d[i].first = v[i], d[i].second = number of node

    for(auto elem : v){
        if(elem > n_nodes) return 0;
        s += elem;
    }

    if(s % 2 == 1) return 0;

    for(int i = 1; i <= n_nodes; ++i)
        d[i - 1] = make_pair(v[i - 1], i);

    sort(d.begin(), d.end(), [](pair<int, int> a, pair<int, int> b) { return a.first > b.first;});

    while(true){
        int k = d[0].first;
        if(d[0].first == 0) {*this = temp; return 1;}

        for(unsigned int i = 1; i < d.size() && d[0].first; ++i){
            d[0].first--;
            d[i].first--;
            if(d[i].first < 0) return 0;
            temp += make_pair(d[0].second, d[i].second);
        }
        d.erase(d.begin());
        h_merge(d, k);

    }

}

void UnorientedGraph :: do_biconex(const int son, const int dad, vector<bool> &viz, vector<int> &disc, vector<int> &low, stack<int> &s, string &response, int &ct) const {

    /// Utility function used by biconex
    /// viz[i] = 1 if i was visited
    /// disc[i] = discovery time of i
    /// low[i] = earliest visited node reachable from subtree rooted in i
    /// s = to store visited nodes
    /// response = string that contains all biconnected components
    /// ct = number of biconnected components

    viz[son] = 1;
    disc[son] = disc[dad] + 1;
    low[son] = disc[son];

    for(auto node : l[son])
        if(node != dad){
            if(!viz[node]){ ///back-edge
                s.push(node);
                do_biconex(node, son, viz, disc, low, s, response, ct);

                low[son] = min(low[son], low[node]);

                if(disc[son] <= low[node]){
                    ct++;
                    s.push(son);

                    while(!s.empty() && s.top() != node){
                        response += (to_string(s.top()) + " ");
                        s.pop();
                    }
                    if(!s.empty()){
                        response += (to_string(s.top()) + " ");
                        s.pop();
                    }
                    response += "\n";
                }
            }
            else
                if(disc[node] < low[son])
                    low[son] = disc[node];
        }
}

void UnorientedGraph :: biconex() const {

    /// Method that prints out the number of biconnected components
    /// and the biconnected components

    stack<int> s;
    vector<int> disc(n + 1, -1);
    vector<int> low(n + 1, -1);
    vector<bool> viz(n + 1, 0);
    string response = "";
    int ct = 0;

    disc[1] = 1;
    do_biconex(2, 1, viz, disc, low, s, response, ct);

    fout<<ct<<"\n";
    fout<<response;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


class OrientedGraph : public Graph{

public:
    OrientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    OrientedGraph(int, int, vector<pair<int, int> >);
    OrientedGraph& operator+=(const pair<int, int>);
    istream& read(istream&) override;
    void topological_sort() const;
    void read_from_file(string);
    OrientedGraph transpose_graph() const;
    void ctc() const;

private:
     inline void do_dfs(const int, vector<bool> &, stack<int> &) const;

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

void OrientedGraph :: read_from_file(string filename){
    /// Method to read a graph from a file filename

    int n, m;

    ifstream input;
    input.open(filename);

    input >> n >> m;

    OrientedGraph g(n, m);
    for(int i = 1; i <= m; i ++){
        int x, y;
        input >> x >> y;
        g += make_pair(x, y);
    }
    *this = g;
    input.close();

}


/// Operators
OrientedGraph& OrientedGraph :: operator+=(const pair<int, int> edge){
    /// add edge to graph
    m++;
    l[edge.first].push_back(edge.second);
    return *this;
}


/// methods
void OrientedGraph :: do_dfs(const int start, vector<bool> &viz, stack<int> &S) const{
    /// Method that does dfs on a graph
    /// S keeps the finishing times of all nodes

    viz[start] = 1;
    for(auto node : l[start])
        if(!viz[node])
            do_dfs(node, viz, S);

    S.push(start);
}

void OrientedGraph :: topological_sort() const {

    /// Method that prints out the nodes in topological order

    vector<bool> viz(n + 1, 0);
    stack<int> S;

    for(int i = 1; i <= n; ++ i)
        if(!viz[i])
            do_dfs(i, viz, S);

    while(!S.empty()){
        fout << S.top() << " "; /// de schimbat in cout cand pun pe github
        S.pop();
    }
}

OrientedGraph OrientedGraph :: transpose_graph() const {

    /// Utility function that returns the transpose graph of graph g
    /// Used in ctc

    OrientedGraph gt(this->n, this->m);

    for(int i = 1; i <= n; i ++)
        for(auto node : l[i])
            gt += make_pair(node, i);

    return gt;
}

void OrientedGraph :: ctc() const {

    /// method that prints out the number of strongly connected components
    /// and the strongly connected components of a graph

    OrientedGraph gt = transpose_graph();
    vector<bool> viz(n + 1, 0);
    stack<int> S;
    int ct = 0;
    int current_node;
    string response = "";

    for(int i = 1; i <= n; ++ i)
        if(!viz[i]) do_dfs(i, viz, S);

    fill(viz.begin(), viz.end(), 0); /// reinitialize viz

    while(!S.empty()){
        current_node = S.top();
        S.pop();
        if(!viz[current_node]){
            stack<int> component;

            ct ++;
            gt.do_dfs(current_node, viz, component);

            while(!component.empty()){
                response += (to_string(component.top()) + " ");
                component.pop();
            }
            response += "\n";
        }
    }

    fout << ct << "\n";
    fout << response;

}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------


/*void Solve(){
    /// number of components in a graph

    UnorientedGraph g;

    g.read_from_file("dfs.in");
    g.num_of_components();

}*/


/*void Solve() {
    /// distance to all reachable nodes

    int n, m, start;
    ifstream fin("bfs.in")

    fin >> n >> m >> start;
    OrientedGraph g(n, m);

    for(int i = 1; i <= m; i ++){
        int x, y;
        fin >> x >> y;
        g += make_pair(x, y);
    }

    fin.close();

    g.dist_from_node(start);
}*/

/*void Solve() {
    /// print nodes in topological order

    OrientedGraph g;

    g.read_from_file("topsort.in");
    g.topological_sort();
}*/

/*void Solve(){
    /// Havel Hakimi problem

    int n;
    vector<int> d;

    cin >> n;

    for(int i = 1; i <= n; i ++){
        int x;
        cin >> x;
        d.push_back(x);
    }

    UnorientedGraph g;
    bool r = g.havel_hakimi(d);

    cout << r << "\n";
    if(r)
        cout << g;
}*/

/*void Solve(){
    UnorientedGraph g;

    g.read_from_file("bridges");
    g.print_bridges;

}*/

void Solve(){
    /// find all biconnected components

    UnorientedGraph g;

    g.read_from_file("biconex.in");
    g.biconex();
}

int main()
{
    Solve();
    fout.close();

    return 0;
}

