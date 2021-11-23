#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>
#include <tuple>
using namespace std;

ofstream fout("biconex.out");

class Graph{

protected:
    int n, m; /// n -> number of nodes / m -> number of edges
    vector<vector<int> > l; /// list of adjacency

    void do_dfs(const int start, vector<bool> &viz, const bool print) const;
    void do_bfs(const int start, vector<int> &dist, const bool print) const;

    virtual istream& read(istream&);
    virtual ostream& print(ostream&) const;

    static int findx (int, vector<int> &);
    static void unionx(const int, const int, vector<int> &, vector<int> &);

public:
    Graph(int n, bool dummy) : n(n) {} /// used to construct graph_with_costs
    Graph(int n = 1, int m = 0);
    Graph(const Graph &);
    Graph& operator= (const Graph &);

    int get_numEdges() const;
    int get_numNodes() const;

    friend istream& operator>>(istream&, Graph&);
    friend ostream& operator<<(ostream&, const Graph&);
    virtual void read_from_file(string) = 0;

    void dfs(const int) const;
    int num_of_components() const;
    void bfs(const int) const;
    vector<int> dist_from_node(const int) const;
    static void disjoint(const string, const string);

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

    this->l = g.l;
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
        l = g.l;
     }

     return *this;
}

istream& operator>>(istream& in, Graph& g){
    g.read(in);
    return in;
}

ostream& operator<<(ostream& out, const Graph& g){
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

vector<int> Graph :: dist_from_node (const int start) const {
    /// this function calculates the distance between
    /// the start node and every other node in the graph
    /// if the node is unreachable -> the distance is -1

    vector<int> dist(n + 1, 0);
    do_bfs(start, dist);

    return dist;
}

int Graph :: findx (int x, vector<int> &sets){
    /// method that returns the set that includes x

    int root = x;
    while(sets[root] != root)  /// go up on the tree that represents the set
        root = sets[root];

    while(sets[x] != root){
        int temp;
        temp = sets[x];
        sets[x] = root;
        x = temp;
    }

    return root;
}

void Graph :: unionx(const int x, const int y, vector<int> &sets, vector<int> &card){
    /// method that unifies 2 disjoint sets

    int root_x = findx(x, sets);
    int root_y = findx(y, sets); /// find the root of each disjoint set

    /// we will append the smallest set to the bigger set
    if(card[root_x] >= card[root_y]){
        card[root_x] += card[root_y];  /// update the size of the bigger set
        sets[root_y] = root_x;        /// now root_y is child of root_x
    }
    else{
        card[root_y] += card[root_x];
        sets[root_x] = root_y;
    }
}

void Graph :: disjoint(const string input_file, const string output_file){
    /// method that reads n operations from a file and prints the output to other file
    /// the operations can be 1 -> unify two disjoint sets
    ///                       2 -> verify if two elements are in the same set

    ifstream fin(input_file);
    ofstream fout(output_file);
    int n, m;

    fin >> n >> m;

    vector<int> sets(n + 1, 0), card(n + 1, 1);  /// sets -> vector of sets, card[i] = cardinal of set i

    for(int i = 1; i <= n; ++i) /// we begin with n sets with 1 element
        sets[i] = i;

    card[0] = 0;

    for(int i = 1; i <= m; ++i){
        int x, y, op;
        fin >> op >> x >> y;

         if(op == 1){
            if(findx(x, sets) != findx(y, sets))
                unionx(x, y, sets, card);   /// unify the 2 disjointed sets
        }
        else{
            if(findx(x, sets) == findx(y, sets))  /// if x and y are from same set print yes
                fout << "DA\n";
            else fout << "NU\n";
        }
    }

    fin.close();
    fout.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


class UnorientedGraph : public Graph{

public:
    UnorientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    UnorientedGraph(int n, int m, vector<pair<int, int> > v);
    UnorientedGraph& operator+=(const pair<int, int> edge);

    void read_from_file(string filename);
    vector<pair<int, int> > bridges() const;
    vector<vector<int> > biconex() const;
    bool havel_hakimi(const vector<int> v);

private:
    void h_merge(vector<pair<int, int> > &, const int ) const;
    void do_dfs_with_bridges(const int, int &, vector<bool> &, vector<int> &, vector<int> &, vector<int> &, vector<pair<int, int> > &) const;
    void do_biconex(const int, const int, vector<bool> &, vector<int> &, vector<int> &, stack<int> &, vector<vector<int> > &) const;
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

    UnorientedGraph g(n, 0);
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
void UnorientedGraph :: do_dfs_with_bridges(const int start, int &time, vector<bool> &viz, vector<int> &disc, vector<int> &low, vector<int> &parent, vector<pair<int, int> >& bridges) const {
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
            do_dfs_with_bridges(node, time, viz, disc, low, parent, bridges);

            low[start] = min(low[start], low[node]);
            if (low[node] > disc[start])
               bridges.push_back(make_pair(start, node));
        }
        else if(node != parent[start])  /// back-edge
            low[start]  = min(low[start], disc[node]);
    }

}

vector<pair<int, int> > UnorientedGraph :: bridges() const {
    /// Method that prints all bridges in a graph

    vector<bool> viz(n + 1, 0);
    vector<int> disc(n + 1);
    vector<int> low(n + 1);
    vector<int> parent(n + 1, -1);
    vector<pair<int, int> > bridges;

    int time = 0;

    for(int i = 1; i <= n; i++)
        if (!viz[i])
            do_dfs_with_bridges(i, time, viz, disc, low, parent, bridges);
    return bridges;
}

void UnorientedGraph :: h_merge (vector<pair<int, int> > &d, const int k) const {
    /// Utility function used by havel_hakimi to keep d
    /// sorted - merges first k elements with last n - k

    int i = 0;
    unsigned int j = k;
    vector<pair<int, int> > temp;

    while(i < k && j < d.size())
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

void UnorientedGraph :: do_biconex(const int son, const int dad, vector<bool> &viz, vector<int> &disc, vector<int> &low, stack<int> &s, vector<vector<int> > &comp) const {

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
            if(!viz[node]){
                s.push(node);
                do_biconex(node, son, viz, disc, low, s, comp);

                low[son] = min(low[son], low[node]);

                if(disc[son] <= low[node]){
                    s.push(son);

                    vector<int> temp;
                    while(!s.empty() && s.top() != node){
                        temp.push_back(s.top());
                        s.pop();
                    }
                    if(!s.empty()){
                        temp.push_back(s.top());
                        s.pop();
                    }
                    comp.push_back(temp);
                }
            }
            else /// back-edge
                if(disc[node] < low[son])
                    low[son] = disc[node];
        }
}

vector<vector<int> > UnorientedGraph :: biconex() const {

    /// Method that prints out the number of biconnected components
    /// and the biconnected components

    stack<int> s;
    vector<int> disc(n + 1, -1);
    vector<int> low(n + 1, -1);
    vector<bool> viz(n + 1, 0);
    vector<vector<int> > comp;


    disc[1] = 1;
    do_biconex(2, 1, viz, disc, low, s, comp);

    return comp;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


class OrientedGraph : public Graph{

public:
    OrientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    OrientedGraph(int n, int m, vector<pair<int, int> > v);
    OrientedGraph& operator+=(const pair<int, int> edge);
    vector<int> topological_sort() const;
    void read_from_file(string filename);
    OrientedGraph transpose_graph() const;
    vector<vector<int> > ctc() const;

private:
     inline void do_dfs(const int start, vector<bool> &viz, stack<int> &S) const;
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

void OrientedGraph :: read_from_file(string filename){
    /// Method to read a graph from a file filename

    int n, m;

    ifstream input;
    input.open(filename);

    input >> n >> m;

    OrientedGraph g(n, 0);
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

vector<int> OrientedGraph :: topological_sort() const {

    /// Method that prints out the nodes in topological order

    vector<bool> viz(n + 1, 0);
    stack<int> S;
    vector<int> order;

    for(int i = 1; i <= n; ++ i)
        if(!viz[i])
            do_dfs(i, viz, S);

    while(!S.empty()){
        order.push_back(S.top()); /// de schimbat in cout cand pun pe github
        S.pop();
    }

    return order;
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

vector<vector<int> > OrientedGraph :: ctc() const {

    /// method that prints out the number of strongly connected components
    /// and the strongly connected components of a graph

    OrientedGraph gt = transpose_graph();
    vector<bool> viz(n + 1, 0);
    stack<int> S;
    vector<vector<int> > sol;

    for(int i = 1; i <= n; ++ i)
        if(!viz[i]) do_dfs(i, viz, S);

    fill(viz.begin(), viz.end(), 0); /// reinitialize viz

    int ct = 0;
    string response = "";

    while(!S.empty()){
        int current_node = S.top();
        S.pop();
        if(!viz[current_node]){
            stack<int> component;

            ct ++;
            gt.do_dfs(current_node, viz, component);

            vector<int> temp;
            while(!component.empty()){
                temp.push_back(component.top());
                component.pop();
            }
            sol.push_back(temp);
        }
    }

    return sol;

}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------

class GraphWithCosts : public Graph {

    vector<vector<pair<int,int> > > l; /// same list of adjacency as in Graph but with costs
    bool isOriented;

private:
    //vector<Edge> edges; /// list of edges

public:
    GraphWithCosts(int n = 1, int m = 0, bool o = 0);
    GraphWithCosts(const GraphWithCosts &);
    GraphWithCosts& operator+=(const tuple<int, int, int>&);
    GraphWithCosts& operator= (const GraphWithCosts &);
    void read_from_file(string);
    vector<int> dijkstra(const int) const;
    vector<int> bellman_ford(const int) const;
    void kruskal(const string) const;


private:
    istream& read(istream&) override;
    virtual ostream& print(ostream&) const override;
    vector<tuple<int, int, int> > get_list_of_edges() const;

};

GraphWithCosts :: GraphWithCosts(int n, int m, bool o) : Graph(n, 0), l(n + 1) {
    Graph::m = m;
    isOriented = o;
};

GraphWithCosts :: GraphWithCosts(const GraphWithCosts &g) {
    /// copy constructor
    this->n = g.n;
    this->m = g.m;
    this->isOriented = g.isOriented;
    l = g.l;
    //edges = g.edges;
}


/// Operators
GraphWithCosts& GraphWithCosts :: operator= (const GraphWithCosts &g) {
    if(this != &g){
        this->n = g.n;
        this->m = g.m;
        this-> isOriented = g.isOriented;

        if(!this->l.empty())
            l.clear();
        l = g.l;

        //if(!this->edges.empty())
          //  edges.clear();
        //edges = g.edges;
     }

     return *this;
}

GraphWithCosts& GraphWithCosts :: operator+=(const tuple<int, int, int> &e) {
    /// add edge to graph

    m++;
    l[get<0>(e)].push_back(make_pair(get<1>(e), get<2>(e)));
    if(!isOriented)
        l[get<1>(e)].push_back(make_pair(get<0>(e), get<2>(e)));

    //edges.push_back(e); ///trebuie sa ma uit sa vad daca trebuie pusa si invers!!!!

    return *this;
}


/// read method
istream& GraphWithCosts :: read(istream& in) {
    Graph::read(in); /// call read from base
    GraphWithCosts temp(n, m, this->isOriented);

    for(int i = 1; i <= m; i ++){
        int x, y, c;
        cin >> x >> y >> c;

        temp.l[x].push_back(make_pair(y, c));
        if(!this->isOriented)
            temp.l[y].push_back(make_pair(x, c));
    }

    operator=(temp);
    return in;
}

void GraphWithCosts :: read_from_file(string filename) {
    /// Method to read a graph from a file filename

    int n, m;

    ifstream input;
    input.open(filename);

    input >> n >> m;

    GraphWithCosts g(n, 0, this->isOriented);
    for(int i = 1; i <= m; i ++){
        tuple<int, int, int> e;
        input >> get<0>(e) >> get<1>(e) >> get<2>(e);
        g += e;
        }

    *this = g;
    input.close();

}


/// methods
ostream& GraphWithCosts :: print(ostream& out) const {
    /// the default format is n m on one line followed
    /// by the number of the current node and a list of every
    /// reachable nodes

    cout << n << " " << m << "\n";

    for(int i = 1; i <= n; ++i){
        cout << i << ": ";
            for(auto node : l[i])
                cout << "(" << node.first << ", " << node.second << ") ";
        cout << "\n";
    }

    return out;
};

vector<int> GraphWithCosts :: dijkstra(const int s) const {
    /// method to find Dijkstra's shortest path using
    /// in the priority queue we store elements of pair
    /// first -> cost from the start node to the current node
    /// second -> node

    const int INF = 0x3f3f3f3f;
    priority_queue< pair<int, int>, vector <pair<int, int> > , greater<pair<int, int> > > pq;
    vector<int> dist(n + 1, INF);
    vector<bool>inHeap(n + 1, false);


    pq.push(make_pair(0, s));  /// the cost to reach start node is always 0
    dist[s] = 0;

    while (!pq.empty()){
        int node = pq.top().second;  /// node is the current node
        pq.pop();

        if (!inHeap[node]) {
            inHeap[node] = true;
            for(auto v : l[node])                             ///  iterate through all vertexes reachable from node
                if (dist[v.first] > dist[node] + v.second){   ///  if there is shorted path to v through node
                    dist[v.first] = dist[node] + v.second;    ///  update distance
                    pq.push(make_pair(dist[v.first], v.first));
            }
        }
    }

    return dist;
}

vector<int> GraphWithCosts :: bellman_ford(const int s) const {
    /// method to find shortest path using Bellman Ford algorithm
    /// optimization with a queue
    /// complexity N * M

    const int INF = 0x3f3f3f3f;
    bool noCycle = 1;                /// noCycle = 1 -> no negative costs cycles
    vector<int> dist(n + 1, INF);
    vector<bool> inQ(n + 1, false);
    vector<int> viz(n + 1, 0);      /// viz[i] = number of times node i was visited if > n => noCycle = 0
    queue<int> q;                   /// this queue is used to store nodes that could still make the cost smaller

    q.push(s);                      /// we begin only with the start node with cost 0
    dist[s] = 0;
    inQ[s] = 1;

    while(!q.empty() && noCycle){
        int node = q.front();      /// current node is q.front()
        q.pop();
        inQ[node] = 0;

        for(auto v : l[node])     /// iterate through all reachable vertexes
            if(dist[v.first] > dist[node] + v.second){   /// if we get a better cost
                dist[v.first] = dist[node] + v.second;   /// relax the edge
                viz[v.first] ++;

                if(!inQ[v.first]){       /// if we used the node we put it in the queue
                    q.push(v.first);
                    inQ[v.first] = 1;
                }

                if(viz[v.first] >= n)   /// test if we have a negative cycle
                    noCycle = 0;

            }
    }

    if(noCycle) return dist;

    vector<int> dummy;
    return dummy;    /// if there is a negative cycle we will return an empty list

}

vector<tuple<int, int, int> > GraphWithCosts :: get_list_of_edges() const {

    vector<tuple<int, int, int> > edges;

    for(int i = 1; i <= n; i ++)
        for(auto elem : l[i])
        if(isOriented == 1 || (!isOriented && i < elem.first)){
            tuple<int, int, int> e;
            get<0>(e) = elem.second;     /// put the cost first so we can sort by cost
            get<1>(e) = i;
            get<2>(e) = elem.first;
            edges.push_back(e);
    }

    return edges;

}

void GraphWithCosts :: kruskal (const string output_file) const {
    /// prints the mst of the graph in the output file using kruskal algorithm

    ofstream fout(output_file);

    vector<tuple<int, int, int> > edges = get_list_of_edges();
    vector<int> sets(n + 1, 0), card(n + 1, 1);       /// sets -> vector of disjointed sets, card[i] = size of set i
    vector<int> sol;                                  /// to store the index of the edges used for mst
    int cost = 0;                                     /// length and cost of solution

    card[0] = 0;
    for(int i = 1; i <= n; ++i)  /// at the start we did not choose any edge
        sets[i] = i;             /// so each node is part of a disjointed set

    sort(edges.begin(), edges.end());

    //for(auto t : edges)
      // cout<<get<0>(t)<<" "<<get<1>(t)<<" "<<get<2>(t)<<"\n";

    for(int i = 0; i < m && (int)sol.size() < n - 1; ++i){
        int x, y;

        x = get<1>(edges[i]);
        y = get<2>(edges[i]);

        if( findx(x,sets) != findx(y, sets) ){ /// if the selected nodes are not in the same set make union on sets
            unionx(x, y, sets, card);
            cost += get<0>(edges[i]);         /// update the cost of the mst
            sol.push_back(i);
        }
    }

    fout << cost << "\n" << sol.size() << "\n";
    for(int i = 0; i < (int)sol.size(); ++i){
        int idx = sol[i];
        fout << get<1>(edges[idx]) << " " << get<2>(edges[idx])<< "\n";
    }

    fout.close();
}


int main()
{
    GraphWithCosts g;
    g.read_from_file("apm.in");
    g.kruskal("apm.out");
    return 0;
}
