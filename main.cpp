#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>
#include <tuple>
#include <map>
using namespace std;

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

    static bool find_augumentation_path(vector<vector<int> >&, vector<int>&, vector<int>&, vector<int>&);
    static bool add_augumentation_path(int, vector<vector<int> >&, vector<int>&, vector<int>&, vector<int>&);

public:
    Graph(int n, bool dummy) : n(n) {} /// used to construct graph_with_costs
    Graph(int n = 1, int m = 0);
    Graph(const Graph &);
    Graph& operator= (const Graph &);

    int get_numEdges() const;
    int get_numNodes() const;
    vector<vector<int> > get_adjacency_matrix() const;

    friend istream& operator>>(istream&, Graph&);
    friend ostream& operator<<(ostream&, const Graph&);
    virtual void read_from_file(string) = 0;

    void dfs(const int) const;
    int num_of_components() const;
    void bfs(const int) const;
    vector<int> dist_from_node(const int) const;
    static void disjoint(const string, const string);
    vector<int> eulerian_cycle();
    static void hopcroft_karp(const string, const string);
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

vector<vector<int> > Graph :: get_adjacency_matrix() const {
    vector<vector<int> > m(n + 1, vector<int>(n + 1, 0));
    for(int i = 1; i <= n; i++)
        for(auto node : l[i])
            m[i][node];
    return m;
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
        cout << start << " ";

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

vector<int> Graph :: eulerian_cycle(){
    /// method that finds eulerian cycle in the given graph
    /// this method returns the cycles as a vector or
    /// a vector with one element (-1) is there is no such cycle

    vector<int> deg(n + 1);             /// deg[i] = degree of node i
    map <pair<int, int>, bool> h; /// h m

    for(int i = 1; i <= n; ++i){
        deg[i] += l[i].size();
        if(deg[i] & 1)
            return vector<int>(1, -1);  /// if the graph has eulerian cycle -> all degrees are even
        for(auto node : l[i])
            if(i <= node)
                h.insert(make_pair(make_pair(i, node), 0));
    }

    stack<int> s;      /// the algorithm is recursive but we want an in place implementation so we use a stack
    vector<int> res;   /// res -> result

    s.push(1);         /// start the cycle from first node
    while(!s.empty()){
        int node = s.top();

        if(l[node].empty()){
            s.pop();
            res.push_back(node);
        }
        else{
            int next = l[node].back(); /// take the last edge in the list
            l[node].pop_back();                     /// remove the edges

            pair<int, int> edge = make_pair(node, next);

            if(node >= next)
                edge.first = next, edge.second = node;
            cout<<edge.first<<"  "<<edge.second<<" "<<h[edge]<<"\n";
            if(!h[edge]){
                h[edge] = 1;
                s.push(next);
            }
        }

    }


    for(auto nodes_list : l)
        if(!nodes_list.empty())
            return vector<int>(1, -1); /// if there are edges left in the graph -> we don t have eulerian cycle

    return res;
}

void Graph :: hopcroft_karp(const string input, const string output){
    /// method that reads a bipartite graph from input file and does hopcroft-karp algorithm
    /// printing in the output file the maximum matching in the graph

    ifstream fin(input);
    ofstream fout(output);

    int n, m, nre;            /// n -> number of nodes in left side / m -> for right side / nre -> number of edges
    fin>>n>>m>>nre;

    vector<vector<int> > l(n + 1);                              /// adjacency list for the bipartite graph
    vector<int> left(n + 1, 0), right(m + 1, 0), dist(n + 1);  /// dist stores the distance of left side nodes -- sau m??
    int cnt = 0;  /// number of nodes in the maximum matching

    for(int i = 1; i <= nre; i ++){
        int x, y;
        fin>>x>>y;
        l[x].push_back(y);  /// we only need to store the edge from left to right
    }
    fin.close();

    while(find_augumentation_path(l, left, right, dist)){
        for(int i = 1; i <= n; i++) /// find a free node
            if(!left[i] && add_augumentation_path(i, l, left, right, dist))
                cnt++;

    }

    fout<<cnt<<"\n";
    for(int i = 1; i <= n; i ++)
        if(left[i])
            fout<<i<<" "<<left[i]<<"\n";

    fout.close();


}

bool Graph :: find_augumentation_path(vector<vector<int> > &l, vector<int> &left, vector<int> &right, vector<int> &dist){

    queue<int> q;
    const int INF = 0x3f3f3f3f;
    /// first layer of nodes with dist = 0
    for(int i = 1; i < (int)left.size(); i++)
        if(!left[i]){
            /// node i is not matched
            dist[i] = 0;
            q.push(i);
        }
        else dist[i] = INF;

    dist[0] = INF;

    while(!q.empty()){
        int node = q.front();
        q.pop();

        if(dist[node] < dist[0])
            for(auto next : l[node])
                if(dist[right[next]] == INF){
                    dist[right[next]] = dist[node] + 1;
                    q.push(right[next]);
                }

    }

    return dist[0] != INF;
}

bool Graph :: add_augumentation_path(int start, vector<vector<int> > &l, vector<int> &left, vector<int> &right, vector<int> &dist){
   const int INF = 0x3f3f3f3f;;
    if(!start) return true;

    for(auto i = l[start].begin(); i != l[start].end(); ++i){
        int v = *i;
        if(dist[right[v]] == dist[start]+1)
            if(add_augumentation_path(right[v], l, left, right, dist)){
                right[v] = start;
                left[start] = v;
                return 1;
                }
        }

        dist[start] = INF;
        return 0;

}



//------------------------------------------------------------------------------------------- ----------------------------------------------------------------------


class UnorientedGraph : public Graph{

public:
    UnorientedGraph(int n = 1, int m = 0) : Graph(n, m) {};
    UnorientedGraph(int n, int m, vector<pair<int, int> > v);
    UnorientedGraph& operator+=(const pair<int, int> edge);

    void read_from_file(string filename);
    vector<pair<int, int> > bridges() const;
    vector<vector<int> > biconex() const;
    bool havel_hakimi(const vector<int> v);
    int tree_diameter();
    void euler_infoarena(string, string);

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

int UnorientedGraph :: tree_diameter(){
    int node = 1;
    int maxi = 0;

    vector<int> dist = dist_from_node(node);

    for(int i = 1; i <= n; i ++)
        if(dist[i] > maxi){
            maxi = dist[i];
            node = i;
        }

    dist = dist_from_node(node);

     for(int i = 1; i <= n; i ++)
        if(dist[i] > maxi)
            maxi = dist[i];
    return maxi;

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
    vector<vector<int> > matrix;       /// matrix of adjacency that will be used only for roy_floyd
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

    void generate_matrix();
    void read_matrix_from_file(string);
    void erase_matrix();
    vector<vector<int> > roy_floyd();


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
     }

     return *this;
}

GraphWithCosts& GraphWithCosts :: operator+=(const tuple<int, int, int> &e) {
    /// add edge to graph

    m++;
    l[get<0>(e)].push_back(make_pair(get<1>(e), get<2>(e)));
    if(!isOriented)
        l[get<1>(e)].push_back(make_pair(get<0>(e), get<2>(e)));

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

void GraphWithCosts :: generate_matrix () {
     matrix = vector<vector<int> >(n, vector<int>(n, 0));

    for(int i = 1; i <= n; i++)
        for(auto edge : l[i])
        matrix[i - 1][edge.first - 1] = edge.second;


}

void GraphWithCosts :: read_matrix_from_file (string filename) {

    ifstream input;
    input.open(filename);

    input>>n;
    matrix = vector<vector<int> >(n, vector<int>(n, 0));

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            input>>matrix[i][j];


}

void GraphWithCosts :: erase_matrix () {
    matrix.clear();
    matrix.resize(0);
}

vector<vector<int> > GraphWithCosts :: roy_floyd () {
        /// this method will perform the Roy Floyd algorithm for the
        /// given graph. If the graph does not have a matrix of adjacency
        /// this method will call generate_matrix

        if(matrix.empty()) generate_matrix();

        const int INF = 0x3f3f3f3f;
        vector<vector<int> > dist(n, vector<int>(n, INF));

        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(matrix[i][j])
                    dist[i][j] = matrix[i][j];
                else if(i == j)
                    dist[i][j] = 0;

        for (auto k = 0; k < n; k++) { /// Pick every pair of vertiges (i, j)
            for (auto i = 0; i < n; i++) { /// try to improve the cost using node k
                for (auto j = 0; j < n; j++){
                    if (dist[i][j] > dist[i][k] + dist[k][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }


    return dist;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------


struct Edge{
    int cost, cap, used, dst;
};

class Network : public Graph{
    vector<vector<Edge> > l;

public:
    Network(int n = 1, int m = 0) : Graph(n, 0), l(n + 1) {}
    Network(const Network &);
    Network& operator+=(const pair<int, Edge>&);
    Network& operator= (const Network &);
    void read_from_file(string);
    int maxflow(int s, int d) const;

private:
    istream& read(istream&) override;
    ostream& print(ostream&) const override;
    bool bfs_with_flow(const int, const int, vector<int> &, vector<int> &, vector<vector<int> > &, vector<vector<int> > &, vector<vector<int> > &) const;

};


Network :: Network(const Network &g) {
    /// copy constructor
    this->n = g.n;
    this->m = g.m;
    l = g.l;
}

Network& Network :: operator= (const Network &g) {
    if(this != &g){
        this->n = g.n;
        this->m = g.m;

        if(!this->l.empty())
            l.clear();
        l = g.l;
     }

     return *this;
}

Network& Network :: operator+=(const pair<int, Edge> &e) {
    /// add edge to graph

    m++;
    l[e.first].push_back(e.second);

    return *this;
}


void Network :: read_from_file(string filename){
    int n, m;

    ifstream input;
    input.open(filename);

    input >> n >> m;

    Network g(n, 0);
    for(int i = 1; i <= m; i ++){
        Edge e;
        int src;
        input >> src >> e.dst >> e.cap;
        g += make_pair(src, e);
        }

    *this = g;
    input.close();


}

istream& Network :: read(istream& in) {
    Graph::read(in); /// call read from base
    Network temp(n, m);

    for(int i = 1; i <= m; i ++){
        int src;
        Edge e;
        cin >> src >> e.dst >> e.cap;

        temp.l[src].push_back(e);
    }

    operator=(temp);
    return in;
}


/// methods
ostream& Network :: print(ostream& out) const {
    /// the default format is n m on one line followed
    /// by the number of the current node and a list of every
    /// reachable nodes

    cout << n << " " << m << "\n";

    for(int i = 1; i <= n; ++i){
        cout << i << ": ";
            for(auto node : l[i])
                cout << "(" << node.dst << ", " << node.cap << ") ";
        cout << "\n";
    }

    return out;
};

bool Network :: bfs_with_flow(const int s, const int d, vector<int>& parent, vector<int> &viz, vector<vector<int> > &res, vector<vector<int> > &capacity, vector<vector<int> > &flux) const{
    /// this does bfs from s and stops when reaching node d
    /// the function returns the flow of the path and
    /// parent vector holds the bfs tree

    //fill(parent.begin(), parent.end(), 0); /// reinitialize the parent vector
    parent[s] = -1;
    parent[d] = 0;

    viz.clear();
    viz.resize(n+1, 0);
    viz[s] = 1;

    queue<int> q;
    q.push(s);

    while(!q.empty() && !parent[d]){
        int node = q.front();
        q.pop();

        for(int next : res[node])
            if(!viz[next] && capacity[node][next] > flux[node][next]){ /// if the node is new and the edge is not saturated
                parent[next] = node;
                viz[next] = 1;
                q.push(next);
            }

    }

    return parent[d];
}

int Network :: maxflow(int s, int d) const{
    /// this is an implementation of Edmonds Karp algorithm
    /// to determine the max flow from s (source) to d (destination)
    /// res - residual graph
    /// capacity - matrix that stores the residual capacity of the network
    /// flux - matrix that stores how much flow actually goes through an edge
    /// parent - vector that stores bfs tree
    /// https://cp-algorithms.com/graph/edmonds_karp.html?fbclid=IwAR2yev5zuz7lSFDE0H9uXDxSsEIc5nRIAQ7b8glHiAVCWUorK4mlfKhiyds

    vector<vector<int> > res(n + 1);
    vector<int> parent(n + 1, 0);
    vector<vector<int> > capacity(n + 1, vector<int>(n + 1, 0));
    vector<vector<int> > flux(n + 1, vector<int>(n + 1, 0));
    vector<int> viz(n + 1, 0);

    const int INF = 0x3f3f3f3f;
    int flow = 0;
    int new_flow;

    /// initializing
    for(int i = 1; i <= n; i++){
        for(auto elem : l[i]){
            int x = i;
            int y = elem.dst;

            capacity[x][y] =  elem.cap;
            res[x].push_back(y);
            res[y].push_back(x);
        }
    }

    /// Edmonds-Karp
    while(bfs_with_flow(s, d, parent, viz, res, capacity, flux))
        for(auto node : res[d])
            if(viz[node] && flux[node][d] != capacity[node][d]){
                parent[d] = node;
                new_flow = INF;
                for(int i = d; i != s; i = parent[i]){
                    int p = parent[i];
                    if(capacity[p][i] - flux[p][i] < new_flow)
                        new_flow = capacity[p][i] - flux[p][i];
                }

                if(new_flow){
                    /// go up the bfs tree to update capacities
                    for(int i = d; i != s; i = parent[i]){
                        int p = parent[i];
                        flux[p][i] += new_flow;
                        flux[i][p] -= new_flow;
                    }
                }
                flow += new_flow;
            }


    return flow;

}

int main()
{
    Graph::hopcroft_karp("cuplaj.in", "cuplaj.out");
    return 0;
}
