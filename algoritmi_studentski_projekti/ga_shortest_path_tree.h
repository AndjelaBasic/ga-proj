#ifndef GA_SHORTEST_PATH_TREE_H
#define GA_SHORTEST_PATH_TREE_H


#include <vector>
#include <map>
#include <set>
#include "algoritambaza.h"

using namespace std;

class Vertex_SPT;
class Face;
class Edge;
class TreeHistory;


/* za cvor nam je dovoljno da imamo koordinate i jednu poluivicu */
class Vertex_SPT {
public:
    int id = 0;
    double x;
    double y;
    Edge* incidentEdge;
    Vertex_SPT() {}
    void printVertex_SPT() const;
    Vertex_SPT(double xx, double yy, Edge* e1) :x(xx), y(yy), incidentEdge(e1) {}
    bool operator == (const Vertex_SPT& v) const;
    bool operator > (const Vertex_SPT& v) const;
    bool operator < (const Vertex_SPT& v) const;
    double dist2(const Vertex_SPT& v) const;
    double norm2() const;
    QPointF toQPointF();

};

/* posto ovo strukturu koristimo za trouglove u Delaunay trijangulaciji
ona znamo da necemo imati rupe i dovoljno nam je da znamo jednu ivicu
sa spoljasnje granice */
class Face {
public:
    Edge* outerComponent;
    Face() {}
    void printFace() const;
    Face(Edge* oc) :outerComponent(oc) {}
    bool circumCircleContains(Vertex_SPT* v);
    bool inTriangle(Vertex_SPT* v) const;
    TreeHistory* my_history;
};

/* struktura za poluivicu, ona je usmerena,
pa je origin cvor od koga krece, incidentFace je lice
koje nam stoji s leve strane kada idemo po usmerenoj ivici,
twin nam je poluivica s druge strane, za istu ivicu,
prev i next nam daju mogucnost obilaska sa obe strane datog lica */
class Edge {
public:
    Vertex_SPT* origin;
    Edge* twin;
    Edge* next;
    Edge* prev;
    Face* incidentFace;
    Edge() {}
    inline Vertex_SPT* getOrigin() const { return origin; }
    void printEdge() const;
    void setEdgeParameters(Vertex_SPT* o, Edge* t, Face* f, Edge* n, Edge* p);
    bool containsPoint(double x, double y) const;
    void flip();
    bool isLeft(Vertex_SPT* v) const;
    bool isBordering() const;

private:
    void help_triangle(Edge* a, Edge* b, Edge* c);
};

class TreeHistory {
public:
    TreeHistory() {}
    Face* pointingFace;
    vector<TreeHistory*> childrenFaces = vector<TreeHistory*>();
    void getDelaunayGraph(map<Vertex_SPT*, set<Vertex_SPT*>> &delaunayGraph);
};


class DCEL_TRIANGLE {
public:
    DCEL_TRIANGLE() {}
    void add_point_in_triangle(Face* t, Vertex_SPT* v);
    void add_point_on_edge(Face* t, Edge* e, Vertex_SPT* v);
    void add_point_general(Face* t, Vertex_SPT* v);
    inline TreeHistory* getRootHistory() const { return root; }
    TreeHistory* root = nullptr;
};

class DelaunayUtil{

public:
    DelaunayUtil() {}
    Face* getSuperTriangle(std::vector<Vertex_SPT*> &in_points);
    Face* containingTriangle(TreeHistory* f, Vertex_SPT* v);
    bool testTriangle(Face* f, Vertex_SPT* v);
    void legalize_edge(Vertex_SPT* v, Edge* e);
    Face* CopyOfTriangle(Face* triangle);
};

/* Klasa koja predstavlja algoritam */
class UnweightedShortestPathTree : public AlgoritamBaza {
public:
    /* Konstruktor i destruktor klase */
    UnweightedShortestPathTree(QWidget *, int,
                        const bool & = false,
                        std::string = "",
                        int = BROJ_SLUCAJNIH_OBJEKATA);
    virtual ~UnweightedShortestPathTree() override;

    /* Virtuelni metodi iz natklase */
    void pokreniAlgoritam() final;
    void crtajAlgoritam(QPainter *) const final;
    void pokreniNaivniAlgoritam() final;
    void crtajNaivniAlgoritam(QPainter *) const final;

    int _n = 0;
    vector<Vertex_SPT*> P = vector<Vertex_SPT*>();
    Vertex_SPT *s = nullptr;
    map<Vertex_SPT*, set<Vertex_SPT*>> delaunayGraph;
    map<Vertex_SPT*, set<Vertex_SPT*>> unitDiscGraph;
    map<Vertex_SPT*, pair <int, Vertex_SPT*>> resultPathGlavni;
    map<Vertex_SPT*, pair <int, Vertex_SPT*>> resultPathNaivni;

private:
    /* Rad sa podacima, inicijalizacija */
    void generisiNasumicnaTemena(int);
    void ucitajPodatkeIzDatoteke(std::string);
    void inicijalizujDeloneGraf();
    void inicijalizujUnitDiscGraf();
};


// kd-tree node structure
struct KdNode {
    Vertex_SPT point;
    KdNode* left;
    KdNode* right;
    int splitDim;

    KdNode(const Vertex_SPT& p, KdNode* l = nullptr, KdNode* r = nullptr, int d = 0)
        : point(p), left(l), right(r), splitDim(d) {}
};


class KdTree {

public:
    KdTree() {}
    KdNode* buildKdTree(vector<Vertex_SPT*>& vertices, int start, int end, int depth);
    Vertex_SPT* findNearestNeighbor(KdNode* root, const Vertex_SPT& query, double& minDist);

};
#endif // GA_SHORTEST_PATH_TREE_H
