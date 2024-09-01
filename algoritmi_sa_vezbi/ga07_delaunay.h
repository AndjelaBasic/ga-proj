#ifndef GA07_DELAUNAY_H
#define GA07_DELAUNAY_H

#include <vector>
using namespace std;

class Vertex;
class Face;
class Edge;
class TreeHistory;
/* za cvor nam je dovoljno da imamo koordinate i jednu poluivicu */
class Vertex{
public:
    int id=0;
    double x;
    double y;
    Edge* incidentEdge;
    Vertex(){}
    void printVertex() const;
    Vertex(double xx, double yy, Edge* e1):x(xx),y(yy),incidentEdge(e1){}
    bool operator == (const Vertex& v) const;
    bool operator > (const Vertex& v) const;
    double dist2(const Vertex& v) const;
    double norm2() const;

};

/* posto ovo strukturu koristimo za trouglove u Delaunay trijangulaciji
ona znamo da necemo imati rupe i dovoljno nam je da znamo jednu ivicu
sa spoljasnje granice */
class Face{
public:
    Edge* outerComponent;
    Face(){}
    void printFace() const;
    Face(Edge* oc):outerComponent(oc){}
    bool circumCircleContains(Vertex* v);
    bool inTriangle(Vertex* v) const;
    //inline std::vector<Face*> getChildrenFaces() const {return childrenFaces;}
    TreeHistory* my_history;


    //std::vector<Face*> childrenFaces;
};

/* struktura za poluivicu, ona je usmerena,
pa je origin cvor od koga krece, incidentFace je lice
koje nam stoji s leve strane kada idemo po usmerenoj ivici,
twin nam je poluivica s druge strane, za istu ivicu,
prev i next nam daju mogucnost obilaska sa obe strane datog lica */
class Edge{
public:
    Vertex* origin;
    Edge* twin;
    Edge* next;
    Edge* prev;
    Face* incidentFace;
    Edge(){}
    inline Vertex* getOrigin() const { return origin;}
    void printEdge() const;
    void setEdgeParameters(Vertex* o, Edge* t, Face* f, Edge* n, Edge* p);
    bool containsPoint(double x, double y) const;
    void flip();
    bool isLeft(Vertex* v) const;
    bool isBordering() const;

private:
    void help_triangle(Edge* a, Edge* b, Edge* c);
};

class TreeHistory{
public:
    TreeHistory(){}
    Face* pointingFace;
    std::vector<TreeHistory*> childrenFaces = std::vector<TreeHistory*>();
};


class DCEL_TRIANGLE{
public:
    DCEL_TRIANGLE(){}
    void add_point_in_triangle(Face* t, Vertex* v);
    void add_point_on_edge(Face* t, Edge* e, Vertex* v);
    void add_point_general(Face* t, Vertex* v);
    //void legalize_edge(Vertex* v, Edge* e);
    inline TreeHistory* getRootHistory() const {return root;}
    //Face* containingTriangle(TreeHistory* f, Vertex* v) const;
    TreeHistory* root=nullptr;
private:

    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;

};


#endif // GA07_DELAUNAY_H
