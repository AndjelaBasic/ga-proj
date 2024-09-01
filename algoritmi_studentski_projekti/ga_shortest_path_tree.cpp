#include "ga_shortest_path_tree.h"
#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <utility>
#include <fstream>
#include <algorithm>

#define epsilon 0.0001

DCEL_TRIANGLE d_object;

DelaunayUtil util;

KdTree kdTree;

using namespace std;

// util functions
void DelaunayUtil::legalize_edge(Vertex_SPT* v, Edge* e) {

    if (e->isBordering())
        return;

    Vertex_SPT* p = (e->prev)->origin;

    Edge* considered_edge = e;
    if ((*p) == (*v)) {
        p = ((e->twin)->prev)->origin;
        considered_edge = e->twin;
    }

    Face* my_triangle = considered_edge->incidentFace;
    Edge* next_to_legalize1 = considered_edge->next;
    Edge* next_to_legalize2 = considered_edge->prev;

    if (my_triangle->circumCircleContains(v)) {
        considered_edge->flip();
        legalize_edge(v, next_to_legalize1);
        legalize_edge(v, next_to_legalize2);
    }
}

bool DelaunayUtil::testTriangle(Face* f, Vertex_SPT* v) {
    Edge* e1 = f->outerComponent;
    Edge* e2 = (f->outerComponent)->next;
    Edge* e3 = ((f->outerComponent)->next)->next;
    return (f->inTriangle(v) || e1->containsPoint(v->x, v->y) || e2->containsPoint(v->x, v->y) || e3->containsPoint(v->x, v->y));
}

Face* DelaunayUtil::containingTriangle(TreeHistory* f, Vertex_SPT* v) {
    if ((f->childrenFaces).size() == 0) {
        return f->pointingFace;
    }
    else if ((f->childrenFaces).size() == 2) {
        if (testTriangle((f->childrenFaces)[0]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[0], v);
        else
            return containingTriangle((f->childrenFaces)[1], v);
    }
    else {
        if (testTriangle((f->childrenFaces)[0]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[0], v);
        else if (testTriangle((f->childrenFaces)[1]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[1], v);
        else
            return containingTriangle((f->childrenFaces)[2], v);
    }
}

Face* DelaunayUtil::getSuperTriangle(std::vector<Vertex_SPT*> &in_points) {

    double miny = in_points[0]->y;
    double minslope = std::numeric_limits<double>::infinity();

    for(int i=1; i<in_points.size(); i++){

        // calculate smallest y
        if (in_points[i]->y < miny)
            miny = in_points[i]->y;

        if (*in_points[i] > *in_points[0]) {

            // swap such that first element is the largest
            Vertex_SPT* tmp = in_points[0];
            in_points[0] = in_points[i];
            in_points[i] = tmp;

            // set ids such that largest element has id = 0
            in_points[0]->id = 0;
            in_points[i]->id = i;
        }
    }
    // calculate min slope
    for (int i = 1; i < in_points.size(); i++) {
        double slope = fabs((in_points[0]->y - in_points[i]->y) / (in_points[0]->x - in_points[i]->x + epsilon));
        if (slope < minslope)
            minslope = slope;
    }

    minslope = 0.8 * minslope; // further reduce minslope by 20%

    double t = miny - 5;  // all points are above y = t

    Vertex_SPT* pminus2 = new Vertex_SPT(in_points[0]->x + (t- in_points[0]->y)/minslope,t, NULL);
    pminus2->id = -2;

    Vertex_SPT* pminus1 = new Vertex_SPT(in_points[0]->x + (t - in_points[0]->y) / (-minslope), t, NULL);
    pminus1->id = -1;

    // set beginning edges and triangle

    Edge* e1 = new Edge();
    Edge* e2 = new Edge();
    Edge* e3 = new Edge();
    Edge* e1_twin = new Edge();
    Edge* e2_twin = new Edge();
    Edge* e3_twin = new Edge();

    Face* beginTriangle = new Face(e1);


    e1->setEdgeParameters(in_points[0], e1_twin, beginTriangle, e2, e3);
    e2->setEdgeParameters(pminus2, e2_twin, beginTriangle, e3, e1);
    e3->setEdgeParameters(pminus1, e3_twin, beginTriangle, e1, e2);
    e1_twin->setEdgeParameters(pminus2, e1, NULL, e3_twin, e2_twin);
    e2_twin->setEdgeParameters(pminus1, e2, NULL, e1_twin, e3_twin);
    e3_twin->setEdgeParameters(in_points[0], e3, NULL, e2_twin, e1_twin);


    return beginTriangle;
}

Face* DelaunayUtil::CopyOfTriangle(Face* triangle) {
    Edge* e1_old = triangle->outerComponent;
    Edge* e2_old = (triangle->outerComponent)->next;
    Edge* e3_old = ((triangle->outerComponent)->next)->next;

    Vertex_SPT* a = e1_old->getOrigin();
    Vertex_SPT* b = e2_old->getOrigin();
    Vertex_SPT* c = e3_old->getOrigin();

    Edge* e1 = new Edge(); Edge* e1_twin = new Edge();
    Edge* e2 = new Edge(); Edge* e2_twin = new Edge();
    Edge* e3 = new Edge(); Edge* e3_twin = new Edge();

    Vertex_SPT* v1_new = new Vertex_SPT(a->x, a->y, e1);
    Vertex_SPT* v2_new = new Vertex_SPT(b->x, b->y, e2);
    Vertex_SPT* v3_new = new Vertex_SPT(c->x, c->y, e3);
    v1_new->id = a->id;
    v2_new->id = b->id;
    v3_new->id = c->id;

    Face* d = new Face(e1);
    e1->setEdgeParameters(v1_new, e1_twin, d, e2, e3);
    e2->setEdgeParameters(v2_new, e2_twin, d, e3, e1);
    e3->setEdgeParameters(v3_new, e3_twin, d, e1, e2);

    e1_twin->setEdgeParameters(v2_new, e1, d, e3_twin, e2_twin);
    e2_twin->setEdgeParameters(v3_new, e2, d, e1_twin, e3_twin);
    e3_twin->setEdgeParameters(v1_new, e3, d, e2_twin, e1_twin);

    return d;
}

// Function to build the kd-tree
KdNode* KdTree::buildKdTree(vector<Vertex_SPT*>& vertices, int start, int end, int depth) {
    if (start > end) return nullptr;

    int mid = (start + end) / 2;
    int splitDim = depth % 2; // 0 for x, 1 for y

    // Sort the vertices based on the current split dimension
    sort(vertices.begin() + start, vertices.begin() + end + 1,
         [splitDim](const Vertex_SPT* a, const Vertex_SPT* b) {
             return splitDim == 0 ? a->x < b->x : a->y < b->y;
         });

    KdNode* root = new KdNode(*vertices[mid], buildKdTree(vertices, start, mid - 1, depth + 1),
                              buildKdTree(vertices, mid + 1, end, depth + 1), splitDim);
    return root;
}


// Function to find the nearest neighbor
Vertex_SPT* KdTree::findNearestNeighbor(KdNode* root, const Vertex_SPT& query, double& minDist) {
    if (!root) return nullptr;

    double dist = root->point.dist2(query);
    Vertex_SPT* bestPoint = nullptr;

    // Check if the current point is closer
    if (dist < minDist) {
        minDist = dist;
        bestPoint = &root->point;
    }

    int splitDim = root->splitDim;
    KdNode* nearSubtree = (splitDim == 0 ? query.x <= root->point.x : query.y <= root->point.y)
                              ? root->left
                              : root->right;
    KdNode* farSubtree = (splitDim == 0 ? query.x <= root->point.x : query.y <= root->point.y)
                             ? root->right
                             : root->left;

    // Recursively search the nearer subtree
    Vertex_SPT* temp = findNearestNeighbor(nearSubtree, query, minDist);
    if (temp != nullptr) {
        bestPoint = temp;
    }

    // Check if we need to search the farther subtree
    double splitDist = splitDim == 0 ? (query.x - root->point.x) : (query.y - root->point.y);
    if (splitDist * splitDist < minDist) {
        temp = findNearestNeighbor(farSubtree, query, minDist);
        if (temp != nullptr) {
            bestPoint = temp;
        }
    }

    return bestPoint;
}

void TreeHistory::getDelaunayGraph(map<Vertex_SPT*, set<Vertex_SPT*>> &delaunayGraph) {

    // map where keys are start vertices and values are vecors of end vertices for edges in Delaunay triangulation

    if ((this->childrenFaces).size() == 0) {

        Vertex_SPT* A1 = ((this->pointingFace)->outerComponent)->origin;
        Vertex_SPT * A2 = (((this->pointingFace)->outerComponent)->twin)->origin;
        if (A1->id >=0 && A2->id >= 0)
            delaunayGraph[A1].insert(A2);

        Vertex_SPT* B1 = (((this->pointingFace)->outerComponent)->next)->origin;
        Vertex_SPT* B2 = ((((this->pointingFace)->outerComponent)->next)->twin)->origin;
        if (B1->id >= 0 && B2->id >= 0)
            delaunayGraph[B1].insert(B2);

        Vertex_SPT* C1 = ((((this->pointingFace)->outerComponent)->next)->next)->origin;
        Vertex_SPT* C2 = (((((this->pointingFace)->outerComponent)->next)->next)->twin)->origin;
        if (C1->id >= 0 && C2->id >= 0)
            delaunayGraph[C1].insert(C2);
    }

    if ((this->childrenFaces).size() == 2) {
        (*(this->childrenFaces)[0]).getDelaunayGraph(delaunayGraph);
        (*(this->childrenFaces)[1]).getDelaunayGraph(delaunayGraph);
    }

    if ((this->childrenFaces).size() == 3) {
        (*(this->childrenFaces)[0]).getDelaunayGraph(delaunayGraph);
        (*(this->childrenFaces)[1]).getDelaunayGraph(delaunayGraph);
        (*(this->childrenFaces)[2]).getDelaunayGraph(delaunayGraph);
    }

}

void DCEL_TRIANGLE::add_point_general(Face* t, Vertex_SPT* v) {

    Face* triangle = t;
    Edge* e1 = triangle->outerComponent;
    Edge* e2 = (triangle->outerComponent)->next;
    Edge* e3 = ((triangle->outerComponent)->next)->next;

    double x = v->x;
    double y = v->y;

    if (e1->isBordering() && e2->isBordering() && e3->isBordering()) {
        add_point_in_triangle(triangle, v);
    }
    else if (e1->containsPoint(x, y)) {
        add_point_on_edge(triangle, e1, v);
    }
    else if (e2->containsPoint(x, y)) {
        add_point_on_edge(triangle, e2, v);
    }
    else if (e3->containsPoint(x, y)) {
        add_point_on_edge(triangle, e3, v);
    }
    else {
        add_point_in_triangle(triangle, v);
    }

}


void DCEL_TRIANGLE::add_point_in_triangle(Face* triangle, Vertex_SPT* v) {

    Face* t = util.CopyOfTriangle(triangle);
    TreeHistory* history = triangle->my_history;

    Edge* e1 = triangle->outerComponent;
    Edge* e2 = (triangle->outerComponent)->next;
    Edge* e3 = ((triangle->outerComponent)->next)->next;

    Vertex_SPT* a = e1->getOrigin();
    Vertex_SPT* b = e2->getOrigin();
    Vertex_SPT* c = e3->getOrigin();


    Edge* a1 = new Edge(); Edge* a2 = new Edge();
    Edge* b1 = new Edge(); Edge* b2 = new Edge();
    Edge* c1 = new Edge(); Edge* c2 = new Edge();

    v->incidentEdge = a2;

    Face* f1 = new Face(a2);
    Face* f2 = new Face(b2);
    Face* f3 = new Face(c2);

    e1->next = b1; e1->prev = a2; e1->incidentFace = f1;
    e2->next = c1; e2->prev = b2; e2->incidentFace = f2;
    e3->next = a1; e3->prev = c2; e3->incidentFace = f3;

    a1->setEdgeParameters(a, a2, f3, c2, e3);
    a2->setEdgeParameters(v, a1, f1, e1, b1);
    b1->setEdgeParameters(b, b2, f1, a2, e1);
    b2->setEdgeParameters(v, b1, f2, e2, c1);
    c1->setEdgeParameters(c, c2, f2, b2, e2);
    c2->setEdgeParameters(v, c1, f3, e3, a1);


    history->pointingFace = t;
    (history->pointingFace)->my_history = history;

    TreeHistory* newFace1 = new TreeHistory();
    newFace1->pointingFace = f1;
    f1->my_history = newFace1;
    (history->childrenFaces).push_back(newFace1);


    TreeHistory* newFace2 = new TreeHistory();
    newFace2->pointingFace = f2;
    f2->my_history = newFace2;
    (history->childrenFaces).push_back(newFace2);


    TreeHistory* newFace3 = new TreeHistory();
    newFace3->pointingFace = f3;
    f3->my_history = newFace3;
    (history->childrenFaces).push_back(newFace3);

    util.legalize_edge(v, e1);
    util.legalize_edge(v, e2);
    util.legalize_edge(v, e3);

}

void DCEL_TRIANGLE::add_point_on_edge(Face* triangle, Edge* e, Vertex_SPT* v) {
    Face* triangle_incident = (e->twin)->incidentFace;

    TreeHistory* history_incident = triangle_incident->my_history;

    TreeHistory* history = triangle->my_history;


    Edge* e_next = e->next;
    Edge* e_prev = e->prev;
    Edge* e_twin = e->twin;

    Edge* e1 = new Edge();
    Edge* e1_twin = new Edge();
    Edge* e2 = new Edge();
    Edge* e2_twin = new Edge();
    Edge* e3 = new Edge();
    Edge* e3_twin = new Edge();
    Edge* e4 = new Edge();
    Edge* e4_twin = new Edge();

    Face* f1 = new Face(e1);
    Face* f2 = new Face(e4);
    Face* f3 = new Face(e2);
    Face* f4 = new Face(e2_twin);

    v->incidentEdge = e2;

    e_next->next = e3_twin;
    e_next->prev = e2;
    e_next->incidentFace = f3;
    e_prev->next = e1;
    e_prev->prev = e3;
    e_prev->incidentFace = f1;

    (e_twin->next)->next = e4;
    (e_twin->next)->prev = e1_twin;
    (e_twin->next)->incidentFace = f2;
    (e_twin->prev)->next = e2_twin;
    (e_twin->prev)->prev = e4_twin;
    (e_twin->prev)->incidentFace = f4;


    e1->setEdgeParameters(e->origin, e1_twin, f1, e3, e_prev);
    e1_twin->setEdgeParameters(v, e1, f2, e_twin->next, e4);
    e2->setEdgeParameters(v, e2_twin, f3, e_next, e3_twin);
    e2_twin->setEdgeParameters(e_twin->origin, e2, f4, e4_twin, e_twin->prev);
    e3->setEdgeParameters(v, e3_twin, f1, e_prev, e1);
    e3_twin->setEdgeParameters(e_prev->origin, e3, f3, e2, e_next);
    e4->setEdgeParameters((e_twin->prev)->origin, e4_twin, f2, e1_twin, e_twin->next);
    e4_twin->setEdgeParameters(v, e4, f4, e_twin->prev, e2_twin);


    history->pointingFace = util.CopyOfTriangle(triangle);
    history_incident->pointingFace = util.CopyOfTriangle(triangle_incident);
    (history->pointingFace)->my_history = history;
    (history_incident->pointingFace)->my_history = history_incident;

    //Face* f = new Face(*f1);
    TreeHistory* newFace1 = new TreeHistory();
    newFace1->pointingFace = f1;
    f1->my_history = newFace1;
    (history->childrenFaces).push_back(newFace1);

    //f = new Face(*f3);
    TreeHistory* newFace3 = new TreeHistory();
    newFace3->pointingFace = f3;
    f3->my_history = newFace3;
    (history->childrenFaces).push_back(newFace3);

    //f = new Face(*f2);
    TreeHistory* newFace2 = new TreeHistory();
    newFace2->pointingFace = f2;
    f2->my_history = newFace2;
    (history_incident->childrenFaces).push_back(newFace2);

    //f = new Face(*f4);
    TreeHistory* newFace4 = new TreeHistory();
    newFace4->pointingFace = f4;
    f4->my_history = newFace4;
    (history_incident->childrenFaces).push_back(newFace4);

    util.legalize_edge(v, e_next);
    util.legalize_edge(v, e_prev);
    util.legalize_edge(v, e_twin->next);
    util.legalize_edge(v, e_twin->prev);

}

void Edge::flip() {

    Edge* a2 = this->twin;

    Edge* e1 = this->next;
    Edge* e2 = this->prev;

    Edge* e3 = a2->next;
    Edge* e4 = a2->prev;

    Vertex_SPT* v1 = e1->origin;
    Vertex_SPT* v2 = e4->origin;
    Vertex_SPT* v3 = this->origin;
    Vertex_SPT* v4 = e2->origin;

    TreeHistory* history = (this->incidentFace)->my_history;
    TreeHistory* history_incident = (twin->incidentFace)->my_history;

    history->pointingFace = util.CopyOfTriangle(this->incidentFace);
    history_incident->pointingFace = util.CopyOfTriangle(a2->incidentFace);
    (history->pointingFace)->my_history = history;
    (history_incident->pointingFace)->my_history = history_incident;

    Face* A = this->incidentFace;
    Face* B = a2->incidentFace;
    this->setEdgeParameters(v2, a2, A, e2, e3);
    a2->setEdgeParameters(v4, this, B, e4, e1);

    v1->incidentEdge = e1;
    v3->incidentEdge = e3;

    A->outerComponent = e2;
    B->outerComponent = e1;

    e1->incidentFace = B;
    e4->incidentFace = B;
    e3->incidentFace = A;
    e2->incidentFace = A;

    help_triangle(this, e2, e3);
    help_triangle(a2, e4, e1);

    TreeHistory* hA = new TreeHistory();
    hA->pointingFace = A;
    A->my_history = hA;

    TreeHistory* hB = new TreeHistory();
    hB->pointingFace = B;
    B->my_history = hB;

    (history->childrenFaces).push_back(hA);
    (history->childrenFaces).push_back(hB);
    (history_incident->childrenFaces).push_back(hA);
    (history_incident->childrenFaces).push_back(hB);
}

void Edge::help_triangle(Edge* a, Edge* b, Edge* c) {
    a->next = b; b->next = c; c->next = a;
    c->prev = b; b->prev = a; a->prev = c;
}

double Vertex_SPT::dist2(const Vertex_SPT& v) const {
    return (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y);
}
double Vertex_SPT::norm2() const {
    return x * x + y * y;
}

QPointF Vertex_SPT::toQPointF(){
    return QPointF(x, y);
}

bool Face::circumCircleContains(Vertex_SPT* v) {
    Edge* e1 = outerComponent;
    Edge* e2 = outerComponent->next;
    Edge* e3 = (outerComponent->next)->next;

    Vertex_SPT* a = e1->getOrigin();
    Vertex_SPT* b = e2->getOrigin();
    Vertex_SPT* c = e3->getOrigin();


    const double ax = a->x;
    const double ay = a->y;
    const double bx = b->x;
    const double by = b->y;
    const double cx = c->x;
    const double cy = c->y;

    double A = 0.5 * ((bx - ax) * (cy - by) - (by - ay) * (cx - bx));

    if (A == 0) return false;
    // points are colinear or not distinct

    double xnum = ((cy - ay) * (by - ay) * (cy - by)) - ((bx * bx - ax * ax) * (cy - by)) + ((cx * cx - bx * bx) * (by - ay));
    double circum_x = xnum / (-4 * A);
    double circum_y = (-1 * (bx - ax) / (by - ay)) * (circum_x - 0.5 * (ax + bx)) + 0.5 * (ay + by);

    Vertex_SPT* circum = new Vertex_SPT(circum_x, circum_y, NULL);

    const double circum_radius = a->dist2(*circum);
    const double dist = v->dist2(*circum);

    return dist <= circum_radius;

}
bool Vertex_SPT::operator == (const Vertex_SPT& v) const {
    return ((fabs(x - v.x) < epsilon) && (fabs(y - v.y) < epsilon));
}

bool Vertex_SPT::operator > (const Vertex_SPT& v) const {
    if (fabs(y - v.y) < epsilon)
        return x < v.x;
    return y > v.y;
}

bool Vertex_SPT::operator < (const Vertex_SPT& v) const {
    if (fabs(y - v.y) < epsilon)
        return x > v.x;
    return y < v.y;
}

int min(int x, int y) {
    return x < y ? x : y;
}

bool Edge::isLeft(Vertex_SPT* v) const {
    Vertex_SPT* e_from = origin;
    Vertex_SPT* e_to = twin->origin;

    if (this->isBordering()) {
        return true;
    }
    double cross_prod = (e_to->x - e_from->x) * (v->y - e_from->y) - (e_to->y - e_from->y) * (v->x - e_from->x);
    return ((cross_prod + epsilon) > 0);
}

bool Face::inTriangle(Vertex_SPT* v) const {
    Edge* e1 = outerComponent;
    Edge* e2 = outerComponent->next;
    Edge* e3 = (outerComponent->next)->next;

    if (e1->isBordering() && e2->isBordering() && e3->isBordering())
        return true;

    return (e1->isLeft(v) && e2->isLeft(v) && e3->isLeft(v));
}


bool Edge::isBordering() const {
    Vertex_SPT* from = origin;
    Vertex_SPT* to = twin->origin;

    if ((from->id >= -2 && from->id <= 0) and (to->id >= -2 && to->id <= 0))
        return true;
    return false;
}

bool Edge::containsPoint(double x, double y) const {
    Vertex_SPT* a = origin;
    Vertex_SPT* b = twin->origin;

    if (this->isBordering())
        return false;

    if (x < std::min(a->x, b->x) || x > std::max(a->x, b->x))
        return false;

    if (y < std::min(a->y, b->y) || x > std::max(a->y, b->y))
        return false;

    if (fabs(b->x - a->x) < epsilon) {
        if (fabs(b->y - a->y) < epsilon)
            return true;
        else
            return false;
    }
    double x_rez = (x - a->x) / (b->x - a->x);
    double y_rez = (y - a->y) / (b->y - a->y);

    if (fabs(x_rez - y_rez) < epsilon)
        return true;
    return false;
}

void Edge::setEdgeParameters(Vertex_SPT* o, Edge* t, Face* f, Edge* n, Edge* p) {
    origin = o;
    twin = t;
    next = n;
    prev = p;
    incidentFace = f;
}


// glavni algoritam
void UnweightedShortestPathTree::generisiNasumicnaTemena(int brojTemena){

    _n  = brojTemena;
    const auto vertices = generisiNasumicneTacke(brojTemena);

    for (int i=0; i< brojTemena; i++)
        P.push_back(new Vertex_SPT(vertices[i].x(), vertices[i].y(), NULL));

    s = P[rand() % brojTemena];

}

void UnweightedShortestPathTree::ucitajPodatkeIzDatoteke(std::string imeDatoteke)
{
    /* Otvaranje zadatog fajla */
    std::ifstream datoteka(imeDatoteke);

    /* Broj temena */
    datoteka >> _n;

    int ind =0;

    /* Indeks pocetnog temena */
    datoteka >> ind;

    /* Ucitavanje svakog temena pojedinacno */
    for (auto i = 0; i < _n; i++) {
        double x, y;
        datoteka >> x >> y;
        P.push_back(new Vertex_SPT(x, y, NULL));
        P[i]->id = i;
    }

    s = P[ind];
}


void UnweightedShortestPathTree::inicijalizujDeloneGraf(){

    Face* beginTriangle = util.getSuperTriangle(P);
    TreeHistory* begining = new TreeHistory();
    begining->pointingFace = beginTriangle;
    beginTriangle->my_history = begining;

    d_object.root = begining;
    for (int i = 1; i < _n; i++)
    {
        Face* locateT = util.containingTriangle(d_object.root, P[i]);
        d_object.add_point_general(locateT, P[i]);
    }
    (*d_object.root).getDelaunayGraph(delaunayGraph);
}


void UnweightedShortestPathTree::inicijalizujUnitDiscGraf() {

    for (int i = 0; i < P.size(); i++)
        for (int j = 0; j < P.size(); j++)
        {
            double dist2 = sqrt(P[i]->dist2(*P[j]));
            if (dist2 <= 100 && dist2 > 0)
                unitDiscGraph[P[i]].insert(P[j]);
        }
}



/* Konstrukcija algoritma */
UnweightedShortestPathTree::UnweightedShortestPathTree(QWidget *pCrtanje,
                                         int pauzaKoraka,
                                         const bool &naivni,
                                         std::string imeDatoteke,
                                         int brojTemena)
    : AlgoritamBaza(pCrtanje, pauzaKoraka, naivni)
{

    /* Inicijalizacija niza temena i pocetnog temena */
    if (imeDatoteke == "")
        generisiNasumicnaTemena(brojTemena);
    else
        ucitajPodatkeIzDatoteke(imeDatoteke);

}

/* Deinicijalizacija algoritma */
UnweightedShortestPathTree::~UnweightedShortestPathTree() {
    // Clean up the vector of Vertex_SPT pointers
    for (auto Vertex_SPT : P) {
        delete Vertex_SPT;
    }
    P.clear();

    // Clean up the single Vertex_SPT pointer
    delete s;

    // Clean up the Delaunay graph map
    for (auto& pair : delaunayGraph) {
        // Iterate through the set of Vertex_SPT pointers and delete them
        for (auto Vertex_SPT : pair.second) {
            delete Vertex_SPT;
        }
        // Optionally, clear the set if needed
        pair.second.clear();
    }
    // Clear the map if needed
    delaunayGraph.clear();

    // Clear the Unit disc graph map
    // Clean up the Delaunay graph map
    for (auto& pair : unitDiscGraph) {
        // Iterate through the set of Vertex_SPT pointers and delete them
        for (auto Vertex_SPT : pair.second) {
            delete Vertex_SPT;
        }
        // Optionally, clear the set if needed
        pair.second.clear();
    }
    // Clear the map if needed
    unitDiscGraph.clear();
}


void UnweightedShortestPathTree::pokreniAlgoritam() {

    inicijalizujDeloneGraf();

    AlgoritamBaza_updateCanvasAndBlock()

    for (int i = 0; i < P.size(); i++)
        resultPathGlavni[P[i]] = make_pair(-1, nullptr);

    vector<std::set<Vertex_SPT*>> distanceSet = std::vector<std::set<Vertex_SPT*>>();

    std::set<Vertex_SPT*> W0;
    W0.insert(s);
    distanceSet.push_back(W0);
    resultPathGlavni.at(s).first = 0;

    int i = 1;

    while (!distanceSet[i - 1].empty()) {

        std::set<Vertex_SPT*> Q(distanceSet[i - 1]);

        std::set<Vertex_SPT*> Wi = std::set<Vertex_SPT*>();

        vector<Vertex_SPT*> Q_vec(Q.begin(), Q.end());

        KdNode* KdRoot = kdTree.buildKdTree(Q_vec, 0, Q_vec.size() - 1, 0);

        while (!Q.empty()) {

            // get first element in a set and erase it from the set
            Vertex_SPT* q = *Q.begin();
            Q.erase(Q.begin());


            std::set<Vertex_SPT*> delaunayAdj(delaunayGraph.at(q));

            for (auto itr1 = delaunayAdj.begin(); itr1 != delaunayAdj.end(); itr1++){
                {
                    Vertex_SPT* p = *itr1;

                    // distance is not set
                    if (resultPathGlavni.at(p).first == -1) {

                        // find closest neighbor w of p in W{i-1}
                        // use kd-trees
                        double minDist2 = numeric_limits<double>::max();
                        Vertex_SPT* w = kdTree.findNearestNeighbor(KdRoot, *p, minDist2);


                        if (sqrt(minDist2) <= 100) {

                            resultPathGlavni.at(p).first = i;
                            resultPathGlavni.at(p).second = new Vertex_SPT((*w).x, (*w).y, nullptr);
                            Q.insert(p);
                            Wi.insert(p);
                        }

                        AlgoritamBaza_updateCanvasAndBlock()
                    }


                }
            }
        }

        distanceSet.push_back(Wi);
        i++;
    }

    /* Obavestavanje pozivaoca o gotovoj animaciji */
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();
}


void UnweightedShortestPathTree::crtajAlgoritam(QPainter *painter) const {
    if (!painter) return;

    auto pen = painter->pen();
    pen.setColor(Qt::red);
    pen.setWidth(6);
    painter->setPen(pen);

    //draw source
    painter->drawPoint((*s).toQPointF());

    /* Okretanje cetkice kako brojevi i slova ne bi bili obrnuti */
    painter->save();
    painter->scale(1, -1);
    painter->translate(0, -2*(*s).y);
    painter->drawText(QPointF((*s).x - 5, (*s).y - 5), "S");
    painter->restore();

    //draw other points
    pen.setColor(Qt::blue);
    painter->setPen(pen);
    for (int i = 0; i < P.size(); i++)
        if(P[i]!=s)
            painter->drawPoint((*P[i]).toQPointF());


    // draw Delaunay Graph
    pen.setColor(Qt::green);
    pen.setWidth(2);
    painter->setPen(pen);
    for (auto itr1 = delaunayGraph.begin(); itr1 != delaunayGraph.end(); itr1++)
        for (auto itr2 = (itr1->second).begin(); itr2 != (itr1->second).end(); itr2++)
            painter->drawLine((itr1->first)->toQPointF(), (*itr2)->toQPointF());



    // draw resultPathGlavni
    pen.setColor(Qt::black);
    pen.setWidth(3);
    painter->setPen(pen);

    for (auto itr = resultPathGlavni.begin(); itr != resultPathGlavni.end(); itr++)
        // if path is calculated
        if(itr->second.first != -1 && itr->second.first != 0 ){

            Vertex_SPT v1 = *(itr->first);
            Vertex_SPT v2 = *itr->second.second;
            painter->save();
            painter->scale(1, -1);
            painter->translate(0, -2*v1.y);
            painter->drawText(v1.x - 5, v1.y - 5, QString::number(itr->second.first));
            painter->restore();
            painter->drawLine(v1.toQPointF(), v2.toQPointF());
        }
}

void UnweightedShortestPathTree::pokreniNaivniAlgoritam() {

    inicijalizujUnitDiscGraf();

    AlgoritamBaza_updateCanvasAndBlock()

    for (int i = 0; i < P.size(); i++)
        resultPathNaivni[P[i]] = make_pair(-1, nullptr);

    vector<std::set<Vertex_SPT*>> distanceSet = std::vector<std::set<Vertex_SPT*>>();
    std::set<Vertex_SPT*> W0;
    W0.insert(s);
    distanceSet.push_back(W0);
    resultPathNaivni.at(s).first = 0;

    int i = 1;

    while (!distanceSet[i - 1].empty()) {

        std::set<Vertex_SPT*> Q(distanceSet[i - 1]);

        std::set<Vertex_SPT*> Wi = std::set<Vertex_SPT*>();

        while (!Q.empty()) {

            // get first element in a set and erase it from the set
            Vertex_SPT* q = *Q.begin();
            Q.erase(Q.begin());


            std::set<Vertex_SPT*> unitDiscAdj(unitDiscGraph.at(q));

            for (auto itr1 = unitDiscAdj.begin(); itr1 != unitDiscAdj.end(); itr1++){
                {
                    Vertex_SPT* p = *itr1;

                    // distance is not set
                    if (resultPathNaivni.at(p).first == -1) {

                        resultPathNaivni.at(p).first = i;
                        resultPathNaivni.at(p).second = new Vertex_SPT((*q).x, (*q).y, nullptr);
                        Wi.insert(p);
                        AlgoritamBaza_updateCanvasAndBlock()

                    }
                }
            }
        }

        distanceSet.push_back(Wi);
        i++;
    }

    /* Obavestavanje pozivaoca o gotovoj animaciji */
    AlgoritamBaza_updateCanvasAndBlock()
    emit animacijaZavrsila();

}

void UnweightedShortestPathTree::crtajNaivniAlgoritam(QPainter *painter) const {

    if (!painter) return;

    auto pen = painter->pen();
    pen.setColor(Qt::red);
    pen.setWidth(6);
    painter->setPen(pen);

    //draw source
    painter->drawPoint((*s).toQPointF());

    /* Okretanje cetkice kako brojevi i slova ne bi bili obrnuti */
    painter->save();
    painter->scale(1, -1);
    painter->translate(0, -2*(*s).y);
    painter->drawText(QPointF((*s).x - 5, (*s).y - 5), "S");
    painter->restore();

    //draw other points
    pen.setColor(Qt::black);
    painter->setPen(pen);
    for (int i = 0; i < P.size(); i++)
        if(P[i]!=s)
            painter->drawPoint((*P[i]).toQPointF());

    // draw Unit Disc Graph
    pen.setColor(Qt::green);
    pen.setWidth(2);
    painter->setPen(pen);
    for (auto itr1 = unitDiscGraph.begin(); itr1 != unitDiscGraph.end(); itr1++)
        for (auto itr2 = (itr1->second).begin(); itr2 != (itr1->second).end(); itr2++)
            painter->drawLine((itr1->first)->toQPointF(), (*itr2)->toQPointF());


    // draw resultPathNaivni
    pen.setColor(Qt::black);
    pen.setWidth(3);
    painter->setPen(pen);

    for (auto itr = resultPathNaivni.begin(); itr != resultPathNaivni.end(); itr++)
        // if path is calculated
        if(itr->second.first != -1 && itr->second.first != 0 ){

            Vertex_SPT v1 = *(itr->first);
            Vertex_SPT v2 = *itr->second.second;
            painter->save();
            painter->scale(1, -1);
            painter->translate(0, -2*v1.y);
            painter->drawText(v1.x - 5, v1.y - 5, QString::number(itr->second.first));
            painter->restore();
            painter->drawLine(v1.toQPointF(), v2.toQPointF());
        }
}
