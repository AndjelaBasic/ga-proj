#include "ga07_delaunay.h"
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>

#define epsilon 0.0001
DCEL_TRIANGLE d_object;
using namespace std;

Face* CopyOfTriangle(Face* triangle){
    Edge* e1_old = triangle->outerComponent;
    Edge* e2_old = (triangle->outerComponent)->next;
    Edge* e3_old = ((triangle->outerComponent)->next)->next;

    Vertex* a = e1_old->getOrigin();
    Vertex* b = e2_old->getOrigin();
    Vertex* c = e3_old->getOrigin();

    Edge* e1 = new Edge(); Edge* e1_twin = new Edge();
    Edge* e2 = new Edge(); Edge* e2_twin = new Edge();
    Edge* e3 = new Edge(); Edge* e3_twin = new Edge();

    Vertex* v1_new = new Vertex(a->x, a->y, e1);
    Vertex* v2_new = new Vertex(b->x, b->y, e2);
    Vertex* v3_new = new Vertex(c->x, c->y, e3);
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

void printLeaves(TreeHistory* root){
    if ((root->childrenFaces).size()==0)
        (root->pointingFace)->printFace();
    if ((root->childrenFaces).size()==2){
        printLeaves((root->childrenFaces)[0]);
        printLeaves((root->childrenFaces)[1]);
    }
    if ((root->childrenFaces).size()==3){
        printLeaves((root->childrenFaces)[0]);
        printLeaves((root->childrenFaces)[1]);
        printLeaves((root->childrenFaces)[2]);
    }
}

void printHistory(TreeHistory* root){
    if ((root->childrenFaces).size()==0){
        cout<<"LEAF"<<endl;
        (root->pointingFace)->printFace();
    }
    if ((root->childrenFaces).size()==2){
        (root->pointingFace)->printFace();
        cout<<"2 of them"<<endl;
        printHistory((root->childrenFaces)[0]);
        printHistory((root->childrenFaces)[1]);
        cout<<"END OF TWO"<<endl;
    }
    if ((root->childrenFaces).size()==3){
        cout<<"3 of them"<<endl;
        (root->pointingFace)->printFace();
        printHistory((root->childrenFaces)[0]);
        printHistory((root->childrenFaces)[1]);
        printHistory((root->childrenFaces)[2]);
        cout<<"END OF TWO"<<endl;
    }
}


// Function to return the next random number
int getNum(vector<int>& v)
{

    // Size of the vector
    int n = v.size();

    // Generate a random number
    srand(time(NULL));

    // Make sure the number is within
    // the index range
    int index = rand() % n;

    // Get random number from the vector
    int num = v[index];

    // Remove the number from the vector
    swap(v[index], v[n - 1]);
    v.pop_back();

    // Return the removed number
    return num;
}

// Function to generate n non-repeating random numbers
std::vector<int> generateRandom(int n)
{
    vector<int> v(n);
    vector<int> result_id;
    // Fill the vector with the values
    // 1, 2, 3, ..., n
    for (int i = 0; i < n; i++)
        v[i] = i + 1;

    // While vector has elements
    // get a random number from the vector and print it
    while (v.size()) {
        //cout << getNum(v) << " ";
        result_id.push_back(getNum(v));
    }
    return result_id;
}

void legalize_edge(Vertex* v, Edge* e){

    if (e->isBordering())
        return;

    Vertex* p = (e->prev)->origin;

    Edge* considered_edge = e;
    cout<<"TEST OF EQUALITY"<<endl;
    cout<<((*p)==(*v))<<endl;
    p->printVertex();
    v->printVertex();
    e->printEdge();
    (e->twin)->printEdge();
    (((e->twin)->prev)->origin)->printVertex();
    if ((*p)==(*v)){
        p = ((e->twin)->prev)->origin;
        considered_edge = e->twin;
    }

    cout<<"LEGALIZE VERTICES"<<endl;
    cout<<"p="<<endl;
    p->printVertex();
    cout<<"v="<<endl;
    v->printVertex();
    cout<<"END OF LEGALIZE VERTICES"<<endl;
    (considered_edge->incidentFace)->printFace();
    ((considered_edge->twin)->incidentFace)->printFace();
    cout<<"WATCH"<<endl;
    Face* my_triangle = considered_edge->incidentFace;
    Edge* next_to_legalize1 = considered_edge->next;
    Edge* next_to_legalize2 = considered_edge->prev;

    int id1 = v->id;
    int id2 = p->id;
    int id3 = (e->origin)->id;
    int id4 = ((e->twin)->origin)->id;
    if (id1>=0 && id2>=0 && id3>=0 && id4>=0){
        if (my_triangle->circumCircleContains(p)){
            considered_edge->flip();
            legalize_edge(v,next_to_legalize1);
            legalize_edge(v,next_to_legalize2);
        }
    } else if (min(id1,id2)>min(id3,id4)) {
        considered_edge->flip();
        legalize_edge(v,next_to_legalize1);
        legalize_edge(v,next_to_legalize2);
    }

}


void DCEL_TRIANGLE::add_point_general(Face* t, Vertex* v){

    Face* triangle = t;
    Edge* e1 = triangle->outerComponent;
    Edge* e2 = (triangle->outerComponent)->next;
    Edge* e3 = ((triangle->outerComponent)->next)->next;

    double x = v->x;
    double y = v->y;

    if (e1->isBordering() && e2->isBordering() && e3->isBordering()){
        add_point_in_triangle(triangle, v);
        cout<<"add point triangle"<<endl;
    }
    else if (e1->containsPoint(x,y)){
        cout<<"e1 has it"<<endl;
        add_point_on_edge(triangle, e1, v);
    } else if (e2->containsPoint(x,y)){
        cout<<"e2 has it"<<endl;
        add_point_on_edge(triangle, e2, v);
    } else if (e3->containsPoint(x,y)){
        cout<<"e3 has it"<<endl;
        add_point_on_edge(triangle, e3, v);
    } else{
        cout<<"add point triangle"<<endl;
        add_point_in_triangle(triangle, v);
        /*legalize_edge(v, e1);
        legalize_edge(v, e2);
        legalize_edge(v, e3);*/
    }

}
void DCEL_TRIANGLE::add_point_in_triangle(Face* triangle, Vertex* v){
    cout<<"ADDING IN TRIANGLE"<<endl;
    triangle->printFace();

    Face* t = CopyOfTriangle(triangle);
    t->printFace();

    TreeHistory* history = triangle->my_history;

    Edge* e1 = triangle->outerComponent;
    Edge* e2 = (triangle->outerComponent)->next;
    Edge* e3 = ((triangle->outerComponent)->next)->next;

    Vertex* a = e1->getOrigin();
    Vertex* b = e2->getOrigin();
    Vertex* c = e3->getOrigin();


    Edge* a1 = new Edge(); Edge* a2 = new Edge();
    Edge* b1 = new Edge(); Edge* b2 = new Edge();
    Edge* c1 = new Edge(); Edge* c2 = new Edge();

    v->incidentEdge = a2;

    Face* f1 = new Face(a2);
    Face* f2 = new Face(b2);
    Face* f3 = new Face(c2);

    e1->next = b1; e1->prev=a2; e1->incidentFace=f1;
    e2->next = c1; e2->prev=b2; e2->incidentFace=f2;
    e3->next = a1; e3->prev=c2; e3->incidentFace=f3;

    a1->setEdgeParameters(a,a2,f3,c2,e3);
    a2->setEdgeParameters(v,a1,f1,e1,b1);
    b1->setEdgeParameters(b,b2,f1,a2,e1);
    b2->setEdgeParameters(v,b1,f2,e2,c1);
    c1->setEdgeParameters(c,c2,f2,b2,e2);
    c2->setEdgeParameters(v,c1,f3,e3,a1);


    history->pointingFace = t;
    (history->pointingFace)->my_history = history;
    cout  <<"***********************************************************************"<<endl;
    cout  <<"***********************************************************************"<<endl;
    cout<<"check"<<endl;
    (history->pointingFace)->printFace();
    cout <<"done here"<<endl;

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

    cout<<"check"<<endl;
    (history->pointingFace)->printFace();
    cout <<"done here"<<endl;

    legalize_edge(v, e1);
    legalize_edge(v, e2);
    legalize_edge(v, e3);

    cout<<"check"<<endl;
    (history->pointingFace)->printFace();
    cout <<"done here"<<endl;
    cout  <<"***********************************************************************"<<endl;
    cout  <<"***********************************************************************"<<endl;

}

void DCEL_TRIANGLE::add_point_on_edge(Face* triangle, Edge* e, Vertex* v){
    cout<<"EDGE ADDING"<<endl;
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

    e1->setEdgeParameters(e->origin, e1_twin, f1, e3, e_prev);
    e1_twin->setEdgeParameters(v, e1, f2, e_twin->next, e4);
    e2->setEdgeParameters(v, e2_twin, f3, e_next, e3_twin);
    e2_twin->setEdgeParameters(e_twin->origin, e2, f4, e4_twin, e_twin->prev);
    e3->setEdgeParameters(v, e3_twin, f1, e_prev, e1);
    e3_twin->setEdgeParameters(e_prev->origin, e3, f3, e2, 	e_next);
    e4->setEdgeParameters((e_twin->prev)->origin, e4_twin, f2, e1_twin,e_twin->next);
    e4_twin->setEdgeParameters(v, e4_twin, f4, e_twin->prev, e2_twin);


    history->pointingFace = CopyOfTriangle(triangle);
    history_incident->pointingFace = CopyOfTriangle(triangle_incident);
    (history->pointingFace)->my_history = history;
    (history_incident->pointingFace)->my_history = history_incident;

    //Face* f = new Face(*f1);
    TreeHistory* newFace1 = new TreeHistory();
    newFace1->pointingFace = f1;
    f1->my_history = newFace1;
    (history->childrenFaces).push_back(newFace1);

    //f = new Face(*f3);
    TreeHistory* newFace3 = new TreeHistory();
    newFace1->pointingFace = f3;
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

    /*
    (triangle->childrenFaces).push_back(f1);
    (triangle->childrenFaces).push_back(f3);
    ((e_twin->incidentFace)->childrenFaces).push_back(f2);
    ((e_twin->incidentFace)->childrenFaces).push_back(f4);*/

    /*f1->printFace();
    f2->printFace();
    f3->printFace();
    f4->printFace();*/

    legalize_edge(v,e_next);
    legalize_edge(v,e_prev);
    legalize_edge(v,e_twin->next);
    legalize_edge(v,e_twin->prev);

}
void Vertex::printVertex() const{
    cout << "vertex"<<endl;
    cout << "x= "<< x <<" "<<"y= "<<y <<endl;
}

void Edge::flip(){
    //Edge* a1 = new Edge();
    //a1->setEdgeParameters(origin, twin, incidentFace,next,prev);

    Edge* a2 = this->twin;

    Edge* e1 = this->next;
    Edge* e2 = this->prev;

    Edge* e3 = a2->next;
    Edge* e4 = a2->prev;

    Vertex* v1 = e1->origin;
    Vertex* v2 = e4->origin;
    Vertex* v3 = this->origin;
    Vertex* v4 = e2->origin;

    //Face* const AA =  a1->incidentFace;
    //Face* const BB =  a2->incidentFace;
    TreeHistory* history = (this->incidentFace)->my_history;
    TreeHistory* history_incident = (twin->incidentFace)->my_history;

    history->pointingFace = CopyOfTriangle(this->incidentFace);
    history_incident->pointingFace = CopyOfTriangle(a2->incidentFace);
    (history->pointingFace)->my_history = history;
    (history_incident->pointingFace)->my_history=history_incident;

    Face* A = this->incidentFace;
    Face* B = a2->incidentFace;
    cout<<"OLD"<<endl;
    A->printFace();
    B->printFace();
    this->setEdgeParameters(v2,a2,A,e2,e3);
    a2->setEdgeParameters(v4,this,B,e4,e1);

    v1->incidentEdge = e1;
    v3->incidentEdge = e3;

    A->outerComponent = e2;
    B->outerComponent = e1;

    e1->incidentFace = B;
    e4->incidentFace = B;
    e3->incidentFace = A;
    e4->incidentFace = A;

    help_triangle(this,e2,e3);
    help_triangle(a2,e4,e1);

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

    cout<<"HELLO FLIP"<<endl;
    cout<<"Fliping egde"<<endl;
    this->printEdge();
    A->printFace();
    B->printFace();
    cout<<"END FLIP"<<endl;
}

void Edge::help_triangle(Edge* a, Edge* b, Edge* c){
    a->next=b; b->next=c; c->next =a;
    c->prev=b; b->prev=a; a->prev =c;
}

double Vertex::dist2(const Vertex& v) const{
    return (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y);
}
double Vertex::norm2() const{
    return x*x + y*y;
}

bool Face::circumCircleContains(Vertex* v){
    Edge* e1 = outerComponent;
    Edge* e2 = outerComponent->next;
    Edge* e3 = (outerComponent->next)->next;

    Vertex* a = e1->getOrigin();
    Vertex* b = e2->getOrigin();
    Vertex* c = e3->getOrigin();

    const double ab = a->norm2();
    const double cd = b->norm2();
    const double ef = c->norm2();

    const double ax = a->x;
    const double ay = a->y;
    const double bx = b->x;
    const double by = b->y;
    const double cx = c->x;
    const double cy = c->y;

    const double circum_x = (ab*(cy-by)+cd*(ay-cy)+ef*(by-ay))/(ax*(cy-by)+bx*(ay-cy)+cx*(by-ay));
    const double circum_y = (ab*(cx-bx)+cd*(ax-cx)+ef*(bx-ax))/(ay*(cx-bx)+by*(ax-cx)+cy*(bx-ax));
    Vertex* circum = new Vertex(circum_x, circum_y, NULL);
    const double circum_radius = a->dist2(*circum);
    const double dist = v->dist2(*circum);
    return dist <= circum_radius;

}
bool Vertex::operator == (const Vertex& v) const{
    return ((fabs(x-v.x)<epsilon) && (fabs(y-v.y)<epsilon));
}

bool Vertex::operator > (const Vertex& v) const{
    if (fabs(y-v.y)<epsilon)
        return x < v.x;
    return y > v.y;
}

bool Edge::isLeft(Vertex* v) const{
    Vertex* e_from = origin;
    Vertex* e_to = twin->origin;

    if (this->isBordering()){
        cout<<"yes is left from bordering"<<endl;
        return true;
    }

    if (e_from->id == -1){
        return	v < e_to;
    } else if (e_to->id==-1){
        return  v > e_from;
    } else if (e_to->id==-2){
        cout << "mark"<<(v < e_from)<<endl;
        return v < e_from;
    } else if (e_from->id==-2) {
        return v > e_to;
    } else {
        double cross_prod=(e_to->x - e_from->x)*(v->y - e_from->y)-(e_to->y - e_from->y)*(v->x-e_from->x);
        cout<<"contains"<<((cross_prod+epsilon)>0)<<endl;
        cout<<cross_prod<<endl;
        return ((cross_prod+epsilon) > 0);
    }
}

bool Face::inTriangle(Vertex* v) const{
    Edge* e1 = outerComponent;
    Edge* e2 = outerComponent->next;
    Edge* e3 = (outerComponent->next)->next;

    if (e1->isBordering() && e2->isBordering() && e3->isBordering())
        return true;

    cout<<"testing vertex"<<endl;
    v->printVertex();
    cout<<"e1 is left"<<e1->isLeft(v)<<endl;
    e1->printEdge();
    cout<<"e2 is left"<<e2->isLeft(v)<<endl;
    e2->printEdge();
    cout<<"e3 is left"<<e3->isLeft(v)<<endl;
    e3->printEdge();

    return (e1->isLeft(v) && e2->isLeft(v) && e3->isLeft(v));
}


bool Edge::isBordering() const{
    Vertex* from = origin;
    Vertex* to = twin->origin;

    if ((from->id >=-2 && from->id<=0) and (to->id >=-2 && to->id<=0))
        return true;
    return false;
}

void Edge::printEdge() const{
    cout << "edge"<<endl;
    cout << "from ";
    origin->printVertex();
    cout << "to ";
    twin->origin->printVertex();
    cout << endl;
}

bool Edge::containsPoint(double x, double y) const{
    Vertex* a = origin;
    Vertex* b = twin->origin;

    if (this->isBordering())
        return false;

    if (fabs(b->x-a->x)<epsilon){
        if (fabs(b->y-a->y)<epsilon)
            return true;
        else
            return false;
    }
    double x_rez =(x-a->x)/(b->x-a->x);
    double y_rez =(y-a->y)/(b->y-a->y);

    if (fabs(x_rez-y_rez)<epsilon)
        return true;
    return false;
}

bool testTriangle(Face* f, Vertex* v){
    Edge* e1 = f->outerComponent;
    Edge* e2 = (f->outerComponent)->next;
    Edge* e3 = ((f->outerComponent)->next)->next;
    return (f->inTriangle(v) || e1->containsPoint(v->x, v->y) || e2->containsPoint(v->x, v->y) || e3->containsPoint(v->x, v->y));
}

int min(int x, int y){
    return x<y ? x:y;
}

void Face::printFace() const{
    cout << "triangle";
    cout << "---------------------------------------------"<<endl;
    outerComponent->printEdge();
    (outerComponent->next)->printEdge();
    ((outerComponent->next)->next)->printEdge();
    cout << "---------------------------------------------"<<endl;
    cout << endl<<endl<<endl;
}

void Edge::setEdgeParameters(Vertex* o, Edge* t, Face* f, Edge* n, Edge* p){
    origin=o;
    twin=t;
    next=n;
    prev=p;
    incidentFace=f;
}

Face* containingTriangle(TreeHistory* f, Vertex* v){
    if ((f->childrenFaces).size() == 0){
        return f->pointingFace;
    } else if ((f->childrenFaces).size() == 2){
        if (testTriangle((f->childrenFaces)[0]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[0],v);
        else
            return containingTriangle((f->childrenFaces)[1],v);
    } else{
        if (testTriangle((f->childrenFaces)[0]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[0],v);
        else if (testTriangle((f->childrenFaces)[1]->pointingFace, v))
            return containingTriangle((f->childrenFaces)[1],v);
        else
            return containingTriangle((f->childrenFaces)[2],v);
    }
}
int main(){
    /*
    DCEL_TRIANGLE d_object;

    Edge* e1 = new Edge();
    Edge* e2 = new Edge();
    Edge* e3 = new Edge();
    Edge* e1_twin = new Edge();
    Edge* e2_twin = new Edge();
    Edge* e3_twin = new Edge();


    Vertex* v1 = new Vertex(1.0,2.0, e1);
    Vertex* v2 = new Vertex(-2.0, 3.4, e2);
    Vertex* v3 = new Vertex(10.0, 8.0, e3);


    Vertex* v = new Vertex(2.0, 4.0, NULL);
    Face* d = new Face(e1);


    e1->setEdgeParameters(v1,e1_twin, d, e3,e2);
    e2->setEdgeParameters(v2,e2_twin, d, e1, e3);
    e3->setEdgeParameters(v3,e3_twin, d, e2, e1);
    e1_twin->setEdgeParameters(v3, e1, NULL, e2_twin, e3_twin);
    e2_twin->setEdgeParameters(v1, e2, NULL, e3_twin, e1_twin);
    e3_twin->setEdgeParameters(v2, e3, NULL, e1_twin, e2_twin);
    d_object.root = d;
    d_object.add_point_general(d, v);*/
    int n;
    double x,y;
    cin>>n;
    n = n-1;
    std::vector<Vertex*> help_points(n+1);
    for (int i=0;i<n+1;i++)
        help_points[i]=NULL;
    int poz_max=0;
    for (int i=0;i<n+1;i++){
        cin >>x>>y;
        Vertex* newVertex = new Vertex(x,y, NULL);
        help_points[i] = newVertex;
        //newVertex->printVertex();
        //help_points[i]->printVertex();
        if ((*newVertex) > (*help_points[poz_max]))
            poz_max = i;
    }
    Vertex* tmp = help_points[0];
    help_points[0] = help_points[poz_max];
    help_points[poz_max] = tmp;
    //help_points[n]->printVertex();
    //help_points[poz_max]->printVertex();
    for (int i=0;i<n+1;i++){
        help_points[i]->printVertex();
    }


    std::vector<Vertex*> in_points(n+1);
    for (int i=0;i<n+1;i++)
        in_points[i]=NULL;
    in_points[0] = help_points[0];
    in_points[0]->id = 0;

    std::vector<int> rez=generateRandom(n);

    for (int i=1;i<n+1;i++){

        help_points[i]->id = rez[i-1];
        in_points[rez[i-1]] = help_points[i];
        help_points[i]->printVertex();
        cout<<help_points[i]->id<<endl;
    }

    for (int i=0;i<n+1;i++){
        in_points[i]->printVertex();
        cout << in_points[i]->id <<endl;
    }

    Vertex* pminus1 = new Vertex(-1,-1,NULL);
    pminus1->id = -1;
    Vertex* pminus2 = new Vertex(-2,-2,NULL);
    pminus2->id = -2;


    Edge* e1 = new Edge();
    Edge* e2 = new Edge();
    Edge* e3 = new Edge();
    Edge* e1_twin = new Edge();
    Edge* e2_twin = new Edge();
    Edge* e3_twin = new Edge();

    Face* beginTriangle = new Face(e1);


    e1->setEdgeParameters(in_points[0],e1_twin, beginTriangle, e2,e3);
    e2->setEdgeParameters(pminus2,e2_twin, beginTriangle, e3, e1);
    e3->setEdgeParameters(pminus1,e3_twin, beginTriangle, e1, e2);
    e1_twin->setEdgeParameters(pminus2, e1, NULL, e3_twin, e2_twin);
    e2_twin->setEdgeParameters(pminus1, e2, NULL, e1_twin, e3_twin);
    e3_twin->setEdgeParameters(in_points[0], e3, NULL, e2_twin, e1_twin);

    TreeHistory* begining = new TreeHistory();
    begining->pointingFace = beginTriangle;
    beginTriangle->my_history = begining;

    d_object.root = begining;
    //beginTriangle->printFace();
    /*for (int i=1;i<n+1;i++)
        cout << beginTriangle->inTriangle(in_points[i])<<endl;

    Face* locateT = d_object.containingTriangle(d_object.root, in_points[1]);
    locateT->printFace();*/
    //printHistory(d_object.root);
    for (int i=1; i<n+1; i++)
    {
        if (i<=1){
            Face* locateT = containingTriangle(d_object.root, in_points[i]);
            cout<<"In T"<<endl;
            locateT->printFace();
            cout<<"HERE"<<endl;
            d_object.add_point_general(locateT, in_points[i]);
        }
        if (i==2){
            cout<<"found"<<endl;
            Face* locateT = containingTriangle(d_object.root, in_points[i]);
            locateT->printFace();
            d_object.add_point_general(locateT, in_points[i]);

        }
    }
    cout<<"HELLO HISTORY"<<endl;
    printHistory(d_object.root);
    cout<<"hello root"<<endl;
    ((d_object.root)->pointingFace)->printFace();
    return 0;
}
