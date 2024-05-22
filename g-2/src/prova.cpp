#include<iostream>
using namespace std;
struct prova{

    double x;
    double y;
};

class MyClass{

	public :
	double p5;

    double pici(){
                return 0;
                };
};


int main(){

struct prova p1 = {1,2};
struct prova *p2 = &p1;

cout << p2->x << endl;

MyClass a;
a.p5 = 5;

MyClass *myobject= &a;
cout << myobject -> pici() << endl;



return 0;

}