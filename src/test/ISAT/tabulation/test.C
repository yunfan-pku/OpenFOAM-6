#include "ISATleaf.H"
#include "ISATNode.H"
#include "ISATbinaryTree.H"
#include "dictionary.H"
using namespace Foam;
int main()
{
    //ISATleaf A(1), B(2), * pl;
    ISATleaf * pl;
    //ISATNode C(&A, &B);
    ISATbinaryTree T;
    scalarList L(2);
    for (int i=0;i<=5000;i++)
    {
        L[0]=rand() % 100;
        L[1]=rand() % 100;
        T.insertNewLeaf(L);

    }

    //T.insertNewLeaf(10);
    //T.insertNewLeaf(12);
    //T.insertNewLeaf(11);

    //Info<<A.value()<<endl;
    T.print(Info);
    Info<<endl;
    Info<<T<<endl;
    Info<<T.size()<<" "<<T.depth()<<endl;
    L[0]=27.5;
    L[1]=10.5;
    T.binaryTreeSearch(L, T.root(), pl);
    Info<<(*pl).value()<<endl;
    T.clear();
    Info<<T.size()<<endl;

}