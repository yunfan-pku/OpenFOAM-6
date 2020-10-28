#include "ISATleaf.H"
#include "ISATNode.H"
#include "ISATmanager.H"
#include "ISATbinaryTree.H"
#include "dictionary.H"
#include<fstream>
using namespace Foam;

struct function
{
    void  value(const scalarList& in, scalarList& out)
    {
        out.resize(1);
        out[0] = in[0]*in[0]+ in[1]* in[1]+20;//+2000*::Foam::sin(in[1]/10);
    }
    void  derive(const scalarList& in, scalarRectangularMatrix& out)
    {
        scalarList temp = in;
        scalarList tempv;
        value(in, tempv);
        scalarList tempv2;
        scalarList dd;
        scalar dx = 1e-5;
        for (int i = 0;i < in.size();i++)
        {
            temp[i] += dx;
            value(temp, tempv2);
            dd = (tempv2 - tempv) / dx;
            for (int j = 0;j < tempv2.size();j++)
            {
                out[j][i] = dd[j];
            }
            temp[i] -= dx;
        }
    }
} funcxy;
int main()
{
    //ISATleaf A(1), B(2), * pl;
    ISATleaf* pl;
    //ISATNode C(&A, &B);
    ISATmanager<function> a(2, 1, funcxy);
    ISATbinaryTree& T = a.tableTree();
    a.epsilon() = 1e-2;
    a.relepsilon() = 1e-2;
    scalarList L(2);
    scalarList R(1);
    ISATleaf* pleaf;
    /*
        for (int i = 0;i < 1;i++)
        {
            L[0] = rand() % 1000 / 10.0;
            L[1] = rand() % 1000 / 10.0;
            //Info<<L<<endl;
            a.add(L);
            //funcxy.value(L,R);
            //R[0] = L[0] * L[1];
            //pleaf=T.insertNewLeaf(L, R);
            //funcxy.derive(L,pleaf->A());
        }
        */
        /*
        pl=a.search(L);
        bool bt;
        scalarList ret;
        //std::ofstream fout("out.csv");
        L[0] = 20;
        L[1] = 20;
        a.call(L,ret);
        L[0] = 20.1;
        L[1] = 20.1;
        a.call(L,ret);
        */

    string sout("out"), scsv(".csv");
    std::ofstream fout;
    scalarList ret;
    /*
    for (double i = 0;i <= 400;i+=0.1)
    {
        for (double j = 20;j <= 40;j+=0.1)
        {
            L[0] = i;
            L[1] = j;
            a.call(L,ret);
            fout << i << "," << j << "," << ret[0] << std::endl;
        }
    } */
    for (int i = 0;i <= 40000;i++)
    {
        L[0] = rand() % 1000 / 100.0 - 5;
        L[1] = rand() % 1000 / 100.0 - 5;
        a.call(L, ret);
        Info << a.treesize() << endl;
        if (i % 1 == -1)
        {
            fout.open(sout + std::to_string(i) + scsv);
            for (double j = -5;j <= 5;j += 0.5)
            {
                for (double k = -5;k <= 5;k += 0.5)
                {
                    //a.tablevalue(L, ret);
                    L[0] = j;
                    L[1] = k;
                    a.tablevalue(L, ret);
                    fout << L[0] << "," << L[1] << "," << ret[0] << std::endl;
                }
            }
        }
        fout.close();
    }


    Info << a.treesize() << endl;

    /*
    for (int i = 0;i <= 100;i += 10)
    {
        for (int j = 0;j <= 100;j += 10)
        {
            L[0] = i;
            L[1] = j;
            T.insertNewLeaf(L, L[0] * L[1]);
        }
    }
*/
//T.insertNewLeaf(10);
//T.insertNewLeaf(12);
//T.insertNewLeaf(11);

//Info<<A.value()<<endl;
/*
    T.print(Info);
    Info << endl;
    Info << T << endl;
    Info << T.size() << " " << T.depth() << endl;
    L[0] = 27.5;
    L[1] = 10.5;
    T.binaryTreeSearch(L, T.root(), pl);
    Info << (*pl).value() << endl;
    //T.clear();
    Info << T.size() << endl;
    scalarList ret;
    T.eval(L, ret);
    Info << ret[0] << endl;
    std::ofstream fout("out.csv");
    for (int i = 50;i <= 100;i++)
    {
        for (int j = 50;j <= 100;j++)
        {
            L[0] = i;
            L[1] = j;
            T.eval(L, ret);
            fout << i << "," << j << "," << ret[0] << std::endl;
        }
    }
    */
}