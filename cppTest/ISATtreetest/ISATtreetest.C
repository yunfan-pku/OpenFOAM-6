#include "ISATmanager.H"
#include "scalarList.H"
#include<iostream>
struct myfunc
{
    void value(const Foam::scalarList& in, Foam::scalarList& out) const
    {
        out.resize(1);
        out[0] = in[0] * in[0] - in[1] * in[1];
    }

    void  derive(const Foam::scalarList& in, Foam::scalarRectangularMatrix& out)
    {
        out[0][0] = 2 * in[0];
        out[1][1] = -2 * in[1];
        out[1][0] = 0;
        out[0][1] = 0;
    }
};
using namespace Foam;


int main()
{
    myfunc f;
    Foam::ISATmanager<myfunc> ISATtest(2, 1, f, "test");
    Foam::scalarList in(2), out(1);
    ISATtest.epsilon() = 1e-2;
    ISATtest.scaleFactor()[0][0] = 1;
    ISATtest.relepsilon() = 0;// 5e-2;
    ISATtest.init_elp()[0][0] = 1 / 0.001;
    ISATtest.init_elp()[1][1] = 1 / 0.001;
    ISATtest.tableTree().setmaxNLeafs(60);
    in[0] = 0.1;
    in[1] = 0.1;
    ISATtest.call(in, out);
    in[0] = 0.2;
    in[1] = 0.2;
    ISATtest.call(in, out);
    in[0] = 0.11;
    in[1] = 0.11;
    ISATtest.call(in, out);
    //auto p = ISATtest.search(in);
    in[0] = 0.21;
    in[1] = 0.21;
    ISATtest.call(in, out);
    in[0] = 0.11;
    in[1] = 0.11;
    ISATtest.call(in, out);
    //std::cout << p->numRetrieve() << std::endl;
    //ISATtest.tableTree().deleteLeaf(p);
    ISATtest.call(in, out);


    for (int i = 0;i < 20;i++)
    {
        ISATtest.newTimeStep();
        for (int j = 0;j < 10;j++)
        {
            in[0] = (rand() % 100) / 200.0;
            in[1] = (rand() % 100) / 200.0;
            ISATtest.call(in, out);
            Info << in[0] * in[0] - in[1] * in[1] << "," << in[0] * in[0] - in[1] * in[1] - out[0] << endl;
        }
    }
    ISATtest.newTimeStep();
    for (int j = 0;j < 10;j++)
    {
        in[0] = (rand() % 100) / 200.0;
        in[1] = (rand() % 100) / 200.0;
        ISATtest.call(in, out);
        Info << in[0] * in[0] - in[1] * in[1] << "," << in[0] * in[0] - in[1] * in[1] - out[0] << endl;
    }
    std::cout << ISATtest.treesize() << std::endl;
    ISATtest.tableTree().timeTagList().print();
    return 0;
}