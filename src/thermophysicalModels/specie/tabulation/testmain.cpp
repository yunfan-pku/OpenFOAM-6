//#include"hybridtablulation.H"
//#include"hybridtablulation.C"
#include"delaunator.C"
#include"RTree.H"
#include"table.H"
#include"table.C"
#include"tabulation.cpp"
#include <chrono>
#include <iomanip>
#define TimerEnd    auto end = std::chrono::steady_clock::now();elapsed += end - start
#define TimerInit 	std::chrono::duration<double, std::micro> elapsed
#define TimerShow	std::cout<< "time: "  << elapsed.count()/1000000 << "s" << std::endl
#define TimerStart  auto start = std::chrono::steady_clock::now()

int main()
{
    RTree<int, double, 2,double> tree;
    double min[3],max[3],point[3]= {  715.599,  89178.3, 0.999}; // {101,2.1e6,0.0510};
    min[0]=1;
    max[0]=2;
    min[1]=1;
    max[1]=2;
    min[2]=1;
    max[2]=2;
    tree.Insert(min,max,100);
    int output=-1;
    int *pint=&output;

    std::cout<<tree.Search(point,pint)<<std::endl;
    std::cout<<output<<std::endl;
    TimerInit;
    TimerStart;
    tabulation ta("ThermoTable.txt");
    TimerEnd;
    TimerShow;
    double input[3]= {142.24600000,6009373.69200000,0.05090000},aoutput[3];//input argument
    //std::cout<<a.interpolation(input,1)<<std::endl;
    //std::cout<<a.interpolation(input,2)<<std::endl;

    // std::cout<<a.Search(point,pint)<<std::endl;
    std::cout<<output<<std::endl;
    typename RTree<coord_int, double, 3,double>::Point a,b,c,d,a2,b2,c2,d2,p;
    typename tabulation::Quad_2d quad_2d;
    typename tabulation::Hex hex;
    a.m[0]=1,a.m[1]=1,a.m[2]=0;
    b.m[0]=-1,b.m[1]=1,b.m[2]=0;
    c.m[0]=-1,c.m[1]=-1,c.m[2]=0;
    d.m[0]=1,d.m[1]=-1,d.m[2]=0;
    a2.m[0]=1,a2.m[1]=0,a2.m[2]=1;
    b2.m[0]=0,b2.m[1]=1,b2.m[2]=1;
    c2.m[0]=-1,c2.m[1]=0,c2.m[2]=1;
    d2.m[0]=0,d2.m[1]=-1,d2.m[2]=1;
    p.m[0]=1,p.m[1]=0.5,p.m[2]=0.4;
    quad_2d.vertex[0]=a;
    quad_2d.vertex[1]=b;
    quad_2d.vertex[2]=c;
    quad_2d.vertex[3]=d;
    hex.vertex[0]=a;
    hex.vertex[1]=b;
    hex.vertex[2]=c;
    hex.vertex[3]=d;
    hex.vertex[4]=a2;
    hex.vertex[5]=b2;
    hex.vertex[6]=c2;
    hex.vertex[7]=d2;
    aoutput[2]=input[2];
    std::cout<<(aoutput[0]=ta.interpolation(input,1))<<" "<<(aoutput[1]=ta.interpolation(input,2))<<std::endl;
    std::cout<<tabulation::point_in_tab_hex(p,  hex)<<std::endl;
    c=1.7*a-2.6*b;
    ta.Search_(aoutput);
    //std::cout<<ta.indv[0][ta.ret_coord.m[0]]<<" "<<ta.indv[1][ta.ret_coord.m[1]]<<" "<<ta.indv[2][ta.ret_coord.m[2]]<<std::endl;

	std::ofstream fout("ctest.txt");
    for(double i=100; i<=1597; i+=21.123)
        for(double j=2e6; j<=1.18e7; j+=100234.3423)
            for(double k=0.0509; k<=0.999; k+=0.03435)
            {
                input[0]=i;
                input[1]=j;
                input[2]=k;
                aoutput[0]=ta.interpolation(input,1);
                aoutput[1]=ta.interpolation(input,2);
                aoutput[2]=input[2];
                ta.interpolation_tab(aoutput[0],aoutput[1],aoutput[2],1);
                fout<<std::setiosflags(std::ios::fixed)<<std::setprecision(8)<<std::fixed<<input[0]<<" "<<input[1]<<" "<<input[2]<<" "<<ta.interpolation_tab(aoutput[0],aoutput[1],aoutput[2],1)<<" "<<ta.interpolation_tab(aoutput[0],aoutput[1],aoutput[2],2)<<std::endl;
            }



    ///**/ta.interpolation_tab(aoutput[0],aoutput[1],aoutput[2],1);
    fout.close();
    return 0;

}
