#include<fstream>
#include<iomanip>
#include<iostream>

template<class T>
void csvwrite(std::string name,T &a,int flag=0)
{
    const int sz=10000;
    char mv[sz];
    std::ifstream ex(name);
    if(!ex||flag==1)
    {
        std::ofstream csvfile(name);
        forAll(a, i)
        {
            csvfile<<a[i]<<std::endl;
        }
        csvfile.close();
    }
    else
    {
        std::ofstream  tempcsvfile("__temp.csv");
        forAll(a,i)
        {
            ex.getline(mv,sz);
            tempcsvfile<<mv<<","<<a[i]<<std::endl;
        }
        ex.close() ;
        tempcsvfile.close();
        std::string mv("mv ");
        std::string tempName("__temp.csv ");
        std::string com=mv+tempName+name;
        system(com.data());

    }

    return ;
}
int main()
{
    double a[100];
    for(int i=0; i<100; i++)
        a[i]=i;
    const int sz=10000;
    char mv[sz];
    csvwrite("test.csv",a);
}

