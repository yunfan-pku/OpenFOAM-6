#include <fstream>
#include <iomanip>
#include <iostream>

template <class T>
void csvwrite(std::string name, T &a, int flag = 0)
{
    const int sz = 10000;
    char mv[sz];
    std::ifstream ex(name);
    if (!ex || flag == 1)
    {
        std::ofstream csvfile(name);
        forAll(a, i)
        {
            csvfile << a[i] << std::endl;
        }
        csvfile.close();
    }
    else
    {
        std::ofstream tempcsvfile("__temp.csv");
        forAll(a, i)
        {
            ex.getline(mv, sz);
            tempcsvfile << mv << "," << a[i] << std::endl;
        }
        ex.close();
        tempcsvfile.close();
        std::string mv("mv ");
        std::string tempName("__temp.csv ");
        std::string com = mv + tempName + name;
        auto re = system(com.data());
    }

    return;
}

template <class T>
void csvwrite(std::string name, T &a,const char * title, int flag = 0)
{
    const int sz = 10000;
    char mv[sz];
    std::ifstream ex(name);
    if (!ex || flag == 1)
    {
        std::ofstream csvfile(name);
        csvfile << title << std::endl;
        forAll(a, i)
        {
            csvfile << a[i] << std::endl;
        }
        csvfile.close();
    }
    else
    {
        std::ofstream tempcsvfile("__temp.csv");
        ex.getline(mv, sz);
        tempcsvfile << mv << "," << title << std::endl;
        forAll(a, i)
        {
            ex.getline(mv, sz);
            tempcsvfile << mv << "," << a[i] << std::endl;
        }
        ex.close();
        tempcsvfile.close();
        std::string mv("mv ");
        std::string tempName("__temp.csv ");
        std::string com = mv + tempName + name;
        auto re = system(com.data());
    }

    return;
}
