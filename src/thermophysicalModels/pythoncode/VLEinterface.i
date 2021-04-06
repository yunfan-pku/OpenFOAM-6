%module VLE

%include <std_vector.i>
%include <std_string.i>
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;
%naturalvar solver::X;
%naturalvar solver::specie;
%naturalvar solver::equalconstant;
%naturalvar solver::comp_liq;
%naturalvar solver::comp_gas;
%typemap(out) std::vector<double> %{
    $result = PyList_New($1.size());
    for (int i = 0; i < $1.size(); ++i)
        PyList_SET_ITEM($result,i,PyFloat_FromDouble($1[i]));
%}

%naturalvar solver_new::comp;
%naturalvar solver_new::equalconstant;
%naturalvar solver_new::comp_liq;
%naturalvar solver_new::comp_gas;
%naturalvar solver_new::ret;
%naturalvar solver_new::m_specie;
%{
  #include"VLEinterface.h"
  /*
struct solver
{
    std::string path;
    double P;
    double T;
    std::vector<double> X;
    std::vector<std::string> specie;
    //Foam::speciesTable* s;
    //Foam::dictionary dict;
    void update();
    //Foam::PengRobinsonM<Foam::specie>* PR;
    solver(std::string);
    void solve();
    double vaporfra;
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
};
*/
%}
%include"VLEinterface.h"
/*
struct solver
{
    std::string path;
    double P;
    double T;
    std::vector<double> X;
    std::vector<std::string> specie;
    //speciesTable* s;
    void update();
    solver(std::string);
    void solve();
    double vaporfra;
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
};
*/