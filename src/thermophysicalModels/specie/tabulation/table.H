
#ifndef TABULATION_H
#define TABULATION_H

#include<fstream>
#include<iostream>
#include<sstream>
#include<cstring>
#include<algorithm>

/*

class tabulation by Hongyuan,
computing interpolation from a table, this version can only handle three independent variable.
sample:
int main() {

	tabulation a("ThermoTable.txt"); //"ThermoTable.txt":the path of table 
	double input[3]= {105,2,0.999};//input argument
	std::cout<<a.interpolation(input,1)<<std::endl;//1 is the index of dependent variable, the index starts from 1.
	return 0;
}
*/
struct point_find {
	int prev,next;
};
class simple_tabulation {
	public:
	//protected:
		//rownum: the number of rows (ignoring the first line)
		//colnum: the number of colnums,depvnum:the number of input colnum,idndepvnum:the number of output colnum
		
		int rownum,colnum,depvnum,indepvnum;  
		int dnum[3];// the number of diferent value in colnum 0,1,2
		char filename_[1000];//input file path 
		double indv[3][50000];//diferent value in colnum 0,1,2. The size can be changed 
		double **ttable,****table;//ttable: temporary table, table:table to store output values, input values are stored in indv
		//find_diff find different value in a[],then store in dv, a_num:num of a[], dv_num:num of dv
		// used to get indv[][]
		void find_diff(double a[],int a_num,double dv[],int & dv_num);
		//calculate rownum,and colnum 
		void initalfile(const char*);
		//read table
		void read(const char*);
		//reshape the input table to build  ****table
		void reshape();
   //use binary search to find value "in" in a[l],a[l+1],...,a[r-1],a[r],
// reture a point_find, prev is the index of lower bound, next is index of upper bound
//no lower bound  : prev=-1
//no upper bound :next=-1
   point_find find(double in, const double a[],int l,int r)const; 
	public:
		simple_tabulation(const char*filename,int va=3);
		//get tabulation
		double 	interpolation(double a[],int index)const ;
		~ simple_tabulation();
};

//find value "b" in a[], a[] has max_ variables 
int find_reshape(double *a,double b,int max_) ;



#ifdef NoRepository
    #include "table.C"
#endif



#endif

