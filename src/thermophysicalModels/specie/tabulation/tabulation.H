#ifndef TABULATION2_H
#define TABULATION2_H
#include"table.H"
#include"RTree.H"

//Hongyuan:if there are overlapping in data, inerpolation will be irreversible in these region
struct coord_int
{
	int m[3];
};

class tabulation:public simple_tabulation,public RTree<coord_int, double,3,double>
{
	public:
		struct Tetr;
		struct Hex;
		struct Tri_2d;
		struct Quad_2d;
		tabulation(const char*filename,int va=3);
		int Search_(const double a[]);//Search point (T,p,yi),return value is useless, the position of value can be getten from ret_coord
		void getBox(double min[],double max[],Point arr_point[]);//calculate upper and lower bound of point set
		friend bool SearchCallback(coord_int  id, void* arg);//used to search box
		double interpolation_tab(double x,double y,double z,int index);//reverse interpolation,index=1(T),2(p),only pure interpolation, for extrapolation return 1e308
		static inline double mixed_product(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC,const Point& a_pointD);//mixed product
		static inline double abs(const double& a);//absolute value 
		static inline double cross_2d(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC);//cross product
		static bool point_in_tab_quad_2d(const Point& a_point, const Quad_2d& a_quad_2d);//this and following functions judge whether a point is in a shape 
		bool point_in_tab_quad_2d(const Point& a_point, const Quad_2d& a_quad_2d,coord_int  id);
		static bool point_in_tab_hex(const Point& a_point, const Hex& a_hex);
		bool point_in_tab_hex(const Point& a_point, const Hex& a_hex,coord_int  id);
		static bool point_in_tab_tri_2d(const Point& a_point, const Tri_2d& a_tri_2d);
		struct Tetr
		{
			Tetr() {};
			Tetr(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC,const Point& a_pointD)
			{
				vertex[0]=a_pointA;
				vertex[1]=a_pointB;
				vertex[2]=a_pointC;
				vertex[3]=a_pointD;
			};
			Point vertex[4];
		};

		struct Hex
		{
			Point vertex[8];
		};
		struct Tri_2d
		{
			Point vertex[3];
			Tri_2d(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC)
			{
				vertex[0]=a_pointA;
				vertex[1]=a_pointB;
				vertex[2]=a_pointC;
			};
		} ;
		struct Quad_2d
		{
			Point vertex[4];
		} ;
		static bool point_in_hex(const Point& a_point, const Hex& a_hex);
		static	bool point_in_tetr(const Point& a_point,const Tetr &a_tetr);

//	private:
		Point s_point;
		coord_int ret_coord;

} ;
const double eps = 1e-7;

bool SearchCallback(coord_int  id, void* arg);

#ifdef NoRepository
    #include "tabulation.C"
#endif


#endif