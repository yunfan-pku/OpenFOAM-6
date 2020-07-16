#include"tabulation.H"


void tabulation::getBox(double min[],double max[],Point arr_point[])
{
	for(int iDim=0; iDim<3; iDim++)
	{
		min[iDim]=arr_point[0].m[iDim];
		max[iDim]=arr_point[0].m[iDim];
		for(int iP=1; iP<8; iP++)
		{
			if(min[iDim]>arr_point[iP].m[iDim])
				min[iDim]=arr_point[iP].m[iDim];
			if(max[iDim]<arr_point[iP].m[iDim])
				max[iDim]=arr_point[iP].m[iDim];
		}
	}
}

//((b-a)x(c-a)).(d-a)
inline double tabulation::mixed_product(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC,const Point& a_pointD)
{
	return (a_pointB.m[0]-a_pointA.m[0])*(a_pointC.m[1]-a_pointA.m[1])*(a_pointD.m[2]-a_pointA.m[2])\
	       +(a_pointB.m[1]-a_pointA.m[1])*(a_pointC.m[2]-a_pointA.m[2])*(a_pointD.m[0]-a_pointA.m[0])\
	       +(a_pointB.m[2]-a_pointA.m[2])*(a_pointC.m[0]-a_pointA.m[0])*(a_pointD.m[1]-a_pointA.m[1])\
	       -(a_pointB.m[1]-a_pointA.m[1])*(a_pointC.m[0]-a_pointA.m[0])*(a_pointD.m[2]-a_pointA.m[2])\
	       -(a_pointB.m[2]-a_pointA.m[2])*(a_pointC.m[1]-a_pointA.m[1])*(a_pointD.m[0]-a_pointA.m[0])\
	       -(a_pointB.m[0]-a_pointA.m[0])*(a_pointC.m[2]-a_pointA.m[2])*(a_pointD.m[1]-a_pointA.m[1]);

}

inline double tabulation::abs(const double& a)
{
	return a>0?a:-a;
}

bool tabulation::point_in_tetr(const Point& a_point,const Tetr &a_tetr)
{
	double sum=0,vol=0,t;
	t=mixed_product(a_tetr.vertex[0],a_tetr.vertex[1],a_tetr.vertex[2],a_tetr.vertex[3]);
	vol=abs(t);
	t=mixed_product(a_point,a_tetr.vertex[1],a_tetr.vertex[2],a_tetr.vertex[3]);
	sum+=abs(t);
	t=mixed_product(a_tetr.vertex[0],a_point,a_tetr.vertex[2],a_tetr.vertex[3]);
	sum+=abs(t);
	t=mixed_product(a_tetr.vertex[0],a_tetr.vertex[1],a_point,a_tetr.vertex[3]);
	sum+=abs(t);
	t=mixed_product(a_tetr.vertex[0],a_tetr.vertex[1],a_tetr.vertex[2],a_point);
	sum+=abs(t);
	if(sum-eps>vol)
		return false;
	return true;
}
inline double tabulation::cross_2d(const Point& a_pointA,const Point& a_pointB,const Point& a_pointC)
{
	return (a_pointB.m[0]-a_pointA.m[0])* (a_pointC.m[1]-a_pointA.m[1])-(a_pointB.m[1]-a_pointA.m[1])* (a_pointC.m[0]-a_pointA.m[0]);
}

bool tabulation::point_in_tab_tri_2d(const Point& a_point, const Tri_2d& a_tri_2d)
{
	double sum=0,area=abs(cross_2d(a_tri_2d.vertex[0],a_tri_2d.vertex[1],a_tri_2d.vertex[2])),t;
	t=cross_2d(a_point,a_tri_2d.vertex[1],a_tri_2d.vertex[2]);
	sum+=abs(t);
	t=cross_2d(a_tri_2d.vertex[0],a_point,a_tri_2d.vertex[2]);
	sum+=abs(t);
	t=cross_2d(a_tri_2d.vertex[0],a_tri_2d.vertex[1],a_point);
	sum+=abs(t);
	if(sum-eps>area)
		return false;
	return true;
}

bool tabulation::point_in_tab_quad_2d(const Point& a_point, const Quad_2d& a_quad_2d)
{
	if(!point_in_tab_tri_2d( a_quad_2d.vertex[2], Tri_2d(a_quad_2d.vertex[0],a_quad_2d.vertex[1],a_quad_2d.vertex[3])))
	{

		if(point_in_tab_tri_2d( a_point, Tri_2d(a_quad_2d.vertex[0],a_quad_2d.vertex[1],a_quad_2d.vertex[3])))
			return true;
		if(point_in_tab_tri_2d( a_point, Tri_2d(a_quad_2d.vertex[1],a_quad_2d.vertex[2],a_quad_2d.vertex[3])))
			return true;
	}
	else
	{
		if(point_in_tab_tri_2d( a_point, Tri_2d(a_quad_2d.vertex[0],a_quad_2d.vertex[2],a_quad_2d.vertex[3])))
			return true;
		if(point_in_tab_tri_2d( a_point, Tri_2d(a_quad_2d.vertex[0],a_quad_2d.vertex[1],a_quad_2d.vertex[2])))
			return true;
	}
	return false;
}


bool tabulation::point_in_tab_quad_2d(const Point& a_point, const Quad_2d& a_quad_2d,coord_int  id)
{
	//double ip[3],ti[3];
	long double x1,x2,y1,y2,v11,v12,v21,v22,w11,w12,w21,w22,vol;
	x1=indv[0][id.m[0]];
	x2=indv[0][id.m[0]+1];
	y1=indv[1][id.m[1]];
	y2=indv[1][id.m[1]+1];
	vol=(x2-x1)*(y2-y1);
	//std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
	//ti[2]=z;
	//ti[0]=x1;
	//ti[1]=y1;
	v11=a_quad_2d.vertex[0].m[0]/vol;
	w11=a_quad_2d.vertex[0].m[1]/vol;
	//ti[0]=x1;
	//ti[1]=y2;
	v12=a_quad_2d.vertex[3].m[0]/vol;
	w12=a_quad_2d.vertex[3].m[1]/vol;
	//ti[0]=x2;
	//ti[1]=y1;
	v21=a_quad_2d.vertex[1].m[0]/vol;
	w21=a_quad_2d.vertex[1].m[1]/vol;
	//ti[0]=x2;
	//ti[1]=y2;
	v22=a_quad_2d.vertex[2].m[0]/vol;
	w22=a_quad_2d.vertex[2].m[1]/vol;

	long double A1,A2,A3,A4,B1,B2,B3,B4;
	//A1xy+A2x+A3y+A4=0
	//interpolation(ip,int index)
	A1=v11-v12-v21+v22;
	B1=w11-w12-w21+w22;
	A2=y1*(v12-v22)+y2*(v21-v11);
	B2=y1*(w12-w22)+y2*(w21-w11);
	A3=x1*(v21-v22)+x2*(v12-v11);
	B3=x1*(w21-w22)+x2*(w12-w11);
	A4=v22*x1*y1-v12*x2*y1-v21*x1*y2+v11*x2*y2-a_point.m[0];
	B4=w22*x1*y1-w12*x2*y1-w21*x1*y2+w11*x2*y2-a_point.m[1];
	long double a1,a2,a3;
	a1=A2*B1-A1*B2;
	a2=A4*B1- A3*B2+A2*B3-A1*B4;
	a3=A4*B3-A3*B4;
	long double x,y,delta=a2*a2-4*a1*a3;
	if(delta<0)
		return false;

	if(fabs(a1)<(10e-12))
	{
		x=-a3/a2;
		if(x>=x1-eps&&x<=x2+eps)
		{
			y=-(A4+A2*x)/(A3+A1*x);
			if(y>=y1-eps&&y<=y2+eps)
				return true;
		}
	}


	x=(-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);
	if(x<x1-eps||x>x2+eps)
		x=(-a2+sqrt(a2*a2-4*a1*a3))/(2*a1);
	if(x<x1-eps||x>x2+eps)
		return false;

	y=-(A4+A2*x)/(A3+A1*x);
	if((y-y1)/y1<-eps||(y-y2)/y2>eps)
		return false;
	return true;

}

bool tabulation::point_in_tab_hex(const Point& a_point, const Hex& a_hex)
{
	double z1=a_hex.vertex[0].m[2],z2=a_hex.vertex[4].m[2],zp=a_point.m[2];
	double r1=(z2-zp)/(z2-z1),r2=(zp-z1)/(z2-z1);
	Quad_2d quad;
	for(int i=0; i<4; i++)
		quad.vertex[i]=r1*a_hex.vertex[i]+r2*a_hex.vertex[4+i];
	return point_in_tab_quad_2d(a_point, quad);
}

bool tabulation::point_in_tab_hex(const Point& a_point, const Hex& a_hex,coord_int  id)
{
	double z1=a_hex.vertex[0].m[2],z2=a_hex.vertex[4].m[2],zp=a_point.m[2];
	double r1=(z2-zp)/(z2-z1),r2=(zp-z1)/(z2-z1);
	Quad_2d quad;
	for(int i=0; i<4; i++)
		quad.vertex[i]=r1*a_hex.vertex[i]+r2*a_hex.vertex[4+i];
	return point_in_tab_quad_2d(a_point, quad,id);
}

bool tabulation::point_in_hex(const Point& a_point, const Hex& a_hex)
{
	if(point_in_tetr(a_point,Tetr(a_hex.vertex[0],a_hex.vertex[1],a_hex.vertex[2],a_hex.vertex[4])))
		return true;
	if(point_in_tetr(a_point,Tetr(a_hex.vertex[1],a_hex.vertex[4],a_hex.vertex[5],a_hex.vertex[7])))
		return true;
	if(point_in_tetr(a_point,Tetr(a_hex.vertex[1],a_hex.vertex[2],a_hex.vertex[3],a_hex.vertex[7])))
		return true;
	if(point_in_tetr(a_point,Tetr(a_hex.vertex[2],a_hex.vertex[4],a_hex.vertex[6],a_hex.vertex[7])))
		return true;
	if(point_in_tetr(a_point,Tetr(a_hex.vertex[1],a_hex.vertex[2],a_hex.vertex[4],a_hex.vertex[7])))
		return true;
	return false;

}
/*
                 7------6
			    /��    /��
			  4------5  ��
              �� 3---��-2
              ��/    ��/
              0 -----1


*/

tabulation::tabulation(const char*filename,int va):simple_tabulation(filename, va)
{
	Point t_Point[8];
	double min[3],max[3];
	int a[8][3]= {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
	int i[3];//ti[3];
	coord_int id;
	for( i[0]=0; i[0]<dnum[0]-1; i[0]++)
		for( i[1]=0; i[1]<dnum[1]-1; i[1]++)
			for( i[2]=0; i[2]<dnum[2]-1; i[2]++)
			{

				for(int ia=0; ia<8; ia++)
				{
					for(int id=0; id<2; id++)
					{

						t_Point[ia].m[id]=table[i[0]+a[ia][0]][i[1]+a[ia][1]][i[2]+a[ia][2]][id];
					}
					t_Point[ia].m[2]=indv[2][i[2]+a[ia][2]];
				}
				id.m[0]=i[0];
				id.m[1]=i[1];
				id.m[2]=i[2];

				getBox(min,max,t_Point);
				//if(id==0)
				//	std::cout<< min[0]<<" "<<min[1]<<" "<<min[2]<<" "<<max[0]<<" "<<max[1]<<" "<<max[2]<<std::endl;
				int t=0;
				if(id.m[0]==500&&id.m[1]==14&&id.m[2]==10)
					t++;
				Insert(min,max,id);// zhanshishanchu!!!!!!!!!!!!!!!!!!!!
				//id++;

			}
} ;


bool SearchCallback(coord_int  id, void* arg)
{
	typename tabulation::Hex hex;
	tabulation* This=reinterpret_cast<tabulation*>(arg);
	int a[8][3]= {0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
	for(int i=0; i<8; i++)
	{
		hex.vertex[i].m[0]=This->table[id.m[0]+a[i][0]][id.m[1]+a[i][1]][id.m[2]+a[i][2]][0];
		hex.vertex[i].m[1]=This->table[id.m[0]+a[i][0]][id.m[1]+a[i][1]][id.m[2]+a[i][2]][1];
		hex.vertex[i].m[2]=This->indv[2][id.m[2]+a[i][2]];
	}

	if(This->point_in_tab_hex(This->s_point, hex,id))
	{
		This->ret_coord=id;

		return false;
	}
	return true; // keep going
}

int tabulation::Search_(const double a[])
{
	s_point.m[0]=a[0];
	s_point.m[1]=a[1];
	s_point.m[2]=a[2];
	ret_coord.m[0]=-1;
	ret_coord.m[1]=-1;
	ret_coord.m[2]=-1;
	return 	Search(a, reinterpret_cast<void*>(this), SearchCallback);
}

double tabulation::interpolation_tab(double ix,double iy,double z,int index)
{
	double ip[3],ti[3];
	long double x1,x2,y1,y2,v11,v12,v21,v22,w11,w12,w21,w22,vol;
	ip[0]=ix;
	ip[1]=iy;
	ip[2]=z;
	Search_(ip);
	x1=indv[0][ret_coord.m[0]];
	x2=indv[0][ret_coord.m[0]+1];
	y1=indv[1][ret_coord.m[1]];
	y2=indv[1][ret_coord.m[1]+1];
	vol=(x2-x1)*(y2-y1);
	//std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
	ti[2]=z;
	ti[0]=x1;
	ti[1]=y1;
	v11=interpolation(ti,1)/vol;
	w11=interpolation(ti,2)/vol;
	ti[0]=x1;
	ti[1]=y2;
	v12=interpolation(ti,1)/vol;
	w12=interpolation(ti,2)/vol;
	ti[0]=x2;
	ti[1]=y1;
	v21=interpolation(ti,1)/vol;
	w21=interpolation(ti,2)/vol;
	ti[0]=x2;
	ti[1]=y2;
	v22=interpolation(ti,1)/vol;
	w22=interpolation(ti,2)/vol;

	long double A1,A2,A3,A4,B1,B2,B3,B4;
	//A1xy+A2x+A3y+A4=0
	//interpolation(ip,int index)
	A1=v11-v12-v21+v22;
	B1=w11-w12-w21+w22;
	A2=y1*(v12-v22)+y2*(v21-v11);
	B2=y1*(w12-w22)+y2*(w21-w11);
	A3=x1*(v21-v22)+x2*(v12-v11);
	B3=x1*(w21-w22)+x2*(w12-w11);
	A4=v22*x1*y1-v12*x2*y1-v21*x1*y2+v11*x2*y2-ix;
	B4=w22*x1*y1-w12*x2*y1-w21*x1*y2+w11*x2*y2-iy;
	long double a1,a2,a3;
	a1=A2*B1-A1*B2;
	a2=A4*B1- A3*B2+A2*B3-A1*B4;
	a3=A4*B3-A3*B4;
	double x,y;
	if(abs(a1)<(10e-12))
	{
		x=-a3/a2;
		if(x>=x1-eps&&x<=x2+eps)
		{
			y=-(A4+A2*x)/(A3+A1*x);
			if((y>=y1-eps&&y<=y2+eps)||a1==0)
			{
				if(index==1)
					return x;
				return y;
			}
		}
	}
	if(a2*a2-4*a1*a3>=0)
	{

		x=(-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);
		if(x<x1-eps||x>x2+eps)
			x=(-a2+sqrt(a2*a2-4*a1*a3))/(2*a1);
	}
	else
	{
		x=(-a2)/(2*a1);
	}
	
	y=-(A4+A2*x)/(A3+A1*x);
	//std::cout<<x<<" "<<y<<std::endl;
	if(index==1)
		return x;
	return y;



}

