
#include"table.H" 



simple_tabulation::simple_tabulation(const char*filename,int va):rownum(0),colnum(0),indepvnum(va) {
	strcpy(filename_,filename);
	initalfile(filename_);
	read(filename_);
	depvnum=colnum-va;
	reshape();
}
void simple_tabulation::initalfile(const char* filename) {
	char s[10000];
	std::ifstream fin(filename);

	while(!fin.eof()) {
		fin.getline(s,10000);
		if(strlen(s)==0)
			break;
		rownum++;
		if(rownum==2) {
			std::stringstream line(s);
			double t;
			while(!line.eof()) {
				line>>t;
				colnum++;
			}
		}
	}
	fin.close();
	rownum--;
	return ; 

}
void simple_tabulation::read(const char* filename)
{
	//TimerInit;

	FILE * fp;
	fp = fopen (filename, "r");
	//double te;
	ttable =new double*[colnum];
	for (int i=0; i<colnum; i++)
		ttable[i]=new double[rownum];

	//std::ifstream fin(filename);

	char s[10000];
	//fin.getline(s,10000);
  if(fgets(s,10000,fp)!=NULL)
	for (int i=0; i<rownum; i++)
	{

		//fin.getline(s,10000);
		//// std::stringstream line(s);
		for(int j=0; j<colnum; j++)
		{
			//TimerStart;
			//fin>>te;

			if(fscanf(fp, "%lf", &ttable[j][i])==EOF)
      break;
			//TimerEnd;
			//ttable[j][i]=te;

		}
	}


	//TimerShow;
	fclose(fp);
	//fin.close();
	return ;
}
int find_reshape(double *a,double b,int max_) {
	for(int i =0; i<max_; i++)
		if (a[i]==b)
			return i;
      return -1;
}
void simple_tabulation::find_diff(double a[],int a_num,double dv[],int & dv_num) {
	double *tempa=new double [a_num];
	for(int i=0; i<a_num; i++)
		tempa[i]=a[i];

	std::sort(tempa,tempa+a_num);

	dv_num=1;
	dv[0]=tempa[0];
	for (int i=1; i<a_num; i++) {
		if(	tempa[i]!=	tempa[i-1]) {
			dv[dv_num]=tempa[i];
			dv_num++;
		}
	}
	delete [] tempa;
	return ; 

}
void simple_tabulation::reshape() {
	for(int i=0; i<3; i++) {
		find_diff(ttable[i],rownum,	indv[i],dnum[i]);
	}

	table= new double***[dnum[0]];
	for(int i=0; i<dnum[0]; i++) {
		table[i]=new double**[dnum[1]];
		for(int j=0; j<dnum[1]; j++) {
			table[i][j]=new double*[dnum[2]];
			for(int k=0; k<dnum[2]; k++) {
				table[i][j][k]=new double[depvnum];
			}
		}
	}

	int x[3];
	for (int i=0; i<rownum; i++) {
		for (int j=0; j<3; j++) {
			x[j]=find_reshape(indv[j],ttable[j][i],dnum[j]);
		}

		for(int k=0; k<depvnum; k++)
			table[x[0]][x[1]][x[2]][k]=ttable[3+k][i];
	}

	for (int i=0; i<colnum; i++)
		delete[] ttable[i];
	delete[] ttable;
	return ; 
}


point_find simple_tabulation::find(double in,const double a[],int l,int r)const {
	point_find re;
	if (in<a[l]) {
		re.prev=-1;
		re.next=0;
		return re;
	}
	if (in>a[r]) {
		re.next=-1;
		re.prev=r;
		return re;
	}
	if(l==r) {
		re.prev=l;
		re.next=r;
		return re;

	}
	if(l+1==r) {
		if(in==a[l])
			r=l;
		if(in==a[r])
			l=r;
		re.prev=l;
		re.next=r;
		return re;
	}
	int m=(l+r)/2;

	if (in<=a[m]) {
		return find(in,a,l,m);
	} else {
		return find(in,a,m,r);
	}
}

double simple_tabulation::interpolation(double a[],int index) const{
	point_find px[3];

	for (int i=0; i<3; i++)
		px[i]=find(a[i],indv[i],0,dnum[i]-1);

	int x[2][3];
	for (int i=0; i<3; i++) {

		if (px[i].prev==-1) {
			x[0][i]=px[i].next;
			x[1][i]=x[0][i]+1;
		} else if(px[i].next==-1) {
			x[1][i]=px[i].prev;
			x[0][i]=x[1][i]-1;
		} else {
			x[0][i]=px[i].prev;
			x[1][i]=px[i].next;
		}
	}
	double fx[3][2];
	int ar[8][3]= {0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,0,1,0,1,1,1,1,1};

	for(int i=0; i<3; i++) {
		if(x[1][i]==x[0][i]) {
			fx[i][1]=0;
			fx[i][0]=1;
		} else {
			fx[i][1]=(a[i]-indv[i][x[0][i]])/(indv[i][x[1][i]]-indv[i][x[0][i]]);
			fx[i][0]=1-fx[i][1];
		}
	}
	double sum=0;
	for(int i=0; i<8; i++)
		sum+=table[x[ar[i][0]][0]][x[ar[i][1]][1]][x[ar[i][2]][2]][-1+index]*fx[0][ar[i][0]]*fx[1][ar[i][1]]*fx[2][ar[i][2]];

	return sum;
}

simple_tabulation::~simple_tabulation() {

	for(int i=0; i<dnum[0]; i++) {

		for(int j=0; j<dnum[1]; j++) {

			for(int k=0; k<dnum[2]; k++) {
				delete[] table[i][j][k];
			}
			delete[] table[i][j];
		}
		delete[] table[i];
	}
	delete[] table;

}


