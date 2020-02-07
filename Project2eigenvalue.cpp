#include <stdio.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

void get_augmented(float** augm,float** mat1,float* mat2,int row){  // dynamically allocated matrices used for parameters
	
   for(int i=0;i<row;i++){
   	  for(int j=0;j<row+1;j++){
   	  	if(j<row)
   	  	augm[i][j]=mat1[i][j];   // first n columns of augmented matrix consists of the columns of a 
   	  	else
   	  	augm[i][j]=mat2[i]; 	  // last column of augmmented matrix is taken from the b	
		}
   }
}

void gaussian_elimination(float** augm,int row){
    float rate;
     int holder;
    
    for(int m=0;m<row;m++){
    
       for(int l=m;l<m;l++){
       	if(abs(augm[l][m])>abs(augm[m][m])){  // for partial pivoting magnitudes of rows are compared
       		for(int p=0;p<m;p++){             
       			holder=augm[l][p];            // interchange of rows, if needed
       			augm[l][p]=augm[m][p];        
       			augm[m][p]= holder;
			   }
		   }
	   }
    	
    	for(int k=m;k<row-1;k++){            // loop for gaussian elimination
    		rate=augm[k+1][m]/augm[m][m];    // calculating ratio of the pivot elements and the entries above them
    		for(int t=0;t<row+1;t++){
      		augm[k+1][t]=augm[k+1][t]-rate*augm[m][t];   // making the entries above pivot 0
      	}	 	
	}
	}
}


void makepivots1 (float** aug,int row){

    for(int i=0;i<row;i++){
    	float divider=aug[i][i];
    	for(int j=0;j<row+1;j++){
    		aug[i][j]=aug[i][j]/divider;   // dividing every column to its pivot to make all pivots 1
		}                                  // for the ease of back substitution
	}
}

void seperate_augmented (float** aug,float** a,float* b,int row){

	for(int i=0;i<row;i++){
		for(int j=0;j<row;j++){ 
			a[i][j]=aug[i][j];      // first n columnns of augmented matrix will be new version of matrix a
    }
	}

	for(int i=0;i<row;i++){
		for(int j=row;j<row+1;j++){
			b[i]=aug[i][j];          // last column of augmented will be held by matrix b
	}
}
}

void back_substitution(float* x,float** a,float* b, int row){
									
	for(int k=row-1;k>=0;k--){	
	x[k]=b[k];                     // since pivots are 1, last entry of x matrix directly will be equal to b
	for(int i=k-1;i>=0;i--){
		b[i]=b[i]-x[k]*a[i][k];    // portions of x calculated above(with its appropriate coefficients) is substracted from other matrix b elements
	}	
}
}

int singularity_test(float** aug,int row,float tolerance){    
	float control=0;
    for(int i=0;i<row;i++){
    	
    		if(aug[i][i]<tolerance){
    			control=control+1;
		}
	}
	return control;
}


class matrix {
public:
  float**A; 
  float*X;
  
  float evalue(float**A, float*x,float t,int row){         // normalized power iteration function

	    int cont=1;
	    int counter=0;
	    float previous=0;
		float max;
		float y[row];  
		
		float eigenvalue;
 
    while(cont==1){

    	for(int i=0;i<row;i++){
			y[i]=0;
		}
	
    	counter=counter+1;
		
	    for(int i=0;i<row;i++){
  	  	for(int j=0;j<row;j++){
  	  		y[i]=y[i]+A[i][j]*x[j];           // matrix multiplication A.x
			}
		} 
	/*	  for(int i=0;i<row;i++){
  	  	for(int j=0;j<row;j++){
  	  		cout<<A[i][j]<<" ";           // matrix multiplication A.x
			}
			cout<<endl;
		}
		
		
		cout<<endl;
			
		for(int i=0;i<row;i++){
			cout<<"y"<<y[i]<<endl;
		}
		
		cout<<endl;*/
    
		max=y[0];
		
		for(int i=0;i<row;i++){
			if(y[i]>max){
				max=y[i];                 // max of y is equal to infinity norm
			}
		}
  	 /*    cout<<" y inf:"<<max<<endl;
	    for(int i=0;i<row;i++){
	    x[i]=y[i]/max;
		}
			for(int i=0;i<row;i++){
			cout<<"x"<<x[i]<<endl;
		}
		
				cout<<endl;*/
		
		if((abs(previous-max)<t && abs(previous-max)!=0)|| counter==1000)             // checking the convergence wrt tolerance
		break; 
		else{
			previous=max;
		}	
	}
	
	
	
	eigenvalue=max;
	return eigenvalue;       
  }
  
  float evaluesmall(float**A, float*x,float t,int row){         // inverse power iteration function
        int cont=1;
	    int counter=0;
	    float previous=0;
		float max;
		
		float eigenvalue;
	
    float temp[row][row];
    
    for( int k=0;k<row;k++){
   	for(int j=0;j<row;j++){ 
	    temp[k][j]=A[k][j];                  
	  }
   } 
   
	while(cont==1){
	float y[row];  
    
	
	float **aug=new float*[row];    // defining dynamically allocated multidimensional array for augmented matrix
    for(int i=0;i<row+1;i++){
   		aug[i]=new float[row];
    }
    
   
     get_augmented(aug,A,x,row);  
     gaussian_elimination(aug,row);
     makepivots1(aug,row); 	
     float control=0;
		control=singularity_test(aug,row,t);   // unchanged control variable indicates singulariy

		if (control!=0){
		cout<<"Error! Matrix is singular";
		exit;
     	}
     	
     seperate_augmented(aug,A,x,row);
	 back_substitution(y,A,x,row);
    
    	counter=counter+1;    
		max=y[0];
		
		for(int i=0;i<row;i++){
			if(y[i]>max){
				max=y[i];                 // max of y is equal to infinity norm
			}
		}
  	     
	    for(int i=0;i<row;i++){
	    x[i]=y[i]/max;
		}
		
		if((abs(previous-max)<t && abs(previous-max)!=0) || counter==1000)             // checking the convergence wrt tolerance
		break; 
		else{
			previous=max;
		}	
		
	for( int k=0;k<row;k++){
   	for(int j=0;j<row;j++){ 
	    A[k][j]=temp[k][j];                  // holding the values of txt file in a matrix
	  }
   } 
   
	}
	eigenvalue=1/max;
	return eigenvalue;       
  }
  
  // eigenvector is the last version of x itself
};



int main( int argc, char *argv[] )  {
	 float tolerance;
	if(argc==4){
    tolerance= atof(argv[2]);
    }
    
    else{
    	cout<<"Error";
	}
	

   string line;
   int n=0;
   ifstream myfile (argv[1]);              // opening the chosen txt fil
   if (myfile.is_open())
   {
     while ( getline (myfile,line) )
     {
       n++;
     }
     myfile.close();

   }

   else {
   cout << "Unable to open file";
}
   float **a=new float*[n];                // creating dynamically allocated memory
	for( int i=0;i<n;i++){
   		a[i]=new float[n];
   }

   ifstream second (argv[1]);

   for( int k=0;k<n;k++){
   	for(int j=0;j<n;j++){ 
	    second>>a[k][j];                  // holding the values of txt file in a matrix
	  }
   }

    float*x;                              // arbitrary vector that would be used in normalized iteration algorithm
    x=new float[n];

    for(int i=0;i<n;i++){
    	x[i]=1/sqrt(n);
	}
	
	float tempx[n];
	int compare=0;
	

	   float eival1;
	   matrix matx;
	  
    	 
	 ofstream mytext;               // opening .txt file to write in the solution
     mytext.open (argv[3]);

        mytext<<"eigenvalue#1: "<<matx.evalue(a,x,tolerance,n)<<endl<<"eigenvector: "<<endl;

     for(int i=0;i<n;i++){
    	mytext<<x[i]<<endl;    // writing values of matrix is one by one to x.txt file
	 }
        mytext<<endl;
	 for(int i=0;i<n;i++){
    	x[i]=1/sqrt(n);
	}
	
	float rate;
	rate=x[1]/tempx[1];
	
	for(int i=0;i<n;i++){
		if(x[i]/tempx[i]!=rate){
			compare=compare+1;
		}
	}
	
	cout<<compare;
	if(compare==0){
		x[1]=x[1]+1;
	}
	
	  mytext<<"eigenvalue#2: "<<matx.evaluesmall(a,x,tolerance,n)<<endl<<"eigenvector: "<<endl;
     for(int i=0;i<n;i++){
    	mytext<<x[i]<<endl;    // writing values of matrix is one by one to x.txt file
	 }
      mytext.close();
     
  


   return 0;
}
