#include <assert.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdlib>

using std::cout;
using std::cin;
using std::endl;
using std::min;
using std::ios;

int h, wid, n1, n2 ;
inline int distance(int x, int y, int x1, int y1)
{return (abs(x-x1)+abs(y-y1));}

double scorer(int wid, int h, const int p1[],const int p2[]);
void pointdist(int x, int y, double data[], const int array[], int manip, int size);
void graph(int x, int y, double score, int h, int wid);
void ntoa (int n, int siz, int x, int y, int max, int a[]);
void adder(int place, int max, int siz, int *a);
void converter(int siz, int x, int y, int *a, int *b);
void simplegraph(int xmax, int ymax, const int *p1, const int *p2, int stores);


int main()
{   n1 = 3; n2 = 3;h=4;wid=4; 
    const int p1[8] = {1,1,4,2,2,3};
    int p2[8];//={0,0,2,1,1,2};
    //scorer(2,2,p1,p2);
    
    
    for (int i =0; i < 100*23; i++)
        {/*printf("%d\n", i)*/; ntoa(i, n2,wid+1,h+1,(h+1)*(wid+1), p2); scorer(wid,h,p1,p2);}
    return 0;
}


double scorer(int wid, int h, const int p1[],const int p2[])
{   double score1, score2, eff;
    score1=0; score2=0; eff=0;
    double mx= (wid+1)*(h+1)+1;
    for (int cnth = h; cnth >= 0; cnth--)
	{   for (int cntw = 0; cntw <= wid; cntw++)
		{   double data[]={0,mx,0,mx,0}; // eff, array, distance, mult

		    pointdist(cntw,cnth,data,p1,0,2*(n1)); 
		    pointdist(cntw,cnth,data,p2,2,2*(n2)); 
		    data[0]=min(data[1], data[3]);
		    eff=data[0]+eff;

		for(int cnt1=1; cnt1 < 4; cnt1= cnt1+2)
			{if (data[cnt1] > data[0])
			    {data[cnt1+1]=0;}
			} 
		assert((data[2]+data[4]));
		score1 = data[2] / (data[2]+data[4])+score1;
		//cout <<data[2]+data[4]<<","<< score1 << std::endl;
		score2 = data[4] / (data[2]+data[4])+score2;
		//graph(cntw,cnth, data[2] / (data[2]+data[4]), h, wid);
		}	
	}
    if (score1 < score2)
    {   printf("\n");
	    cout << "Player 1: " << score1 << endl;
	    cout << "Player 2: " << score2 << endl;
	    for (int i =0; i < ( 2 * n2 ); i++)
		    {cout<< p2[i];}
		printf("\n");
		simplegraph(wid,h,p1,p2,n1);
    }

    return eff;
}

void pointdist(int x, int y, double data[], const int array[], int manip, int size)
{   int mindist= 1 + (wid+1) * (1+h);
    for (int acnt=0; acnt < size; acnt = acnt+2)
	{   if (distance(array[acnt], array[acnt+1], x, y) < mindist)
		{mindist=distance(array[acnt], array[acnt+1], x, y); 
		data[manip+1] = mindist;
		data[manip+2] = 1; 
		}
	else if (distance(array[acnt], array[acnt+1], x, y) == mindist)
		{data[manip+2]=data[manip+2]+1;}
	else{continue;}
	}
}

void graph(int x, int y, double score1, int h, int w)
{	
    cout.setf(ios::showpoint);
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);
    cout << 2 - score1; 
	
    if (x < w)
		{cout << " - ";}

	
    if (x == w && y)
	{   cout << endl;
	    {for (int xcnt=0; xcnt <= w; xcnt++)
			{cout << " |    ";}
		cout << endl;
		}
	}
	if (x == w && !y) {printf("\n");}
} 

void ntoa (int n, int siz, int x, int y, int max, int a[])
{   int place=siz-1;
    int change= max-siz;
    int * temp = (int*) calloc(x*y, sizeof(int));
    for(int i = 0; i < x*y; i++)
        {temp[i] = i+1;}
    while (n--> 0)
        {adder(place, max, siz, temp);}
    
    converter(siz, x, y , a, temp);
    free(temp);
}

void adder(int place, int max, int siz, int *a)
{   
    if (place < 0) {printf("Error"); exit(1);}

    if (a[place] < max - (siz -place - 1))
        { a[place]= a[place]+1;
          for (int i = place+1; i< siz; i++)
            {a[i] = a[place] + i-place;}
        }
    else {adder(place-1, max,siz,a);}
}

void converter(int siz, int x, int y, int a[], int *b)
{
    for (int i =0; i < siz; i++)
        {   a[2*i +1] = (b[i] / x); 
            if (b[i] % x)
                { a[2*i] = (b[i] % x) - 1;}
            else 
                {
                  a[2*i+1] = a[2*i+1] - 1;
                  a[2*i] = x-1;      
                }
        }
}

void simplegraph(int xmax, int ymax, const int *p1,const int *p2, int stores)
{   int store1=0; int store2 =0;
    for(int ycount = ymax ; ycount>=0 ; --ycount)
    {   for (int xcount = 0; xcount <= xmax; xcount++)
        {   for (int i = 0; i < 2*stores; i= i+2)
            {   if (p1[i] == xcount && p1[i+1] == ycount) {store1=store1+1;}
                if (p2[i] == xcount && p2[i+1] == ycount) {store2=store2+1;}
            }
            if (store1 && store2) {printf(" X ");store1=0; store2=0;}
            else if (store1) {printf (" 1 ");store1=0;}
            else if (store2) {printf (" 2 ");store2=0;}
            else {printf(" O ");}
            if (xcount == xmax)
                {   printf("\n");
                    if (ycount)
                        for (int j =0; j <= ymax ;j++)
                            printf(" |    ");
                    printf("\n");
                }
            else printf(" - ");
        }
    }
}
