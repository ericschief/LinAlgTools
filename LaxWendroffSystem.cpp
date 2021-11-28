#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include"tools.h"
using namespace std;
//pde info
const double pi = 4.0*atan(1);
const double a = 0.7;
const double x_start = 0.0;
const double x_end  = 1.0;
const double TT = 1.0;


//need to ensure that stability conditions are met for u_0+c and u_0c
const double u_0 = 0.25;
const double c = 0.30;
const double rho = 1.5;
//matrices needed to simplify computation
const vector<double> Lambda{u_0+c,u_0-c};
const vector<vector<double>>X{{rho*c,-rho*c},{1,1}};
const vector<vector<double>>Xinverse {{1.0/(2.0*rho*c),0.5},{-1.0/(2.0*rho*c),0.5}};





template<class T> T initial_cond1(T x)
{
  return x/2.0;
}
template<class T> T initial_cond2(T x)
{
  return cos(x);
}
template<class T> T initial_cond3(T x)
{
  return x*x;

}
template<class T> T initial_cond4(T x)
{
  return exp(-(pow(x-0.5,2))/.02);
}

template<class T> T initial_cond5(T x)
{
  return cos(2.0*pi*x);
}

template<class T> T boundary_cond1(T x)
{
  return x/2.0;
}
template<class T> T boundary_cond2(T x)
{
  return exp(x);
}



//makes a mesh
template<class T>  vector<T> mesh(int n, T a, T b)
{
  vector<T> vec(n+1);
  T sum  = a;
  T h = (b-a)/n;
  for(int i = 0; i < vec.size(); i++)
  {
    vec[i] = a + h*i;
  }
  return vec;
}
//makes vector of initial condition function on mesh
template<class T> vector<T> init_cond(vector<T> mesh, T (*f)(T))
{
  int n = mesh.size();
  vector<T>v_0(n);
  for (int i = 0; i < n; i++)
  {
    v_0[i] = f(mesh[i]);
  }
  return v_0;
}

template<class T> vector<vector<T>> lax_wendroff_sys(int M, T(*f1)(T),T(*f2)(T))
{
  vector<T> x_mesh = mesh<T>(M,x_start,x_end);
  T h  = (x_end-x_start)/M;

  const int N  = M;

  int bb  = (int)ceil(-Lambda[1]*TT/h)+1;

  vector<T> time = mesh<T>(N,0,TT);
  T k = TT/N;


  //initial conditions
  vector<vector<T>> v_prev {{init_cond<T>(x_mesh,f1)},init_cond<T>(x_mesh,f2)};

  //multiply by Xinverse to make the equations easily solvable
  vector<vector<T>> w_prev(2,vector<T>(M+1));
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < M+1; j++)
    {

      for(int k = 0; k < 2; k++)
      {
        w_prev[i][j] += Xinverse[i][k]*v_prev[k][j];
      }
    }
  }

  T val1 = pow(k,2)/(2.0*pow(h,2));
  T val2 = k/(2.0*h);



  vector<vector<T>> w_cur(2,vector<T>(M+1));
  for(int i = 1; i < N+1; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        for(int k = 1; k < M; k++)
        {

          w_cur[j][k] = w_prev[j][k]+(val1*pow(Lambda[j],2))*(w_prev[j][k+1]-2*w_prev[j][k]+w_prev[j][k-1])-(val2*Lambda[j])*(w_prev[j][k+1]-w_prev[j][k-1]);
        }
      }
      //Now apply numerical boundary conditions before starting next time step
      //left boundary
      w_cur[1][0] = w_prev[1][1];
      w_cur[0][0] =  1 - w_cur[1][0];
      //right boundary
      w_cur[0][M] = w_prev[0][M-1];
      w_cur[1][M] = w_cur[0][M]-(1.0/(rho*c));
      w_prev = w_cur;
    }

    vector<vector<T>> soln(2,vector<T>(M+1));
    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < M+1; j++)
      {

        for(int k = 0; k < 2; k++)
        {
          soln[i][j] += X[i][k]*w_cur[k][j];
        }
      }
    }



  return soln;
}

template<class T>vector<vector<T>> analytical_soln(int M,T(*f1)(T),T(*f2)(T))
{
  vector<T> x_mesh = mesh<T>(M,x_start,x_end);
  T h  = (x_end-x_start)/M;
  const int N  = (int)ceil(Lambda[0]*TT/h)+1;
  vector<T> time = mesh<T>(N,0,TT);
  T k = TT/N;

  vector<vector<T>> soln(2,vector<T>(M+1));

  T tHat{0};

  for(int i =0; i < 1200 ; i++)
  {
    tHat = (1.0/Lambda[0])*(i*h-Lambda[0]*TT);
    soln[0][i] = 1.0-.5*(((-1.0)/(rho*c))*f1(-Lambda[1]*tHat)+f2(-Lambda[1]*tHat));
    soln[1][i] = 0.5*(((-1.0)/(rho*c))*f1(i*h-Lambda[1]*TT)+f2(i*h-Lambda[1]*TT));
  }
  for(int i = 1200;i < 1900; i++)
  {
    soln[0][i]=0.5*(((1.0)/(rho*c))*f1(i*h-Lambda[0]*TT)+f2(i*h-Lambda[0]*TT));
    soln[1][i]=0.5*(((-1.0)/(rho*c))*f1(i*h-Lambda[1]*TT)+f2(i*h-Lambda[1]*TT));
  }

  for(int i = 1900; i < M+1 ; i++)
  {
    tHat = (1.0/Lambda[1])*(1.0-(i*h-Lambda[1]*TT));
    soln[0][i] = 0.5*(((1.0)/(rho*c))*f1(i*h-Lambda[0]*TT)+f2(i*h-Lambda[0]*TT));
    soln[1][i] = 0.5*(((1.0)/(rho*c))*f1(1.0-Lambda[0]*tHat)+f2(1.0-Lambda[0]*tHat))-(1.0/(rho*c));

  }


  vector<vector<T>> soln2(2,vector<T>(M+1));
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < M+1; j++)
    {

      for(int k = 0; k < 2; k++)
      {
        soln2[i][j] += X[i][k]*soln[k][j];
      }
    }
  }



  return soln2;
}

template<class T> vector<T> Lax_Wendroff(int M, T(*f)(T),T(*g)(T))
{
  vector<T> x_mesh = mesh<T>(M,x_start,x_end);
  T h  = (x_end-x_start)/M;
  const int N  = (int)ceil(a*TT/h)+1;
  vector<T> time = mesh<T>(N,0,TT);
  T k = TT/N;

  vector<T>v_prev = init_cond<T>(x_mesh,f);
  vector<T>v_cur(M+1);
  T val1 = a*k/(2*h);
  T val2 = pow(a,2)*pow(k,2)/(2.0*pow(h,2));

  for(int i = 1; i <= N; i++)
    {
      v_cur[0] = g(time[i]);
      for(int j = 1; j < M; j++)
      {
        v_cur[j] = v_prev[j]-val1*(v_prev[j+1]-v_prev[j-1])+val2*(v_prev[j+1]-2.0*v_prev[j]+v_prev[j-1]);
      }
      v_cur[M] = v_prev[M-1];
      v_prev = v_cur;

    }
  return v_cur;
}
template<class T> vector<T> ONED_analytic(int M,T(*f)(T), T(*g)(T))
{
  vector<T> soln(M+1);
  vector<T>x= mesh<T>(M,x_start,x_end);
  T val = -1.0/a;
  for(int i = 0; i < M+1; i++)
  {
    if((x[i]-a*TT)>=0)
    {
      soln[i] = f(x[i]-a*TT);
    }
    else
    {
      soln[i] = g(val*(x[i]-a*TT));
    }
  }
  return soln;
}


int main()
{
  //choose space mesh size, time mesh is created to ensure stability
  const int M = 2000;
  vector<double> x = mesh<double>(M,x_start,x_end);

  vector<vector<double>> numerical = lax_wendroff_sys<double>(M,initial_cond1,initial_cond2);

  vector<vector<double>> analytic = analytical_soln<double>(M,initial_cond1,initial_cond2);

  //write results of analytical and numerical to files
  write_to_file<double,0>(x,analytic[0],"analytic");
  write_to_file<double,1>(x,analytic[1],"analytic");

  write_to_file<double,0>(x,numerical[0],"numerical");
  write_to_file<double,1>(x,numerical[1],"numerical");







  return 0;
}
