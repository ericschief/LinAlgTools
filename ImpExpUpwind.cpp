#include<iostream>
#include<cmath>
#include<vector>
#include<array>
#include<fstream>
#include"tools.h"
using namespace std;
const double pi = 4*arctan(1);



// initial and boundary conditions
template <class T> T initial_cond1(T x)
{
  return sin(2*pi*x);
}
template <class T> T initial_cond2(T x)
{
    return exp(-pow(x-0.5,2)/(.02));
}
template<class T> T initial_cond3(T x)
{
  return cos(2*pi*x);
}
template<class T> T initial_cond4(T x)
{
  return 0.002*exp(7*x);
}
template<class T> T initial_cond5(T x)
{
  if( x < 0.5)
  {
    return 0;
  }
  return 1;
}


template<class T> T boundary_cond1(T x)
{
  return x/2.0+0.5;
}
template<class T> T boundary_cond2(T x)
{
  return sin(4*pi*x);
}
template<class T> T boundary_cond3(T x)
{
  return sin(x);
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

// takes a mesh and an initial condition function and puts out the proper vector
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

template<class T> vector<T> exp_up_wind(vector<T> x_mesh, vector<T> time_steps,vector<T> v_prev,T a, T (*f)(T))
{
  int n = time_steps.size();
  T h = x_mesh[1]-x_mesh[0];
  int m = x_mesh.size();
  T k = time_steps[1]-time_steps[0];
  v_prev[0] = f(x_mesh[0]);
  vector<T> v_cur(m);
  T val = a*k/h;



  for (int i = 1; i < n; i++)
  {
    v_cur[0] = f(time_steps[i]);    //assign boundary condition
    for(int j = 1; j < m; j++)
    {
      v_cur[j] = v_prev[j]-val*(v_prev[j]-v_prev[j-1]); //explicit method
    }
    v_prev = v_cur;
  }
  return v_cur;
}



//solvers used for the system encountered in implicit method
template<class T> void Bi_Diag_solve1(T c_0, T c_1,vector<T> v_prev, vector<T> &v_cur)
{
  int m = v_prev.size();

  for(int i = 1; i < m; i++)
  {
    v_cur[i] = (v_prev[i-1]-c_0*v_cur[i-1])/-c_1;
  }

}

template<class T> void Bi_Diag_solve2(T c_0, T c_1,vector<T> v_prev, vector<T> &v_cur)
{
  int m = v_prev.size();

  for(int i = 1; i < m; i++)
  {
    v_cur[i] = (v_prev[i]+c_0*v_cur[i-1])/c_1;
  }

}



template<class T> vector<T> imp_up_wind(vector<T> x_mesh, vector<T> time_steps,vector<T> v_prev,T a, T (*f)(T))
{
  int n = time_steps.size();
  T h = x_mesh[1]-x_mesh[0];
  int m = x_mesh.size();
  T k = time_steps[1]-time_steps[0];
  v_prev[0] = f(x_mesh[0]);
  vector<T> v_cur(m);

  T val1 = a*k/h;
  T val2 = 1+val1;

  for (int i = 1; i < n; i++)
  {
      v_cur[0] = f(time_steps[i]);
      Bi_Diag_solve2<T>(val1,val2,v_prev,v_cur); //solve the system at every time step
      v_prev = v_cur;
  }

  return v_cur;
}

template <class T> vector<T> dirichlet_solve(vector<T> x,T t, T a, T (*f)(T), T (*g)(T))
{
  int m = x.size();
  vector<T>vec(m);
  T val = -1.0/a;

  for(int i = 0; i < m; i++)
  {

    if((x[i]-a*t) < 0)
    {
      vec[i] = g(val*(x[i]-a*t));
      continue;
    }
    vec[i] = f(x[i]-a*t);
  }


  return vec;
}


int main()
{


  return 0;
}
