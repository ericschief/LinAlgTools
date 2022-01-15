#include<iostream>
#include<ctime>
#include<cmath>
#include<vector>
#include<array>
#include<algorithm>
#include<iomanip>
#include<complex>
#include<fstream>


using namespace std;

//these are a number of helpful functions for expediting work
template <class T> void print_my_vector(vector<T> vec) //code for printing vectors is always useful
{
  for(T item:vec)
  {
    cout << item << " ";
  }
  cout << "\n";
}


template <class T> void print_my_matrix(vector<vector<T>> mat)
{
  for(int i = 0; i<mat.size();i++)
  {
    print_my_vector(mat[i]);
  }
  cout << endl;
}



template <class T, int n> void write_to_file(vector<T> x, vector<T> y,string fname)
{

  fstream fout (fname + to_string(n) + ".txt",ios_base::out);

  for(int i = 0; i < x.size(); i++)
  {
    fout << x[i] << ",";
    fout << y[i] << "\n";
  }
  fout.close();
}

template <class T> T AbsError(T approx, T actual)
{
  T val = abs(approx - actual);
  return val;
}

template <class T> T MaxError(vector<T> vec1, vector<T> vec2)
{
  vector<T> errorvec(vec1.size());
  T error = 0;
  for(int i = 0; i < vec1.size(); i++)
  {
    errorvec[i] = abs(vec1[i]-vec2[i]);
  }

  error = *max_element(errorvec.begin(),errorvec.end());
return error;
}

template<class T> vector<T> vec_subtract(vector<T> x,vector<T> y)
{
  vector<T> z(x.size());

  for (int i = 0; i < x.size(); i++)
  {
    z[i] = x[i]-y[i];
  }
  return z;
}

template<class T> vector<T> vec_addition(vector<T> x, vector<T> y)
{
  vector<T> z(x.size());

  for (int i = 0; i < x.size(); i++)
  {
    z[i] = x[i] + y[i];
  }
  return z;
}
template<class T> T vec_multiply(vector<T>x, vector<T>y)
{
  int n = x.size();
  T val = 0;
  for(int i = 0; i < n; i++)
  {
    val += x[i]*y[i];
  }
  return val;
}

template<class T> vector<T> scalar_distrib(vector<T> x, T alpha)
{
  int n = x.size();
  vector<T> z(n);
  for(int i = 0; i < n; i++)
  {
    z[i] = x[i]*alpha;
  }

  return z;
}

template<class T> vector<complex<T>> vec_division(vector<complex<T>> x, vector<complex<T>> y)
{
  vector<complex<T>> z(x.size());
  for(int i = 0; i < x.size(); i++ )
  {
    z[i] = x[i]/y[i];
  }


  return z;
}

// VECTOR NORMS
template <class T> T vector_one_norm(vector<T> vec)
{
  T val = 0;
  for(int i = 0; i < vec.size(); i++)
  {

    val += abs(vec[i]);
  }
  return val;
}

template<class T> T vector_two_norm(vector<T> vec)
{
  T val = 0;

  for(int i = 1; i < vec.size()-1; i++)
  {

    val += pow(vec[i],2);
  }
  T twonorm = sqrt(val);

  //cout << setprecision(8) << "The vector's two-norm is: " << twonorm << endl;
  return twonorm;
}


//MATRIX NORMS
template<class T> T matrix_one_norm(vector<vector<T>> A)
{
  int n = A.size();
  vector<T> sums(n);
  for(int i = 0; i < n ; i++)
  {
    for(int j = 0; j < n; j++ )
    {
      sums[i] += abs(A[j][i]);
    }
  }
  T max = *max_element(sums.begin(), sums.end());
  return max;
}

template<class T> T frobenius_norm(vector<vector<T>> A)
{
  int n = A.size();
  T val = 0;

  for(int i = 0; i < n ; i++)
  {
    for(int j = 0; j < n; j++ )
    {
      val += abs(pow(A[i][j],2));
    }
  }

  val = sqrt(val);
  return val;
}

template<class T> vector<vector<T>> absval(vector<vector<T>> A)
{
  int n = A.size();
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      A[i][j]= abs(A[i][j]);
    }
  }

  return A;
}
