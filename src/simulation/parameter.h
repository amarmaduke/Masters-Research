#ifndef PARAMETER_H
#define PARAMETER_H

struct parameter
{
  int m;
  int n;
  double beta;
  double len;
  double gamma;
  double epsilon;
  double sigma;

  double sub_h;
  double sub_count;
  double * delta;

  parameter()
  {
    m = 1;
    n = 4;
    beta = 1;
    len = 1;
    gamma = 100;
    epsilon = 1;
    sigma = 1;

    sub_h = 1;
    sub_count = 10;
  };

};

struct triple
{
  double * first;
  double * second;
  double * third;

  triple(double *f, double *s, double *t)
  {
    first = f;
    second = s;
    third = t;
  }
};

#endif

