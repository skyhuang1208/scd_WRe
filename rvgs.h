#pragma once
/* -------------------------------------------------------------
 * Name            : rvgs.h (header file for the library rvgs.c)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 11-03-96
 * -------------------------------------------------------------- 
 */

#if !defined( _RVGS_ )
#define _RVGS_

int Bernoulli(double p);
int Binomial(int n, double p);
int Equilikely(int a, int b);
int Geometric(double p);
int Pascal(int n, double p);
int Poisson(double m);

double Uniform(double a, double b);
double Exponential(double m);
double Erlang(int n, double b);
double Normal(double m, double s);
double Lognormal(double a, double b);
double Chisquare(int n);
double Student(int n);

#endif

