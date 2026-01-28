#include "roots.hpp"
#include <cmath>
#include <iostream>

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    if (!root) return false; // null pointer check

    if (a > b) { double t = a; a = b; b = t; } //  if b < a, swap

    const double tol = 1e-6;
    const int max_iters = 1000000;

    double fa = f(a);
    double fb = f(b);

    // endpoint checks
    if (fabs(fa) <= tol) { *root = a; return true; }
    if (fabs(fb) <= tol) { *root = b; return true; }

    // need opposite signs
    if (fa * fb > 0) return false;

    for (int it = 0; it < max_iters; ++it) 
    {
        double c = (a + b) / 2.0;
        double fc = f(c);

        if (fabs(fc) <= tol) {
            *root = c;
            return true;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }

        if ((b - a) / 2.0 <= tol) { //if the interval is sufficiently small
            *root = (a + b) / 2.0;
            return true;
        }
    }

    return false;
}


bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
        if (!root) return false; // null pointer check

        if (a > b) { double t = a; a = b; b = t; } // if b < a, swap

    const double tol = 1e-6;
    const int max_iters = 1000000;

    double fa = f(a);
    double fb = f(b);

    // endpoint checks
    if (fabs(fa) <= tol) { *root = a; return true; }
    if (fabs(fb) <= tol) { *root = b; return true; }

    // need opposite signs
    if (fa * fb > 0) return false;

    for (int it = 0; it < max_iters; ++it) 
    {
        double denom = fb - fa;
        if (fabs(denom) < 1e-15) return false; // avoid division by zero
        
        double c = a - (fa * (b - a)) / denom; // Regula Falsi formula
        if (c < a || c > b) return false; // ensure c is within [a,b]
        double fc = f(c);

        if (fabs(fc) <= tol) {
            *root = c;
            return true;
        }

        if (fa * fc < 0) { // root is in [a,c] 
            b = c;
            fb = fc;
        } else { // root is in [c,b]
            a = c;
            fa = fc;
        }

    }

    return false;
}


bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    if (!root) return false;

    if (a > b) { double t = a; a = b; b = t; } // if b < a, swap

    if (c < a || c > b) return false;

    const double tol = 1e-6;
    const int max_iters = 1000000;

    double x = c;

    for (int it = 0; it < max_iters; ++it)
    {
        double fx = f(x);

        if (fabs(fx) <= tol) {
            *root = x;
            return true;
        }

        double gx = g(x);

        // derivative too small 
        if (fabs(gx) < 1e-15)
            return false;

        double x_next = x - fx / gx;

        // check if next guess is within [a,b]
        if (x_next < a || x_next > b)
            return false;

        // checks if it is close enough to root
        if (fabs(x_next - x) <= tol) {
            *root = x_next;
            return true;
        }

        x = x_next; // define x for the next iteration
    }

    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    if (!root) return false;

    if (a > b) { double t = a; a = b; b = t; } // if b < a, swap

    if (c < a || c > b) {
        

        return false;
    }

    const double tol = 1e-6;
    const int max_iters = 1000000;

    double x0 = a;
    double x1 = c;

    for (int it = 0; it < max_iters; ++it)
    {
        double f0 = f(x0);
        double f1 = f(x1);

        double denom = f1 - f0;
        if (fabs(denom) < 1e-15){

            return false;
        } // avoid division by zero

        double x2 = x1 - ((f1 * (x1 - x0)) / denom);

        
        double f2 = f(x2);
            if (fabs(f2) <= tol) {
            *root = x2;
            return true;
        }

        // checks if it is close enough to root
        if (fabs(x2 - x1) <= tol) {
            *root = x2;
            return true;
        }

        x0 = x1;
        x1 = x2; // define for next iteration
    }

    return false;
}