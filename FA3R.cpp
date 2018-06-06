#include "FA3R.h"
#include <stdint.h>
#include <time.h>  

const uint8_t shift = 20;
const double d2l = pow(2, (double) shift);
const int64_t shifted1 = ((int64_t) 2) << shift;
const int64_t shifted2 = shift << 1;
const int64_t shifted3 = ((int64_t) 2) << (((int64_t) 3) * shift + 1);
const double l2d = 1.0 / d2l;


void cross(const int64_t &x1, const int64_t &x2, const int64_t &x3, 
           const int64_t &y1, const int64_t &y2, const int64_t &y3, 
           const int64_t &k, int64_t *z1, int64_t *z2, int64_t *z3)
{
    *z1 = (k * (*z1 + ((x2 * y3 - x3 * y2) >> shift))) >> shifted2;
    *z2 = (k * (*z2 + ((x3 * y1 - x1 * y3) >> shift))) >> shifted2;
    *z3 = (k * (*z3 + ((x1 * y2 - x2 * y1) >> shift))) >> shifted2;
}

void cross(const Vector3d &x, const Vector3d &y, const double &k, Vector3d &z)
{
    z(0) = k * (z(0) + x(1) * y(2) - x(2) * y(1));
    z(1) = k * (z(1) + x(2) * y(0) - x(0) * y(2));
    z(2) = k * (z(2) + x(0) * y(1) - x(1) * y(0));
}

void FA3R_int(const vector<Vector3d>* P,
              const vector<Vector3d>* Q,
              Matrix3d * sigma,
              int num,
              Matrix3d * rRes,
              Vector3d * tRes)
{
    int64_t hx1, hx2, hx3,
         hy1, hy2, hy3,
         hz1, hz2, hz3;
    int64_t hx1_, hx2_, hx3_,
         hy1_, hy2_, hy3_,
         hz1_, hz2_, hz3_;

    Matrix3d * sigma_ = sigma;
    Vector3d mean_X, mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
    {
       sigma_ = new Matrix3d();
    
       int n = P->size();
       mean_X.setZero();
       mean_Y.setZero();

       for (int i = 0; i < n; ++i)
       {
          mean_X = mean_X + (*P)[i];
          mean_Y = mean_Y + (*Q)[i];
       }
       mean_X = mean_X / n;
       mean_Y = mean_Y / n;

       sigma_->setZero();

       for (int i = 0; i < n; ++i)
       {
          *sigma_ = *sigma_ + ((*Q)[i] - mean_Y) * (((*P)[i] - mean_X).transpose());
       }
       *sigma_ = *sigma_ / n;
    }

    hx1 = (int64_t)((*sigma_)(0, 0) / max * d2l);  hx2 = (int64_t)((*sigma_)(0, 1) / max * d2l);  hx3 = (int64_t)((*sigma_)(0, 2) / max * d2l);
    hy1 = (int64_t)((*sigma_)(1, 0) / max * d2l);  hy2 = (int64_t)((*sigma_)(1, 1) / max * d2l);  hy3 = (int64_t)((*sigma_)(1, 2) / max * d2l);
    hz1 = (int64_t)((*sigma_)(2, 0) / max * d2l);  hz2 = (int64_t)((*sigma_)(2, 1) / max * d2l);  hz3 = (int64_t)((*sigma_)(2, 2) / max * d2l);

    for(int i = 0; i < num; ++i)
    {
        hx1_ = hx1; hx2_ = hx2; hx3_ = hx3;
        hy1_ = hy1; hy2_ = hy2; hy3_ = hy3; 
        hz1_ = hz1; hz2_ = hz2; hz3_ = hz3;

        int64_t k = shifted3   / (((hx1_ * hx1_ + hx2_ * hx2_ + hx3_ * hx3_ +
                                     hy1_ * hy1_ + hy2_ * hy2_ + hy3_ * hy3_ +
                                     hz1_ * hz1_ + hz2_ * hz2_ + hz3_ * hz3_) >> shift) + shifted1);

        cross(hx1_, hx2_, hx3_, hy1_, hy2_, hy3_, k, &hz1, &hz2, &hz3);
        cross(hz1_, hz2_, hz3_, hx1_, hx2_, hx3_, k, &hy1, &hy2, &hy3);
        cross(hy1_, hy2_, hy3_, hz1_, hz2_, hz3_, k, &hx1, &hx2, &hx3);
    }

    Vector3d Hx(((double) hx1) * l2d, ((double) hy1) * l2d, ((double) hz1) * l2d), 
             Hy(((double) hx2) * l2d, ((double) hy2) * l2d, ((double) hz2) * l2d), 
             Hz(((double) hx3) * l2d, ((double) hy3) * l2d, ((double) hz3) * l2d);
    Hx.normalize();
    Hy.normalize();
    Hz.normalize();

    (*rRes)(0, 0) = Hx(0);  (*rRes)(0, 1) = Hy(0);  (*rRes)(0, 2) = Hz(0);
    (*rRes)(1, 0) = Hx(1);  (*rRes)(1, 1) = Hy(1);  (*rRes)(1, 2) = Hz(1);
    (*rRes)(2, 0) = Hx(2);  (*rRes)(2, 1) = Hy(2);  (*rRes)(2, 2) = Hz(2);
    *tRes = mean_X - (*rRes).transpose() * mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
        delete sigma_;
}


void FA3R_double(const vector<Vector3d>* P,
	             const vector<Vector3d>* Q,
                 Matrix3d * sigma,
                 int num,
		         Matrix3d * rRes,
		         Vector3d * tRes)
{
    Matrix3d * sigma_ = sigma;
    Vector3d mean_X, mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
    {
       sigma_ = new Matrix3d();
    
       int n = P->size();
       mean_X.setZero();
       mean_Y.setZero();

       for (int i = 0; i < n; ++i)
       {
          mean_X = mean_X + (*P)[i];
          mean_Y = mean_Y + (*Q)[i];
       }
       mean_X = mean_X / n;
       mean_Y = mean_Y / n;

       sigma_->setZero();

       for (int i = 0; i < n; ++i)
       {
          *sigma_ = *sigma_ + ((*Q)[i] - mean_Y) * (((*P)[i] - mean_X).transpose());
       }
       *sigma_ = *sigma_ / n;
    }
    
    Vector3d hx((*sigma_)(0, 0), (*sigma_)(1, 0), (*sigma_)(2, 0));
    Vector3d hy((*sigma_)(0, 1), (*sigma_)(1, 1), (*sigma_)(2, 1));
    Vector3d hz((*sigma_)(0, 2), (*sigma_)(1, 2), (*sigma_)(2, 2));
    Vector3d hx_, hy_, hz_;
    double k;
    
    for(int i = 0; i < num; ++i)
    {
        k = 2.0 / (hx(0) * hx(0) + hx(1) * hx(1) + hx(2) * hx(2) +
                   hy(0) * hy(0) + hy(1) * hy(1) + hy(2) * hy(2) +
                   hz(0) * hz(0) + hz(1) * hz(1) + hz(2) * hz(2) + 1.0);
        
        hx_ = hx;  hy_ = hy; hz_ = hz;
        
        cross(hx_, hy_, k, hz);
        cross(hz_, hx_, k, hy);
        cross(hy_, hz_, k, hx);
    }

    (*rRes)(0, 0) = hx(0);  (*rRes)(0, 1) = hy(0);  (*rRes)(0, 2) = hz(0);
    (*rRes)(1, 0) = hx(1);  (*rRes)(1, 1) = hy(1);  (*rRes)(1, 2) = hz(1);
    (*rRes)(2, 0) = hx(2);  (*rRes)(2, 1) = hy(2);  (*rRes)(2, 2) = hz(2);
	*tRes = mean_X - (*rRes).transpose() * mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
        delete sigma_;
}
