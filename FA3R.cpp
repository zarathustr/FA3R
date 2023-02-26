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

    double max = 0;
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            auto tmp = fabs((*sigma_)(i, j));
            if(tmp > max)
                max = tmp;
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

    double max = 0;
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            auto tmp = fabs((*sigma_)(i, j));
            if(tmp > max)
                max = tmp;
        }
    *sigma_ /= max;
    
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


void eig3D_eig(const vector<Vector3d>* P,
	const vector<Vector3d>* Q,
	Matrix3d * sigma,
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
	
    Matrix3d A = (*sigma_) - sigma_->transpose();
    Matrix3d tmp;

    Vector3d D(A(1, 2), A(2, 0), A(0, 1));
    Matrix4d QQ;
    QQ(0, 0) = (*sigma_)(0, 0) + (*sigma_)(1, 1) + (*sigma_)(2, 2);
    tmp = (*sigma_) + sigma_->transpose();
    tmp(0, 0) -= QQ(0, 0);    tmp(1, 1) -= QQ(0, 0);    tmp(2, 2) -= QQ(0, 0);
    QQ(0, 1) = D.x();     QQ(0, 2) = D.y();     QQ(0, 3) = D.z();
    QQ(1, 0) = D.x();     QQ(2, 0) = D.y();     QQ(3, 0) = D.z();

    QQ(1, 1) = tmp(0, 0); QQ(1, 2) = tmp(0, 1); QQ(1, 3) = tmp(0, 2);
    QQ(2, 1) = tmp(1, 0); QQ(2, 2) = tmp(1, 1); QQ(2, 3) = tmp(1, 2); 
    QQ(3, 1) = tmp(2, 0); QQ(3, 2) = tmp(2, 1); QQ(3, 3) = tmp(2, 2);

    SelfAdjointEigenSolver<Matrix4d> es(QQ);

    double max_eig = 0.0;
    int max_index = 0;
    for(int i = 0; i < 4; ++i)
        if(max_eig < es.eigenvalues()(i, 0))
            max_index = i;
	
    Quaterniond q(es.eigenvectors().col(max_index));
    q.normalize();

    Matrix3d R = q.toRotationMatrix();
    (*rRes)(0, 0) = - R(2, 2);  (*rRes)(0, 1) =   R(2, 1);  (*rRes)(0, 2) = - R(2, 0);
    (*rRes)(1, 0) =   R(1, 2);  (*rRes)(1, 1) = - R(1, 1);  (*rRes)(1, 2) =   R(1, 0);
    (*rRes)(2, 0) =   R(0, 2);  (*rRes)(2, 1) = - R(0, 1);  (*rRes)(2, 2) =   R(0, 0);
    *tRes = mean_X - (*rRes).transpose() * mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
        delete sigma_;
}


void eig3D_symbolic(const vector<Vector3d>* P,
	const vector<Vector3d>* Q,
	Matrix3d * sigma,
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
	
    Matrix3d A = (*sigma_) - sigma_->transpose();
    Matrix3d tmp;
    Vector3d D(A(1, 2), A(2, 0), A(0, 1));
    Matrix4d QQ;
    QQ(0, 0) = (*sigma_)(0, 0) + (*sigma_)(1, 1) + (*sigma_)(2, 2);
    tmp = (*sigma_) + sigma_->transpose();
    tmp(0, 0) -= QQ(0, 0);    tmp(1, 1) -= QQ(0, 0);    tmp(2, 2) -= QQ(0, 0);
    QQ(0, 1) = D.x();     QQ(0, 2) = D.y();     QQ(0, 3) = D.z();
    QQ(1, 0) = D.x();     QQ(2, 0) = D.y();     QQ(3, 0) = D.z();

    QQ(1, 1) = tmp(0, 0); QQ(1, 2) = tmp(0, 1); QQ(1, 3) = tmp(0, 2);
    QQ(2, 1) = tmp(1, 0); QQ(2, 2) = tmp(1, 1); QQ(2, 3) = tmp(1, 2);
    QQ(3, 1) = tmp(2, 0); QQ(3, 2) = tmp(2, 1); QQ(3, 3) = tmp(2, 2);

    double c = QQ.determinant();
    double b = - 8.0 * sigma_->determinant();
    double a = - 2.0 * ((*sigma_)(0, 0) * (*sigma_)(0, 0) + (*sigma_)(0, 1) * (*sigma_)(0, 1) + (*sigma_)(0, 2) * (*sigma_)(0, 2) + 
                        (*sigma_)(1, 0) * (*sigma_)(1, 0) + (*sigma_)(1, 1) * (*sigma_)(1, 1) + (*sigma_)(1, 2) * (*sigma_)(1, 2) + 
                        (*sigma_)(2, 0) * (*sigma_)(2, 0) + (*sigma_)(2, 1) * (*sigma_)(2, 1) + (*sigma_)(2, 2) * (*sigma_)(2, 2));

    double T0 = 2.0 * a * a * a + 27.0 * b * b - 72.0 * a * c;
    double tt = a * a + 12.0 * c;
    double theta = atan2(sqrt(4.0 * tt * tt * tt - T0 * T0), T0);
    double aT1 = 1.259921049894873 * sqrt(tt) * cos(theta * 0.333333333333333333);
    double T2 = sqrt( - 4.0 * a + 3.174802103936399 * aT1);
    double lambda = 0.204124145231932 * (T2 + sqrt( - T2 * T2 - 12.0 * a - 29.393876913398135 * b / T2));
	
    double G11 = QQ(0, 0) - lambda, G12 = QQ(0, 1), G13 = QQ(0, 2), G14 = QQ(0, 3);		
    double G22 = QQ(1, 1) - lambda, G23 = QQ(1, 2), G24 = QQ(1, 3);
    double G33 = QQ(2, 2) - lambda, G34 = QQ(2, 3);
    double G44 = QQ(3, 3);

    Quaterniond qRes = Quaterniond(
		       G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34,
			   G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34,
	           G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34, 
	           - (G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33));
    qRes.normalize();

    *rRes = qRes.toRotationMatrix();
    *tRes = mean_X - (*rRes).transpose() * mean_Y;

    if(P != nullptr && Q != nullptr && sigma == nullptr)
        delete sigma_;
}
