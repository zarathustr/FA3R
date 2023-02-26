#include "FA3R.h"
#include <random>
#include <Eigen/Core>
#include <Eigen/Dense>

Eigen::MatrixXd randomMatrix(int m, int n, std::random_device *rd)
{
    std::mt19937 gen((*rd)());  //here you could also set a seed
    std::uniform_real_distribution<double> dis(-1000.0, 1000.0);
	return Eigen::MatrixXd::NullaryExpr(m,n,[&](){return dis(gen);});;
}

int main()
{
	std::random_device rd;
	const int len = 1000;

	std::vector<Eigen::Vector3d> P, Q;
	for(int i = 0; i < len; ++i)
    {
	    auto P_ = randomMatrix(3, 1, &rd);
        auto Q_ = randomMatrix(3, 1, &rd);
        P.push_back(P_);
        Q.push_back(Q_);
    }

	Eigen::Matrix3d rRes_FA3Rd, rRes_FA3Ri, rRes_eig;
    Eigen::Vector3d tRes_FA3Rd, tRes_FA3Ri, tRes_eig;
    FA3R_double(&P, &Q, nullptr, 20, &rRes_FA3Rd, &tRes_FA3Rd);
    FA3R_int(&P, &Q, nullptr, 20, &rRes_FA3Ri, &tRes_FA3Ri);
    eig3D_eig(&P, &Q, nullptr, &rRes_eig, &tRes_eig);

    std::cout << "FA3R double:" << std::endl;
    std::cout << "R: " << rRes_FA3Rd << std::endl;
    std::cout << "t: " << tRes_FA3Rd << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "FA3R int:" << std::endl;
    std::cout << "R: " << rRes_FA3Ri << std::endl;
    std::cout << "t: " << tRes_FA3Ri << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Eig:" << std::endl;
    std::cout << "R: " << rRes_eig << std::endl;
    std::cout << "t: " << tRes_eig << std::endl;

    return 0;
}
