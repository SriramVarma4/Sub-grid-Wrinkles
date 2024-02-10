#include <GL/glut.h>
#include <Eigen/Core>
#include <vector>


#ifndef TIME_INTEGRATOR_H

#define TIME_INTEGRATOR_H

class Integrator {
protected:
    std::vector<Eigen::VectorXf> q;
    // std::vector<Eigen::VectorXf> dq_k;
    std::vector<Eigen::VectorXf> p;

public:
    Integrator(const std::vector<Eigen::VectorXf>& Q, /*const std::vector<Eigen::VectorXf>& dq,*/ const std::vector<Eigen::VectorXf>& P): q(Q), /*dq_k(dq),*/ p(P) {
        
    }

    Integrator(std::vector<Eigen::VectorXf>&& Q, std::vector<Eigen::VectorXf>&& P): q(std::move(Q)), /*dq_k(dq),*/ p(std::move(P)) {
        
    }

    // decltype(auto) getq(int i, int j){
    //     return q[i][j];
    // }

    // decltype(auto) getdq(int i, int j){
    //     return dq_k[i][j];
    // }

    // decltype(auto) getp(int i, int j){
    //     return p[i][j];
    // }

    virtual void update() = 0;

    virtual ~Integrator() {}

    // virtual float L_evaluate(const Eigen::VectorXf& q, const Eigen::VectorXf& dq) = 0;

    // virtual Eigen::VectorXf dldq(const Eigen::VectorXf& q, const Eigen::VectorXf& dq) = 0;

    // virtual Eigen::VectorXf dldq_2(const Eigen::VectorXf& q, const Eigen::VectorXf& dq) = 0;

    // virtual Eigen::VectorXf solveD1Ld(const Eigen::VectorXf& q, const Eigen::VectorXf& p) = 0;

    // virtual Eigen::VectorXf solveD2Ld(const Eigen::VectorXf& q, const Eigen::VectorXf& fq) = 0;

    virtual void run() = 0;
};



#endif